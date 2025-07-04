#### Complete Gait Analysis 
#### Subject 2
#### Author: Michael Leighton
#### Date: 31 Dec 2024

library(tidyverse)
library(ggpubr)
library(patchwork)
library(zoo)
library(here)
library(viridis)

#-------------------------------------------------------------------------------
#### SUBJECT-SPECIFIC DETAILS
#-------------------------------------------------------------------------------

#### IMPORTANT!: change these to relevant subject number prior to running code:
# 1. Line ~30 - In dir.create("gait_analysis/plots/SX...)
# 2. Line ~34 - pattern = "SX..."
# 3. Line ~398 - 'ggsave()' near end of analyze_temporal_parameters() function
# 4. Line ~522 - 'ggsave()' on kinematic plots
# 5. Line ~715 - 'ggsave()' on temporal plots
# 6. Line ~556 - save subject ID for kinematic data
# 7. Line ~560 - save subject kinematic data to CSV
# 6. Line ~731 - save subject ID for temporal data
# 7. Line ~735 - save subject temporal data to CSV

# Create directory for complete analysis results
dir.create("gait_analysis/complete_analysis", recursive = TRUE, showWarnings = FALSE)

# Create subject-specific output directory
dir.create("gait_analysis/plots/S2", recursive = TRUE, showWarnings = FALSE)

# Load all processed data for the subject
files <- list.files("gait_analysis/processed_data", 
                    pattern = "S2_.*_processed.rds", 
                    full.names = TRUE)

results_list <- lapply(files, readRDS)

# Extract marker_data from results_list and add speed & condition
marker_data <- bind_rows(
  lapply(results_list, function(result) {
    result$marker_data %>%
      mutate(speed = result$metadata$speed,  # Add speed column
             condition = result$metadata$condition)  # Add condition column
  })
)

#-------------------------------------------------------------------------------
#### KINEMATIC ANALYSIS FUNCTIONS
#-------------------------------------------------------------------------------
analyze_gait_kinematics <- function(results_list) {
  # Process each trial
  kinematic_data <- lapply(results_list, function(result) {
    # Get ankle marker data
    ankle_data <- result$marker_data %>%
      filter(marker_num == 2) %>%
      arrange(frame)
    
    # Calculate total stride length
    total_stride <- max(ankle_data$x) - min(ankle_data$x)
    
    # Create frame-by-frame data with normalized positions
    frame_data <- ankle_data %>%
      group_by(frame) %>%
      summarise(
        ankle_height = y,
        # Normalize x positions to start at 0
        x_pos = x - min(ankle_data$x),  # This shifts everything to start at 0
        stride_length = total_stride,
        .groups = 'drop'
      ) %>%
      # Add knee and ankle angle calculations
      left_join(
        result$marker_data %>%
          group_by(frame) %>%
          summarise(
            knee_angle = if(all(c(2,3,4) %in% marker_num)) {
              hip_x <- x[marker_num == 4]
              hip_y <- y[marker_num == 4]
              knee_x <- x[marker_num == 3]
              knee_y <- y[marker_num == 3]
              ankle_x <- x[marker_num == 2]
              ankle_y <- y[marker_num == 2]
              
              thigh_vector <- c(knee_x - hip_x, knee_y - hip_y)
              shank_vector <- c(ankle_x - knee_x, ankle_y - knee_y)
              
              dot_product <- sum(thigh_vector * shank_vector)
              magnitudes <- sqrt(sum(thigh_vector^2)) * sqrt(sum(shank_vector^2))
              
              # Add safety checks
              if(magnitudes > 0) {
                cos_angle <- dot_product / magnitudes
                # Clamp the value to valid range for acos
                cos_angle <- pmin(pmax(cos_angle, -1), 1)
                180 - (acos(cos_angle) * (180/pi))
              } else {
                NA_real_
              }
            } else {
              NA_real_
            },
            # Add ankle angle calculation
            ankle_angle = if(all(c(1,2,3) %in% marker_num)) {
              knee_x <- x[marker_num == 3]
              knee_y <- y[marker_num == 3]
              ankle_x <- x[marker_num == 2]
              ankle_y <- y[marker_num == 2]
              foot_x <- x[marker_num == 1]
              foot_y <- y[marker_num == 1]
              
              shank_vector <- c(ankle_x - knee_x, ankle_y - knee_y)
              foot_vector <- c(foot_x - ankle_x, foot_y - ankle_y)
              
              dot_product <- sum(shank_vector * foot_vector)
              magnitudes <- sqrt(sum(shank_vector^2)) * sqrt(sum(foot_vector^2))
              
              acos(dot_product / magnitudes) * (180/pi)
            } else {
              NA_real_
            },
            .groups = 'drop'
          ),
        by = "frame"
      ) %>%
      mutate(
        speed = result$metadata$speed,
        condition = result$metadata$condition
      )
    
    return(frame_data)
  })
  
  # Combine all trials
  bind_rows(kinematic_data) %>%
    mutate(
      stride_length = as.numeric(stride_length),
      x_pos = as.numeric(x_pos),
      knee_angle = as.numeric(knee_angle),
      ankle_angle = as.numeric(ankle_angle),
      ankle_height = as.numeric(ankle_height),
      speed = as.numeric(speed),
      condition = as.factor(condition)
    )
}

#-------------------------------------------------------------------------------
#### TEMPORAL ANALYSIS FUNCTIONS
#-------------------------------------------------------------------------------
analyze_temporal_parameters <- function(results_list) {
  # Create lookup tables for thresholds
  threshold_lookup <- tribble(
    ~condition, ~speed, ~stability, ~heel_height, ~velocity, ~toe_prev_vel, ~toe_curr_vel,
    "limp",     2.5,    2.2,       0.6,          -0.04,     0.15,         0.10,
    "normal",   2.5,    2.4,       0.4,          -0.03,     0.22,         0.07,
    "PD",       2.5,    2.0,       0.3,          -0.06,     0.20,         0.10,
    "limp",     4.0,    1.3,       0.8,          -0.06,     0.20,         0.10,
    "normal",   4.0,    1.2,       0.9,          -0.07,     0.18,         0.12,
    "PD",       4.0,    1.5,       0.8,          -0.07,     0.18,         0.12,
    "limp",     5.5,    1.0,       0.8,          -0.08,     0.12,         0.15,
    "normal",   5.5,    1.2,       0.6,          -0.06,     0.15,         0.11,
    "PD",       5.5,    1.0,       0.8,          -0.08,     0.12,         0.15
  )
  
  temporal_data <- lapply(results_list, function(result) {
    # Get ankle and foot data separately
    ankle_data <- result$marker_data %>%
      filter(marker_num == 2) %>%
      arrange(frame)
    
    foot_data <- result$marker_data %>%
      filter(marker_num == 1) %>%
      arrange(frame)
    
    # Combine foot and ankle data
    stance_data <- ankle_data %>%
      select(frame, ankle_x = x, ankle_y = y) %>%
      left_join(
        foot_data %>% select(frame, foot_x = x, foot_y = y),
        by = "frame"
      ) %>%
      # Add metadata first
      mutate(
        speed = result$metadata$speed,
        condition = result$metadata$condition
      ) %>%
      # Join with threshold lookup table
      mutate(
        # Get the index for the matching condition and speed
        lookup_idx = match(
          paste(condition, speed),
          paste(threshold_lookup$condition, threshold_lookup$speed)
        ),
        # Get thresholds using the index
        stability_threshold = threshold_lookup$stability[lookup_idx],
        heel_height_threshold = threshold_lookup$heel_height[lookup_idx],
        velocity_threshold = threshold_lookup$velocity[lookup_idx],
        toe_prev_threshold = threshold_lookup$toe_prev_vel[lookup_idx],
        toe_curr_threshold = threshold_lookup$toe_curr_vel[lookup_idx]
      ) %>%
      select(-lookup_idx) %>%  # Remove the temporary index column
      # Add other computations
      mutate(
        # Adjust window size based on speed
        stability_window = case_when(
          speed <= 3.0 ~ 5,    # Longer window for slow speeds
          speed <= 4.5 ~ 4,    # Medium window
          TRUE ~ 3             # Short window for fast speeds
        ),
        
        # Minimal smoothing
        ankle_smooth = rollapply(ankle_y, width = 3, mean, fill = NA, align = "center"),
        foot_smooth = rollapply(foot_y, width = 3, mean, fill = NA, align = "center"),
        
        # Calculate velocities
        ankle_y_vel = -c(NA, diff(ankle_smooth)),
        foot_y_vel = -c(NA, diff(foot_smooth)),

        # Check for stable ankle position
        ankle_stable = {
          curr_window <- stability_window[1]
          abs(ankle_y - lag(ankle_y, 1, default = first(ankle_y))) < stability_threshold &
            abs(ankle_y - lag(ankle_y, 2, default = first(ankle_y))) < stability_threshold &
            abs(ankle_y - lag(ankle_y, curr_window, default = first(ankle_y))) < stability_threshold
        },
        
        # Detect heel strike using lookup thresholds
        heel_strike = !is.na(ankle_y) &
          ankle_stable &
          ankle_y > mean(ankle_y) + heel_height_threshold &
          !lag(ankle_stable, 1, default = FALSE) &
          lag(ankle_y_vel, 2, default = 0) < velocity_threshold
      ) %>%
      # Add toe-off detection using lookup thresholds
      mutate(
        # First define min_hs_to_to_frames
        min_hs_to_to_frames = case_when(
          condition == "PD" & speed >= 5.5 ~ 8,    # Short for fast PD
          condition == "PD" & speed >= 4.0 ~ 10,   # Medium for medium PD
          condition == "PD" ~ 12,                  # Longer for slow PD
          speed >= 5.5 ~ 10,                       # Short for fast normal/limp
          speed >= 4.0 ~ 12,                       # Medium for medium normal/limp
          TRUE ~ 15                                # Longer for slow normal/limp
        ),
        
        # Then use it in potential_toe_off
        potential_toe_off = !is.na(foot_y_vel) &
          abs(lag(foot_y_vel, 1, default = 0)) < toe_prev_threshold &
          foot_y_vel > toe_curr_threshold &
          lead(foot_y_vel, 1, default = 0) > toe_curr_threshold &
          lead(foot_y_vel, 2, default = 0) > toe_curr_threshold &
          # Add minimum frames since last heel strike check
          (frame - lag(cummax(if_else(heel_strike, frame, 0)))) > min_hs_to_to_frames,
        
        # Then define toe_off
        toe_off = potential_toe_off &
          !lag(potential_toe_off, 1, default = FALSE) &
          !lag(potential_toe_off, 2, default = FALSE) &
          !lag(potential_toe_off, 3, default = FALSE)
      )
    
    # Initialize in_stance to FALSE for all frames
    stance_data$in_stance <- FALSE
    cat("DEBUG: Initialized in_stance\n")  # Add this line
    
    # First, identify heel strikes and toe-offs
    heel_strikes <- which(stance_data$heel_strike == TRUE)
    cat("DEBUG: Created heel_strikes\n")   # Add this line
    toe_offs <- which(stance_data$toe_off == TRUE)
    cat("DEBUG: Created toe_offs\n")       # Add this line
    
    # Initialize first_toe_off globally to ensure it exists for the state machine
    first_toe_off <- if(length(toe_offs) > 0) toe_offs[1] else nrow(stance_data)
    
    # Print debug info
    cat(sprintf("\nProcessing %s at %.1f km/h\n", 
                result$metadata$condition, 
                result$metadata$speed))
    cat("Detected heel strikes at frames:", paste(heel_strikes, collapse = ", "), "\n")
    cat("Detected toe-offs at frames:", paste(toe_offs, collapse = ", "), "\n")
    
    # Handle the initial stance phase logic
    if (length(heel_strikes) > 0 && length(toe_offs) > 0) {
      first_heel_strike <- heel_strikes[1]
      first_toe_off <- toe_offs[1]
      
      if (first_heel_strike < first_toe_off) {
        # Start stance phase at the first heel strike
        stance_data$in_stance[first_heel_strike:first_toe_off] <- TRUE
        cat("Initial stance phase set from heel strike at frame", first_heel_strike, 
            "to toe-off at frame", first_toe_off, "\n")
      } else {
        # No heel strike before toe-off: assume stance phase starts from frame 1
        stance_data$in_stance[1:first_toe_off] <- TRUE
        cat("Initial stance phase set from frame 1 to first toe-off at frame", first_toe_off, "\n")
      }
    } else if (length(toe_offs) > 0) {
      # No heel strike detected, assume stance phase starts from frame 1
      stance_data$in_stance[1:toe_offs[1]] <- TRUE
      cat("No heel strikes detected, initial stance phase set from frame 1 to first toe-off at frame", 
          toe_offs[1], "\n")
    } else {
      # No heel strikes or toe-offs detected, default to full stance phase
      stance_data$in_stance[] <- TRUE
      warning("No heel strikes or toe-offs detected, assuming full stance phase for this trial.")
    }
    
    # Debugging output to confirm initial stance phase
    cat("Initial stance phase (first 100 frames):\n")
    print(head(stance_data[1:100, c("frame", "in_stance", "heel_strike", "toe_off")]))
    
    
    # State machine for stance phase detection (after the first toe-off)
    current_stance <- FALSE
    stance_frame_count <- 0
    
    for (i in seq_along(stance_data$frame)) {
      # Skip frames before the first toe-off (already handled above)
      if (i <= first_toe_off) {
        next
      }
      
      # Calculate stance frame limits based on speed and condition
      min_stance_frames <- case_when(
        stance_data$speed[1] <= 3.0 ~ 15,    # Longer minimum for slow speeds
        stance_data$speed[1] <= 4.5 ~ 12,    # Medium for medium speeds
        TRUE ~ 8                            # Original for fast speeds
      )
      
      max_stance_frames <- case_when(
        stance_data$condition[1] == "PD" ~ 60,  # Longer maximum for PD
        stance_data$speed[1] <= 3.0 ~ 65,       # Longer for slow speeds
        stance_data$speed[1] <= 4.5 ~ 50,       # Medium for medium speeds
        TRUE ~ 45                               # Original for fast speeds
      )
      
      if (!is.na(stance_data$heel_strike[i]) && stance_data$heel_strike[i]) {
        if (!current_stance) {  # Only start new stance if not already in stance
          current_stance <- TRUE
          stance_frame_count <- 0
        }
      } else if (!is.na(stance_data$toe_off[i]) && stance_data$toe_off[i] && current_stance) {
        # Only end stance if we've had a minimum duration
        if (stance_frame_count >= min_stance_frames) {
          current_stance <- FALSE
        }
      } else if (current_stance && stance_frame_count >= max_stance_frames) {
        # Force end stance if we exceed maximum duration
        current_stance <- FALSE
      }
      
      if (current_stance) {
        stance_frame_count <- stance_frame_count + 1
      }
      
      stance_data$in_stance[i] <- current_stance
    }
    
    # Create individual diagnostic plots for Y Position
    p <- ggplot(stance_data, aes(x = frame)) +
      geom_line(aes(y = ankle_y, color = "Ankle"), size = 0.8) +
      geom_line(aes(y = foot_y, color = "Foot"), size = 0.8) +
      geom_rect(data = filter(stance_data, in_stance),
                aes(xmin = frame, xmax = frame + 1, ymin = -Inf, ymax = Inf, fill = "Stance Phase"),
                alpha = 0.2) +
      geom_point(data = filter(stance_data, heel_strike),
                 aes(y = ankle_y, color = "Heel Strike"), size = 2) +
      geom_point(data = filter(stance_data, toe_off),
                 aes(y = foot_y, color = "Toe Off"), size = 2) +
      scale_color_manual(values = c("Ankle" = "#20EAABFF", 
                                    "Foot" = "#3E9BFEFF",
                                    "Heel Strike" = "#F6C33AFF", 
                                    "Toe Off" = "#F05B12FF")) +
      scale_fill_manual(values = c("Stance Phase" = "#B7EC7A")) +
      scale_y_continuous(
        trans = "reverse",
        breaks = function(limits) {
          y_max <- ceiling(max(abs(limits)) / 5) * 5
          y_min <- floor(min(abs(limits)) / 5) * 5
          seq(y_max, y_min, by = -5)  # Generate breaks from high to low
        },
        labels = function(breaks) {
          rev(breaks)  # Reverse the labels
        }
      ) +
      labs(title = paste(result$metadata$condition, result$metadata$speed, "km/h"),
           y = "Y Position (pixels)",
           x = "Frame Number",
           fill = NULL,
           color = NULL) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal"
      )
    
    # Store plot
    result$plot <- p
    
    # Calculate temporal metrics
    stance_percentage <- mean(stance_data$in_stance, na.rm = TRUE) * 100
    
    # Calculate transitions from FALSE to TRUE in in_stance to get number of steps
    num_actual_strides <- sum(diff(as.numeric(stance_data$in_stance)) == 1, na.rm = TRUE)
    
    # Find frame numbers where stance phase begins
    stance_transitions <- which(diff(as.numeric(stance_data$in_stance)) == 1)
    
    # Calculate raw stride time (in seconds)
    stride_times_raw <- diff(stance_transitions) / 60
    
    # Convert speed to m/s
    speed_ms <- result$metadata$speed / 3.6  # Convert km/h to m/s
    
    # Calculate stride lengths
    stride_lengths <- stride_times_raw * speed_ms  # This gives stride length in meters
    
    
    # Temporal metrics
    temporal_metrics <- tibble(
      speed = result$metadata$speed,
      condition = result$metadata$condition,
      stance_phase = stance_percentage,
      num_strides = num_actual_strides,
      mean_stride_time = mean(stride_times_raw, na.rm = TRUE),
      sd_stride_time = sd(stride_times_raw, na.rm = TRUE),
      mean_stride_length = mean(stride_lengths, na.rm = TRUE),
      sd_stride_length = sd(stride_lengths, na.rm = TRUE),
      stride_times = list(stride_times_raw),
      stride_lengths = list(stride_lengths)
    )
    
    return(list(metrics = temporal_metrics, plot = p))
  })
  
  # Extract metrics and plots
  metrics <- lapply(temporal_data, function(x) x$metrics) %>% bind_rows()
  plots <- lapply(temporal_data, function(x) x$plot)
  
  # Create combined plot
  combined_plot <- wrap_plots(plots, ncol = 3) +  # 3 columns for 9 plots
    plot_layout(guides = "collect") +  # Collect all legends into one
    plot_annotation(
      title = "Stance Phase Analysis: All Speeds and Conditions",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Title formatting
        legend.position = "bottom",  # Place legend at the bottom
        legend.direction = "horizontal",  # Align legend items horizontally
        legend.title = element_blank(),  # Remove the legend titles
        legend.text = element_text(size = 10)  # Adjust text size if needed
      )
    )
  
  # Save the combined plot
  ggsave("gait_analysis/plots/S2/S2_stance_raw_plots_y_adjusted_long.png", combined_plot, 
         width = 10.66, height = 13.33, dpi = 300)
  
  # Display the combined plot
  print(combined_plot)
  
  # Return the metrics
  return(metrics)
}


#-------------------------------------------------------------------------------
#### DATA MANAGEMENT
#-------------------------------------------------------------------------------

save_analysis_results <- function(results, base_path = "gait_analysis") {
  # Create directory structure
  dir.create(file.path(base_path, "raw_data"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(base_path, "processed_data"), recursive = TRUE, showWarnings = FALSE)
  
  # Create filename with underscore instead of decimal
  filename_speed <- gsub("\\.", "_", sprintf("%.1f", results$metadata$speed))
  # Remove trailing zero for whole numbers (4.0 becomes 4)
  filename_speed <- sub("_0$", "", filename_speed)
  
  # Save raw marker data as CSV
  write_csv(results$marker_data, 
            file.path(base_path, "raw_data", 
                      paste0(results$metadata$subject_id, "_", 
                             filename_speed, "_",
                             results$metadata$condition, "_raw.csv")))
  
  # Save processed data as RDS
  saveRDS(results, 
          file.path(base_path, "processed_data", 
                    paste0(results$metadata$subject_id, "_", 
                           filename_speed, "_",
                           results$metadata$condition, "_processed.rds")))
}


#-------------------------------------------------------------------------------
#### ANALYSIS WORKFLOW
#-------------------------------------------------------------------------------

create_all_marker_plots <- function(data, max_frame = 560) {
  # Define marker labels
  marker_names <- c("Hip", "Knee", "Ankle", "Foot")
  markers <- c(4, 3, 2, 1)  # Ensure correct order
  
  # Define condition colors
  condition_colors <- c(
    "normal" = "#377EB8",
    "limp" = "#4DAF4A",
    "PD" = "#E41A1C"
  )
  
  # Prepare data for faceted plot
  plot_data <- data %>%
    filter(marker_num %in% markers, frame <= max_frame) %>%
    group_by(frame, speed, condition, marker_num) %>%
    summarise(y = mean(y), .groups = 'drop') %>%
    mutate(
      condition = factor(condition, levels = c("normal", "limp", "PD")),
      speed = factor(paste0(speed, " km/h"), levels = paste0(sort(unique(speed)), " km/h")),  # Ensure order
      marker = factor(
        case_when(
          marker_num == 4 ~ "Hip",
          marker_num == 3 ~ "Knee",
          marker_num == 2 ~ "Ankle",
          marker_num == 1 ~ "Foot"
        ),
        levels = marker_names  # Ensure descending order
      )
    )
  
  # Create faceted plot
  p <- ggplot(plot_data, aes(x = frame, y = y, color = condition, group = condition)) +
    geom_line(linewidth = 4, alpha = 0.8) +
    scale_color_manual(values = condition_colors) +
    scale_y_continuous(
      trans = "reverse",
      breaks = function(limits) {
        y_min <- floor(min(limits) / 5) * 5
        y_max <- ceiling(max(limits) / 5) * 5
        seq(y_min, y_max, by = 5)  # Ensure increments of 5
      }
    ) +
    facet_grid(rows = vars(marker), cols = vars(speed), scales = "free_y") +
    theme_minimal() +
    labs(x = "Frame Number", y = "Y Position (pixels)", color = "Condition") +
    theme(
      legend.position = "bottom",
      strip.text.x = element_text(size = 14),  # Speed labels (top)
      strip.text.y = element_text(size = 14, angle = 0),  # Marker labels (right)
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray95"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Overall title
  final_plot <- p + plot_annotation(
    title = "Separated Marker Trajectories - All Conditions",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )
  
  return(final_plot)
}

# Create and display the updated plot
all_plots <- create_all_marker_plots(marker_data)
print(all_plots)


# Save the plot
ggsave("gait_analysis/plots/S2/S2_marker_trajectories.png", all_plots, 
       width = 10.66, height = 13.33 , dpi = 300)

# Combine markers in each plot per speed
create_combined_marker_plots <- function(data, max_frame = 560) {
  speeds <- unique(data$speed)
  marker_names <- c("Hip", "Knee", "Ankle", "Foot")
  plots <- list()
  
  condition_colors <- c(
    "normal" = "#377EB8",
    "limp" = "#4DAF4A",
    "PD" = "#E41A1C"
  )
  
  format_speed <- function(s) {
    paste(s, "km/h")
  }
  
  for(s in speeds) {
    message("Processing speed: ", s)
    
    plot_data <- data %>% 
      filter(speed == s, marker_num %in% c(1, 2, 3, 4), frame <= max_frame) %>%
      group_by(frame, condition, marker_num) %>%
      summarise(y = mean(y, na.rm = TRUE), .groups = 'drop') %>%
      mutate(marker = factor(case_when(
        marker_num == 4 ~ "Hip",
        marker_num == 3 ~ "Knee",
        marker_num == 2 ~ "Ankle",
        marker_num == 1 ~ "Foot"
      ), levels = c("Hip", "Knee", "Ankle", "Foot")),
      # **Reorder condition factor levels**
      condition = factor(condition, levels = c("normal", "limp", "PD"))
      )
    
    if(nrow(plot_data) == 0) {
      message("No data for speed ", s, ". Skipping.")
      next  # Skip this iteration
    }
    
    p <- ggplot(plot_data, aes(x = frame, y = y, color = condition, group = interaction(condition, marker))) + 
      geom_line(size = 1.5, alpha = 0.8) +
      scale_color_manual(values = condition_colors, name = "Condition") +
      scale_y_reverse(breaks = seq(40, 120, by = 40),
                      labels = seq(120, 40, by = -40)) +
      theme_minimal() +
      labs(title = format_speed(s), x = "Frame Number", y = "Y Position (pixels)") +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Title formatting
        legend.position = "bottom",  # Ensure legend is at the bottom
        legend.direction = "horizontal",  # Arrange legend items horizontally
        legend.title = element_blank(),  # Remove the legend title
        legend.text = element_text(size = 10)  # Adjust legend text size
      ) +
      guides(
        color = guide_legend(order = 1)  # Ensure condition order is respected
      )
    
    plots[[as.character(s)]] <- p
  }
  
  if(length(plots) == 0) {
    stop("No valid plots were generated. Check marker_data for missing values.")
  }
  
  combined_plots <- wrap_plots(plots, ncol = 3) +
    plot_layout(guides = "collect") +  # **Forces legend to be collected at the bottom**
    plot_annotation(
      title = "Marker Trajectories - All Conditions",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Title formatting
        legend.position = "bottom",  # Force legend at the bottom
        legend.direction = "horizontal",  # Align legend items horizontally
        legend.title = element_blank(),  # Remove the legend title
        legend.text = element_text(size = 10)  # Adjust text size if needed
        
      )
    )
  
  return(combined_plots)
}

# Run the function
all_combined_plots <- create_combined_marker_plots(marker_data)
print(all_combined_plots)


# Print structure of the plots list

# Save the plot
ggsave("gait_analysis/plots/S2/S2_marker_trajectories_combined_long.png", all_combined_plots, 
       width = 10.66, height = 13.33, dpi = 300)


#-------------------------------------------------------------------------------
#### KINEMATIC ANALYSIS SUMMARY
#-------------------------------------------------------------------------------

# Process kinematic data
kinematics <- analyze_gait_kinematics(results_list)

# Create kinematic summary
kinematic_summary <- kinematics %>%
  group_by(condition, speed) %>%
  summarise(
    # Knee angle calculations
    knee_ROM = max(knee_angle, na.rm = TRUE) - min(knee_angle, na.rm = TRUE),
    peak_flexion = max(knee_angle, na.rm = TRUE),
    peak_extension = min(knee_angle, na.rm = TRUE),
    mean_knee_angle = mean(knee_angle, na.rm = TRUE),
    sd_knee_angle = sd(knee_angle, na.rm = TRUE),
    cv_knee_angle = (sd_knee_angle / mean_knee_angle) * 100,
    
    # Ankle angle calculations
    ankle_ROM = max(ankle_angle, na.rm = TRUE) - min(ankle_angle, na.rm = TRUE),
    peak_ankle_flexion = max(ankle_angle, na.rm = TRUE),
    peak_ankle_extension = min(ankle_angle, na.rm = TRUE),
    mean_ankle_angle = mean(ankle_angle, na.rm = TRUE),
    sd_ankle_angle = sd(ankle_angle, na.rm = TRUE),
    cv_ankle_angle = (sd_ankle_angle / mean_ankle_angle) * 100,
    
    # CV for ROM
    cv_knee_ROM = (sd(knee_angle, na.rm = TRUE) / knee_ROM) * 100,
    cv_ankle_ROM = (sd(ankle_angle, na.rm = TRUE) / ankle_ROM) * 100,
    .groups = 'drop'
  )

print(kinematic_summary)

kinematic_summary %>% 
  arrange(condition, speed) %>%
  print(n = Inf, width = Inf)

# Create Subject ID
kinematic_summary$subject_id <- "S2"  

# Save kinematic summary
write_csv(kinematic_summary, 
          "gait_analysis/complete_analysis/S2_kinematic_data.csv")


#-------------------------------------------------------------------------------
#### KINEMATIC ANALYSIS VISUALIZATION
#-------------------------------------------------------------------------------

# Reorder the conditions so 'normal' comes first
kinematic_summary$condition <- factor(kinematic_summary$condition, 
                                      levels = c("normal", "limp", "PD"))

kinematics$condition <- factor(kinematics$condition, 
                               levels = c("normal", "limp", "PD"))

# Boxplot of knee angle RoM
k1 <- ggplot(kinematics, 
             aes(x = factor(speed), y = knee_angle, fill = condition)) +
  geom_boxplot(position = position_dodge(0.9)) +
  labs(
       x = "Speed (km/h)",
       y = "Knee Angle (degrees)",
       fill = "Condition") +
  theme_minimal() +
  scale_fill_manual(values = c("#382A54FF", "#357BA2FF", "#60CEACFF")) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# Boxplot of ankle angle RoM
k2 <- ggplot(kinematics, 
             aes(x = factor(speed), y = ankle_angle, fill = condition)) +
  geom_boxplot(position = position_dodge(0.9)) +
  labs(
       x = "Speed (km/h)",
       y = "Ankle Angle (degrees)",
       fill = "Condition") +
  theme_minimal() +
  scale_fill_manual(values = c("#382A54FF", "#357BA2FF", "#60CEACFF")) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# Combine all both plots
combined_kinematic_plots <- (k1 + k2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Raw Joint Angle: Knee and Ankle",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 10)
    )
  )

print(combined_kinematic_plots)


ggsave("gait_analysis/plots/S2/S2_kinematic_parameters.png", combined_kinematic_plots, 
       width = 15, height = 15, dpi = 300)

#-------------------------------------------------------------------------------
#### KINEMATIC ANALYSIS: ENHANCED STATISTICAL ANALYSIS
#-------------------------------------------------------------------------------

# 1. Effect of speed within each condition
speed_effects <- kinematics %>%
  group_by(condition) %>%
  summarise(
    # Stride length
    stride_speed_f = summary(aov(stride_length ~ factor(speed)))[[1]]["F value"][[1]][1],
    stride_speed_p = summary(aov(stride_length ~ factor(speed)))[[1]]["Pr(>F)"][[1]][1],
    
    # Knee angle
    knee_speed_f = summary(aov(knee_angle ~ factor(speed)))[[1]]["F value"][[1]][1],
    knee_speed_p = summary(aov(knee_angle ~ factor(speed)))[[1]]["Pr(>F)"][[1]][1],
    .groups = 'drop'
  ) %>%
  mutate(across(ends_with("_p"), ~format.pval(., digits = 3)))

# 2. Effect of condition at each speed
condition_effects <- kinematics %>%
  group_by(speed) %>%
  summarise(
    # Stride length
    stride_cond_f = summary(aov(stride_length ~ condition))[[1]]["F value"][[1]][1],
    stride_cond_p = summary(aov(stride_length ~ condition))[[1]]["Pr(>F)"][[1]][1],
    
    # Knee angle
    knee_cond_f = summary(aov(knee_angle ~ condition))[[1]]["F value"][[1]][1],
    knee_cond_p = summary(aov(knee_angle ~ condition))[[1]]["Pr(>F)"][[1]][1],
    .groups = 'drop'
  ) %>%
  mutate(across(ends_with("_p"), ~format.pval(., digits = 3)))

# 3. Calculate percent differences from normal walking
percent_diff <- kinematics %>%
  group_by(speed) %>%
  summarise(
    PD_stride_diff = (mean(stride_length[condition == "PD"]) - 
                        mean(stride_length[condition == "normal"])) / 
      mean(stride_length[condition == "normal"]) * 100,
    
    limp_stride_diff = (mean(stride_length[condition == "limp"]) - 
                          mean(stride_length[condition == "normal"])) / 
      mean(stride_length[condition == "normal"]) * 100,
    
    PD_knee_diff = (mean(knee_angle[condition == "PD"]) - 
                      mean(knee_angle[condition == "normal"])) / 
      mean(knee_angle[condition == "normal"]) * 100,
    
    limp_knee_diff = (mean(knee_angle[condition == "limp"]) - 
                        mean(knee_angle[condition == "normal"])) / 
      mean(knee_angle[condition == "normal"]) * 100,
    .groups = 'drop'
  ) %>%
  mutate(across(-speed, round, 1))

# Print results
cat("\nEffect of Speed on Kinematics within each Condition:\n")
print(speed_effects, width = Inf)

cat("\nEffect of Condition at each Speed:\n")
print(condition_effects, width = Inf)

cat("\nPercent Differences from Normal Walking:\n")
print(percent_diff, width = Inf)



#-------------------------------------------------------------------------------
#### TEMPORAL ANALYSIS: SUMMARY
#-------------------------------------------------------------------------------

# After loading your results
temporal_results <- analyze_temporal_parameters(results_list)

# Print summary statistics
temporal_summary <- temporal_results %>%
  group_by(condition, speed) %>%
  summarise(
    mean_stance = mean(stance_phase, na.rm = TRUE),  # Already a percentage
    num_steps = mean(num_strides, na.rm = TRUE),
    mean_stride_time = mean(mean_stride_time, na.rm = TRUE),
    sd_stride_time = mean(sd_stride_time, na.rm = TRUE),
    # Assuming each video is ~600 frames (10 seconds at 60fps)
    duration_seconds = 10,  # Set to actual video duration
    cadence = (num_steps / duration_seconds) * 60,  # Steps per minute
    stance_normalized = mean_stance / 100,  # Convert percentage to decimal
    stride_frequency = 60 / mean_stride_time,
    mean_stride_length = mean_stride_length,
    sd_stride_length = sd_stride_length,
    # Add Coefficient of Variation (CV) calculations
    cv_stride_time = (sd_stride_time / mean_stride_time) * 100,
    cv_stride_length = (sd_stride_length / mean_stride_length) * 100,
    .groups = 'drop'
  ) %>%
  mutate(
    gait_speed = (cadence * mean_stride_length) /60
  )

# Print nicely formatted results
print(temporal_summary %>%
        arrange(speed, condition, n = Inf, width = Inf)) 

temporal_summary %>% 
  arrange(speed, condition) %>%
  print(n = Inf, width = Inf)

# Create subject ID column
temporal_summary$subject_id <- "S2" 

# Save temporal summary
write_csv(temporal_summary, 
          "gait_analysis/complete_analysis/S2_temporal_data.csv") 

#-------------------------------------------------------------------------------
#### TEMPORAL ANALYSIS: VISUALIZATION
#-------------------------------------------------------------------------------

# Reorder the conditions so 'normal' comes first
temporal_summary$condition <- factor(temporal_summary$condition, 
                                     levels = c("normal", "limp", "PD"))

# Reorder the conditions in stride_times_df to match temporal_summary
stride_times_df$condition <- factor(stride_times_df$condition, 
                                    levels = c("normal", "limp", "PD"))

# 1. Line plot showing stance phase changes across speeds
t1 <- ggplot(temporal_summary, 
             aes(x = speed, y = mean_stance, color = condition, group = condition)) +
  geom_line(size = 1.2, show.legend = FALSE) +
  geom_point(size = 3, show.legend = FALSE) +
  labs(title = "Stance Phase",
       x = "Speed (km/h)",
       y = "Stance Phase (%)",
       color = "Condition") +
  theme_minimal() +
  scale_color_manual(values = c("#382A54FF", "#357BA2FF", "#60CEACFF")) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# 2. Grouped bar plot for stance phase changes
t2 <- ggplot(temporal_summary, 
             aes(x = factor(speed), y = mean_stance, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Stance Phase",
       x = "Speed (km/h)",
       y = "Stance Phase (%)",
       fill = "Condition") +
  theme_minimal() +
  scale_fill_manual(values = c("#382A54FF", "#357BA2FF", "#60CEACFF")) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# 3. Box plot for stride time
t3 <- ggplot(stride_times_df, 
             aes(x = factor(speed), y = stride_time, fill = condition)) +
  geom_boxplot(position = position_dodge(0.9), width = 0.7, show.legend = FALSE) +
  labs(title = "Stride Times",
       x = "Speed (km/h)",
       y = "Stride Time (s)",
       fill = "Condition") +
  theme_minimal() +
  scale_fill_manual(values = c("#382A54FF", "#357BA2FF", "#60CEACFF")) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# 4. Box plot for stride length
t4 <- ggplot(temporal_summary, 
             aes(x = factor(speed), y = mean_stride_length, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_stride_length - sd_stride_length,
                    ymax = mean_stride_length + sd_stride_length),
                position = position_dodge(0.9),
                width = 0.25) +
  labs(title = "Stride Length",
       x = "Speed (km/h)",
       y = "Stride Length (m)",
       fill = "Condition") +
  theme_minimal() +
  scale_fill_manual(values = c("#382A54FF", "#357BA2FF", "#60CEACFF")) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# 5. Bar chart for cadence
t5 <- ggplot(temporal_summary, 
             aes(x = factor(speed), y = cadence, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Cadence",
       x = "Speed (km/h)",
       y = "Cadence (steps/min)",
       fill = "Condition") +
  theme_minimal() +
  scale_fill_manual(values = c("#382A54FF", "#357BA2FF", "#60CEACFF")) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# 6. Bar chart for gait speed
t6 <- ggplot(temporal_summary, 
             aes(x = factor(speed), y = gait_speed, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Gait Speed",
       x = "Treadmill Speed (km/h)",
       y = "Gait Speed (m/s)",
       fill = "Condition") +
  theme_minimal() +
  scale_fill_manual(values = c("#382A54FF", "#357BA2FF", "#60CEACFF")) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# Combine all plots
combined_temp_plots <- (t1 + t3) / (t2 + t4) / (t5 + t6) +
  plot_layout(guides = "collect") &
  plot_annotation(
    title = "Spatiotemporal Parameters",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 10)
    )
  )

print(combined_temp_plots)

# Save the combined plot
ggsave("gait_analysis/plots/S2/S2_temporal_parameters.png", combined_temp_plots, 
       width = 15, height = 15, dpi = 300)



#-------------------------------------------------------------------------------
#### TEMPORAL ANALYSIS: ENHANCED STATISTICAL ANALYSIS
#-------------------------------------------------------------------------------

# Overall condition differences (across all speeds)
condition_stats <- temporal_results %>%
  summarise(
    # ANOVA test for differences between conditions
    cadence_anova = anova(lm(cadence ~ condition))$`Pr(>F)`[1],
    step_time_anova = anova(lm(mean_step_time ~ condition))$`Pr(>F)`[1],
    
    # Mean values and percent differences
    pd_vs_normal_cadence = (mean(cadence[condition == "PD"]) - 
                              mean(cadence[condition == "normal"])) / 
      mean(cadence[condition == "normal"]) * 100,
    limp_vs_normal_cadence = (mean(cadence[condition == "limp"]) - 
                                mean(cadence[condition == "normal"])) / 
      mean(cadence[condition == "normal"]) * 100
  )

# Speed dependency analysis
speed_dependency <- temporal_results %>%
  group_by(condition) %>%
  summarise(
    # How much does speed affect cadence?
    cadence_speed_correlation = cor(speed, cadence),
    step_time_speed_correlation = cor(speed, mean_step_time),
    
    # Average values across speeds
    mean_cadence = mean(cadence),
    sd_cadence = sd(cadence),
    cv_cadence = (sd_cadence / mean_cadence) * 100,
    
    mean_step_time = mean(mean_step_time),
    sd_step_time = sd(mean_step_time),
    cv_step_time = (sd_step_time / mean_step_time) * 100
  )

# Print results
cat("\nOverall Condition Statistics:\n")
print(condition_stats)

cat("\nSpeed Dependency Analysis:\n")
print(speed_dependency)

# Create visualization of speed effects
speed_effect_plot <- ggplot(temporal_results, 
                            aes(x = speed, y = cadence, color = condition)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(title = "Speed Effects on Cadence by Condition",
       x = "Speed (km/h)",
       y = "Cadence (steps/min)")

print(speed_effect_plot)

# Speed-Dependency Analysis with visualization
speed_dependency <- temporal_results %>%
  group_by(condition) %>%
  summarise(
    cadence_slope = coef(lm(cadence ~ speed))[2],
    step_time_slope = coef(lm(mean_step_time ~ speed))[2],
    # R-squared values to show goodness of fit
    cadence_r2 = summary(lm(cadence ~ speed))$r.squared,
    step_time_r2 = summary(lm(mean_step_time ~ speed))$r.squared
  )

print(speed_dependency)

# Normalized parameters with simple visualization
speed_effect_plot <- ggplot(temporal_results, 
                            aes(x = speed, y = cadence, color = condition)) +
  geom_point(size = 3) +
  geom_line() +
  theme_minimal() +
  labs(title = "Speed Effects on Cadence by Condition",
       x = "Speed (km/h)",
       y = "Cadence (steps/min)")

print(speed_effect_plot)








