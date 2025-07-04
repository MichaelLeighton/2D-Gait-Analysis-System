#### Gait Analysis Framework
library(tidyverse)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(zoo)
library(here)

#-------------------------------------------------------------------------------
#### KINEMATIC ANALYSIS
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
              
              acos(dot_product / magnitudes) * (180/pi)
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
#### TEMPORAL ANALYSIS
#-------------------------------------------------------------------------------
analyze_temporal_parameters <- function(results_list) {
  # Create lookup tables for thresholds
  threshold_lookup <- tribble(
    ~condition, ~speed, ~stability, ~heel_height, ~velocity, ~toe_prev_vel, ~toe_curr_vel,
    "limp",     2.5,    2.2,       0.6,          -0.04,     0.22,         0.08,
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
        potential_toe_off = !is.na(foot_y_vel) &
          abs(lag(foot_y_vel, 1, default = 0)) < toe_prev_threshold &
          foot_y_vel > toe_curr_threshold &
          lead(foot_y_vel, 1, default = 0) > toe_curr_threshold &
          lead(foot_y_vel, 2, default = 0) > toe_curr_threshold,
        
        
        toe_off = potential_toe_off &
          !lag(potential_toe_off, 1, default = FALSE) &
          !lag(potential_toe_off, 2, default = FALSE) &
          !lag(potential_toe_off, 3, default = FALSE) &
          !lag(potential_toe_off, 4, default = FALSE)
      )
    
    # Initialize state machine variables
    current_stance <- FALSE
    stance_frame_count <- 0
    
    # Add debugging output
    cat(sprintf("\nProcessing %s at %.1f km/h\n", 
                result$metadata$condition, 
                result$metadata$speed))
    cat("Number of frames:", nrow(stance_data), "\n")
    cat("Number of heel strikes detected:", sum(stance_data$heel_strike, na.rm = TRUE), "\n")
    cat("Number of toe-offs detected:", sum(stance_data$toe_off, na.rm = TRUE), "\n")
    
    # State machine for stance phase
    for(i in 1:nrow(stance_data)) {
      # Calculate stance frame limits based on speed and condition
      min_stance_frames <- case_when(
        stance_data$speed[1] <= 3.0 ~ 15,    # Longer minimum for slow speeds
        stance_data$speed[1] <= 4.5 ~ 12,    # Medium for medium speeds
        TRUE ~ 10                            # Original for fast speeds
      )
      
      max_stance_frames <- case_when(
        stance_data$condition[1] == "PD" ~ 60,  # Longer maximum for PD
        stance_data$speed[1] <= 3.0 ~ 55,       # Longer for slow speeds
        stance_data$speed[1] <= 4.5 ~ 50,       # Medium for medium speeds
        TRUE ~ 45                               # Original for fast speeds
      )
      
      if(!is.na(stance_data$heel_strike[i]) && stance_data$heel_strike[i]) {
        if(!current_stance) {  # Only start new stance if not already in stance
          current_stance <- TRUE
          stance_frame_count <- 0
        }
      } else if(!is.na(stance_data$toe_off[i]) && 
                stance_data$toe_off[i] && 
                current_stance) {
        # Only end stance if we've had a minimum duration
        if(stance_frame_count >= min_stance_frames) {
          current_stance <- FALSE
        }
      } else if(current_stance && stance_frame_count >= max_stance_frames) {
        # Force end stance if we exceed maximum duration
        current_stance <- FALSE
      }
      
      if(current_stance) {
        stance_frame_count <- stance_frame_count + 1
      }
      
      stance_data$in_stance[i] <- current_stance
    }
    
    # Create individual diagnostic plots for Y Position
    p <- ggplot(stance_data, aes(x = frame)) +
      geom_line(aes(y = ankle_y, color = "Ankle"), size = 0.8) +
      geom_line(aes(y = foot_y, color = "Foot"), size = 0.8) +
      geom_rect(data = filter(stance_data, in_stance),
                aes(xmin = frame, xmax = frame + 1, 
                    ymin = -Inf, ymax = Inf),
                fill = "green", alpha = 0.2) +
      geom_point(data = filter(stance_data, heel_strike),
                 aes(y = ankle_y, color = "Heel Strike"), size = 2) +
      geom_point(data = filter(stance_data, toe_off),
                 aes(y = foot_y, color = "Toe Off"), size = 2) +
      scale_color_manual(values = c("Ankle" = "red", "Foot" = "blue",
                                    "Heel Strike" = "green", "Toe Off" = "purple")) +
      scale_y_reverse() +
      labs(title = paste(result$metadata$condition, result$metadata$speed, "km/h"),
           y = "Y Position (pixels)",
           x = "Frame Number") +
      theme_minimal()
    
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
    plot_layout(guides = "collect") &  # Collect all legends into one
    theme(legend.position = "bottom")  # Put legend at bottom
  
  # Save the combined plot
  ggsave("all_position_plots.png", combined_plot, 
         width = 20, height = 15, dpi = 300)
  
  # Display the combined plot
  print(combined_plot)
  
  # Return the metrics
  return(metrics)
}

# Function to plot temporal parameters
plot_temporal_parameters <- function(temporal_results) {
  # Create subplots
  p1 <- ggplot(temporal_results, aes(x = factor(speed), y = cadence, fill = condition)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    labs(title = "Cadence Across Speeds and Conditions",
         x = "Speed (km/h)",
         y = "Cadence (steps/min)")
  
  p2 <- ggplot(temporal_results, aes(x = factor(speed), y = mean_step_time, fill = condition)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    labs(title = "Step Time Across Speeds and Conditions",
         x = "Speed (km/h)",
         y = "Step Time (s)")
  
  p3 <- ggplot(temporal_results, aes(x = factor(speed), y = stance_phase, fill = condition)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    labs(title = "Stance Phase Analysis",
         x = "Speed (km/h)",
         y = "Stance Phase (%)")
  
  # Combine plots using patchwork
  combined_plot <- p1 / p2 / p3 +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  return(combined_plot)
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
