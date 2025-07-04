#### Gait Analysis
### 20 Dec 2024
### Michael Leighton
### Subject 2, 4 km/h, normal

#### PACKAGES
library(magick)
library(imager)
library(tidyverse)
library(av)

#-------------------------------------------------------------------------------
#### UTILITY FUNCTIONS

# Function to extract numerical vector from image
extract_numeric_vector <- function(img) {
  if(is.character(img)){
    img <- image_read(img)
  }
  img_vector <- img %>%
    image_convert("gray") %>%
    image_data(channels = "gray") %>%
    as.numeric()
  return(img_vector)
}

# Function to standardize image to 0-1 scale
standardize_image <- function(img_vector) {
  standardized <- (img_vector - min(img_vector)) / (max(img_vector) - min(img_vector))
  return(standardized)
}

# Function to calculate new dimensions while maintaining aspect ratio
calculate_new_dimensions <- function(width, height, target_size = 100) {
  aspect_ratio <- width / height
  if (width > height) {
    new_width <- target_size
    new_height <- round(target_size / aspect_ratio)
  } else {
    new_height <- target_size
    new_width <- round(target_size * aspect_ratio)
  }
  return(list(width = new_width, height = new_height))
}

# Function to resize image
resize_image <- function(image_path, target_size = 100) {
  img <- image_read(image_path)
  original_width <- image_info(img)$width
  original_height <- image_info(img)$height
  
  new_dims <- calculate_new_dimensions(original_width, original_height, target_size)
  
  resized_img <- image_resize(img, paste0(new_dims$width, "x", new_dims$height))
  return(resized_img)
}

# Function to find the brightest spot in the image
find_brightest_spot <- function(img_matrix) {
  brightest_spot <- which(img_matrix == max(img_matrix), arr.ind = TRUE)[1,]
  return(brightest_spot)
}


#-------------------------------------------------------------------------------
#### CORE GAIT ANALYSIS FUNCTIONS

# Function to measure average brightness of a block
measure_block_brightness <- function(img_matrix, x, y, block_size) {
  half_size <- floor(block_size / 2)
  x_range <- max(1, x - half_size):min(ncol(img_matrix), x + half_size)
  y_range <- max(1, y - half_size):min(nrow(img_matrix), y + half_size)
  block <- img_matrix[y_range, x_range]
  return(mean(block))
}

# Function to check marker size
check_marker_size <- function(img_matrix, position, min_radius = 2, max_radius = 5, threshold = 0.5) {
  # Get dimensions of image
  height <- nrow(img_matrix)
  width <- ncol(img_matrix)
  
  # Create a mask for checking sizes
  y_center <- position[1]
  x_center <- position[2]
  
  # Count bright pixels within max radius
  pixel_count <- 0
  total_pixels <- 0
  
  for(dy in -max_radius:max_radius) {
    for(dx in -max_radius:max_radius) {
      y <- y_center + dy
      x <- x_center + dx
      
      # Check if pixel is within image bounds
      if(y > 0 && y <= height && x > 0 && x <= width) {
        # Calculate distance from center
        distance <- sqrt(dx^2 + dy^2)
        
        if(distance <= max_radius) {
          total_pixels <- total_pixels + 1
          if(img_matrix[y, x] > threshold) {
            pixel_count <- pixel_count + 1
          }
        }
      }
    }
  }
  
  # Calculate area ratio
  area_ratio <- pixel_count / total_pixels
  
  # Check if size is within acceptable range
  min_area_ratio <- pi * min_radius^2 / (pi * max_radius^2)
  max_area_ratio <- 1.2 * pi * max_radius^2 / (pi * max_radius^2)  # 20% tolerance
  
  return(area_ratio >= min_area_ratio && area_ratio <= max_area_ratio)
}

# Function to check marker surroundings
check_marker_surroundings <- function(img_matrix, position, inner_radius = 5, outer_radius = 10, verbose = FALSE) {
  inner_brightness <- measure_block_brightness(img_matrix, position[2], position[1], inner_radius * 2)
  outer_brightness <- measure_block_brightness(img_matrix, position[2], position[1], outer_radius * 2) - 
    measure_block_brightness(img_matrix, position[2], position[1], (inner_radius * 2 - 1))
  if(verbose) {
    print(paste("Inner brightness:", inner_brightness, "Outer brightness:", outer_brightness))
  }
  
  result <- inner_brightness > 0.05 && inner_brightness > (outer_brightness * 1.1)
  if(verbose) {
    print(paste("Brightness check result:", result))
  }
  return(result)
}

# Function to eliminate nearby solutions
erase_marker <- function(img_matrix, marker_position, radius = 5) {
  x_range <- max(1, marker_position[2] - radius):min(ncol(img_matrix), marker_position[2] + radius)
  y_range <- max(1, marker_position[1] - radius):min(nrow(img_matrix), marker_position[1] + radius)
  img_matrix[y_range, x_range] <- 0
  return(img_matrix)
}

# Function to rank solutions by spot brightness or contrast using a Threshold setting
find_brightness_threshold <- function(img_matrix, n_top = 10, verbose = FALSE) {
  sorted_brightness <- sort(img_matrix, decreasing = TRUE)[1:n_top]
  # Tweak values below for spot brightness
  threshold <- max(sorted_brightness[1] * 0.91, 0.8)
  
  # Verbose parameter for testing process
  if(verbose) {
    print("Top brightness values:")
    print(sorted_brightness)
    print(paste("Chosen threshold:", threshold))
    plot(sorted_brightness, type = 'b', main = "Brightness Values")
    abline(h = threshold, col = "red", lty = 2)
  }
  
  return(threshold)
}

# Enhanced circularity check function
check_circularity <- function(img_matrix, position, radius = 5, threshold = 0.3, verbose = FALSE, is_lower_leg = FALSE) {
  # Get dimensions
  height <- nrow(img_matrix)
  width <- ncol(img_matrix)
  
  # Create vectors to store boundary points
  boundary_points_x <- numeric()
  boundary_points_y <- numeric()
  
  # Adjust minimum points required for lower leg/foot markers
  min_points <- if(is_lower_leg) 6 else 10
  
  # Find boundary points
  for(angle in seq(0, 2*pi, length.out = 36)) {
    found_edge <- FALSE
    for(r in seq(1, radius * 1.5, 0.5)) {
      x <- round(position[2] + r * cos(angle))
      y <- round(position[1] + r * sin(angle))
      
      if(x > 0 && x <= width && y > 0 && y <= height) {
        if(img_matrix[y, x] < threshold) {
          # Simplified edge detection for lower leg/foot markers
          if(is_lower_leg) {
            found_edge <- TRUE
            boundary_points_x <- c(boundary_points_x, x)
            boundary_points_y <- c(boundary_points_y, y)
            break
          } else {
            # Original stricter edge detection for other markers
            if(r < radius * 1.4) {
              next_x <- round(position[2] + (r + 0.5) * cos(angle))
              next_y <- round(position[1] + (r + 0.5) * sin(angle))
              if(next_x > 0 && next_x <= width && next_y > 0 && next_y <= height) {
                if(img_matrix[next_y, next_x] < threshold) {
                  found_edge <- TRUE
                  boundary_points_x <- c(boundary_points_x, x)
                  boundary_points_y <- c(boundary_points_y, y)
                  break
                }
              }
            }
          }
        }
      }
    }
  }
  
  if(verbose) {
    cat(sprintf("Found %d boundary points\n", length(boundary_points_x)))
    if(length(boundary_points_x) < min_points) {
      cat("Not enough boundary points found - returning 0\n")
    }
  }
  
  if(length(boundary_points_x) < min_points) {
    return(0)
  }
  
  # Calculate center of found points
  center_x <- mean(boundary_points_x)
  center_y <- mean(boundary_points_y)
  
  # Calculate distances from center to boundary points
  distances <- sqrt((boundary_points_x - center_x)^2 + (boundary_points_y - center_y)^2)
  
  # Calculate circularity as variation in radius
  circularity <- 1 - (sd(distances) / mean(distances))
  
  return(circularity)
}

# Function to check the brightness of a ring around the marker (to ensure it is dark around the marker)
check_surrounding_ring <- function(img_matrix, position, inner_radius, outer_radius) {
  height <- nrow(img_matrix)
  width <- ncol(img_matrix)
  y_center <- position[1]
  x_center <- position[2]
  
  ring_pixels <- numeric()
  
  for(dy in -outer_radius:outer_radius) {
    for(dx in -outer_radius:outer_radius) {
      y <- y_center + dy
      x <- x_center + dx
      
      if(y > 0 && y <= height && x > 0 && x <= width) {
        distance <- sqrt(dx^2 + dy^2)
        # Only consider pixels in the ring between inner and outer radius
        if(distance >= inner_radius && distance <= outer_radius) {
          ring_pixels <- c(ring_pixels, img_matrix[y, x])
        }
      }
    }
  }
  
  # Return average brightness of the ring
  return(mean(ring_pixels))
}

# Function to ensure marker is isolated
check_marker_isolation <- function(img_matrix, position, radius, isolation_factor = 2, threshold = 0.3) {
  
  wide_radius <- radius * isolation_factor
  total_white_area <- 0
  marker_area <- pi * radius^2
  
  height <- nrow(img_matrix)
  width <- ncol(img_matrix)
  
  # Count white pixels in the area
  for(dy in -wide_radius:wide_radius) {
    for(dx in -wide_radius:wide_radius) {
      y <- position[1] + dy
      x <- position[2] + dx
      
      if(y > 0 && y <= height && x > 0 && x <= width) {
        if(img_matrix[y, x] > threshold) {
          total_white_area <- total_white_area + 1
        }
      }
    }
  }
  
  # Calculate area ratio and set threshold
  area_ratio <- marker_area / total_white_area
  return(area_ratio > 0.5)
}

# Function to check for circular dark border
check_circular_border <- function(img_matrix, position, radius, verbose = FALSE) {
  # Get center brightness first
  center_brightness <- measure_block_brightness(img_matrix, position[2], position[1], radius)
  
  # Check ring around the marker
  num_points <- 36
  ring_values <- numeric()
  ring_radius <- radius * 1.1
  
  for(angle in seq(0, 2*pi, length.out = num_points)) {
    x <- round(position[2] + ring_radius * cos(angle))
    y <- round(position[1] + ring_radius * sin(angle))
    
    if(x > 0 && x <= ncol(img_matrix) && y > 0 && y <= nrow(img_matrix)) {
      ring_values <- c(ring_values, img_matrix[y, x])
    }
  }
  
  ring_mean <- mean(ring_values)
  ring_std <- sd(ring_values)
  contrast_ratio <- center_brightness / ring_mean
  
  if(verbose) {
    print(paste("Centre brightness:", center_brightness))
    print(paste("Ring mean brightness:", ring_mean))
    print(paste("Ring std dev:", ring_std))
    print(paste("Contrast ratio:", contrast_ratio))
  }
  
  # Center should be much brighter than ring (contrast_ratio > 1.5)
  # Ring should be consistent (std < 0.3)
  return(contrast_ratio > 1.5 && ring_std < 0.3)
}

# Contrast check function
check_contrast <- function(img_matrix, position, radius = 3) {
  center <- measure_block_brightness(img_matrix, position[2], position[1], radius)
  # Check surrounding ring instead of just above
  surround <- 0
  for(angle in seq(0, 2*pi, length.out = 8)) {
    x <- round(position[2] + radius * cos(angle))
    y <- round(position[1] + radius * sin(angle))
    surround <- surround + img_matrix[y, x]
  }
  surround <- surround / 8
  
  # Different contrast thresholds based on position
  contrast_threshold <- if(position[1] > 140) {  # Foot region
    1.15  # More lenient for foot
  } else if(position[1] >= 120) {  # Ankle region
    1.2
  } else {
    1.3  # Standard threshold for hip and knee
  }
  
  return(center / surround > contrast_threshold)
}

# Marker validation function
validate_marker <- function(img_matrix, position, min_radius = 3, max_radius = 8, 
                            verbose = FALSE) {
  if(verbose) print("\n=== Marker Validation Details ===")
  
  # Define anatomical zones
  hip_zone <- position[1] >= 10 && position[1] <= 50 && 
    position[2] >= 125 && position[2] <= 200
  
  knee_zone <- position[1] >= 75 && position[1] <= 100 && 
    position[2] >= 110 && position[2] <= 198
  
  lower_leg_zone <- position[1] >= 120 && position[1] <= 158 && 
    position[2] >= 80 && position[2] <= 213
  
  # Print zone checks if verbose
  if(verbose) {
    print(paste("Position: x =", position[2], "y =", position[1]))
    print(paste("In hip zone:", hip_zone))
    print(paste("In knee zone:", knee_zone))
    print(paste("In lower leg zone:", lower_leg_zone))
  }
  
  # Check if marker is in any valid anatomical zone
  if(!(hip_zone || knee_zone || lower_leg_zone)) {
    if(verbose) {
      print(paste("Rejected - position (x:", position[2], "y:", position[1], 
                  ") not in any valid anatomical zone"))
    }
    return(FALSE)
  }
  
  # Get measurements
  center_brightness <- measure_block_brightness(img_matrix, position[2], position[1], block_size = 4)
  if(verbose) print(paste("Center brightness:", center_brightness))
  
  # Brightness thresholds - can be specific to each zone
  # Tweak for center brightness
  min_brightness <- if(lower_leg_zone) {
    0.22  # More lenient for ankle/foot
  } else if(knee_zone) {
    0.50
  } else {
    0.54  # Stricter for hip
  }
  
  if(center_brightness < min_brightness) {
    if(verbose) print("Failed brightness check")
    return(FALSE)
  }
  
  # Size check - can be specific to each zone
  size_multiplier <- if(position[1] > 140) { # Specifically for foot
    0.55  # More lenient for foot
  } else if(lower_leg_zone) {
    0.7 # for ankle
  } else {
    0.85 # Standard for hip and knee
  }
  
  size_valid <- check_marker_size(img_matrix, position, 
                                  min_radius * size_multiplier, 
                                  max_radius * 1.2)
  if(verbose) print(paste("Size check result:", size_valid))
  if(!size_valid) {
    if(verbose) print("Failed size check")
    return(FALSE)
  }
  
  # Circularity check - can be specific to each zone
  circularity <- check_circularity(img_matrix, position, max_radius, verbose = verbose, is_lower_leg = lower_leg_zone)
  if(verbose) print(paste("Circularity value:", circularity))
  
  # Tweak circularity threshold
  circularity_threshold <- if(lower_leg_zone) {
    0.44  # More lenient for ankle/foot
  } else if(knee_zone) {
    0.6   # Slightly more lenient for knee
  } else {
    0.6   # Stricter for hip
  }
  
  if(circularity < circularity_threshold) {
    if(verbose) {
      print("Failed circularity check")
      print(paste("Circularity threshold was:", circularity_threshold))
    }
    return(FALSE)
  }
  
  # Contrast check
  if(!check_contrast(img_matrix, position)) {
    if(verbose) print("Failed contrast check")
    return(FALSE)
  }
  
  return(TRUE)
}

# Find markers in a frame using adaptive threshold
find_markers <- function(img_matrix, brightness_threshold, marker_size = 10, 
                         max_iterations = 100, edge_margin = 5, verbose = FALSE) {
  markers <- list()
  iteration <- 0
  min_radius <- marker_size/5
  max_radius <- marker_size/2
  
  if(verbose) {
    print(paste("Using min_radius:", min_radius))
    print(paste("Using max_radius:", max_radius))
    print(paste("Initial brightness threshold:", brightness_threshold))
  }
  
  while(max(img_matrix) > brightness_threshold && iteration < max_iterations) {
    iteration <- iteration + 1
    brightest_spot <- which(img_matrix == max(img_matrix), arr.ind = TRUE)[1,]
    
    if(verbose) {
      print(paste("\nIteration:", iteration))
      print(paste("Examining spot at row:", brightest_spot[1], "col:", brightest_spot[2]))
      print(paste("Spot brightness:", max(img_matrix)))
    }
    
    # Adjust edge margin based on Y position
    if(brightest_spot[1] > 160) { # For foot markers
      edge_margin <- 3
    } else {
      edge_margin <- 5
    }
    
    # Edge margin check
    if(brightest_spot[1] > edge_margin && 
       brightest_spot[1] < (nrow(img_matrix) - edge_margin) &&
       brightest_spot[2] > edge_margin && 
       brightest_spot[2] < (ncol(img_matrix) - edge_margin)) {
      
      if(verbose) print("\nStarting marker validation...")
      
      # Validate marker using all criteria
      if(validate_marker(img_matrix, brightest_spot, 
                         min_radius = min_radius, 
                         max_radius = max_radius,
                         verbose = verbose)) {
        markers <- c(markers, list(brightest_spot))
        if(verbose) {
          print("MARKER ACCEPTED")
          print(paste("Marker added at row:", brightest_spot[1], "col:", brightest_spot[2]))
        }
      } else {
        if(verbose) {
          print("MARKER REJECTED")
        }
      }
    } else {
      if(verbose) print("Spot rejected - too close to edge")
    }
    
    # Erase the examined area
    img_matrix <- erase_marker(img_matrix, brightest_spot, radius = max_radius * 2)
  }
  
  return(markers)
}

# Diagnostic visualization function
visualize_marker_validation <- function(img_matrix, position, 
                                        min_radius = 3, max_radius = 8) {
  # Create validation results
  size_valid <- check_marker_size(img_matrix, position, min_radius, max_radius)
  circularity <- check_circularity(img_matrix, position, max_radius)
  inner_brightness <- measure_block_brightness(img_matrix, position[2], position[1], max_radius * 2)
  outer_brightness <- measure_block_brightness(img_matrix, position[2], position[1], max_radius * 4) -
    measure_block_brightness(img_matrix, position[2], position[1], max_radius * 2)
  
  # Create visualization
  p <- ggplot() +
    # Original image
    geom_raster(data = expand.grid(x = 1:ncol(img_matrix), y = 1:nrow(img_matrix)) %>%
                  mutate(value = as.vector(img_matrix)),
                aes(x = x, y = y, fill = value)) +
    scale_fill_gradient(low = "black", high = "white") +
    geom_circle(aes(x0 = position[2], y0 = position[1], r = min_radius),
                color = "red", fill = NA) +
    geom_circle(aes(x0 = position[2], y0 = position[1], r = max_radius),
                color = "blue", fill = NA) +
    coord_fixed() +
    scale_y_reverse() +
    theme_minimal() +
    labs(title = sprintf("Marker Validation Results\nSize valid: %s\nCircularity: %.2f\nContrast ratio: %.2f",
                         size_valid, circularity, inner_brightness/outer_brightness))
  
  return(p)
}

# Process a single frame
process_frame <- function(frame_path, target_size = 100, marker_size = 10, verbose = FALSE) {
  # Resized image while maintaining aspect ratio
  img <- resize_image(frame_path, target_size)
  dims <- image_info(img)
  # Tweak blur (sigma = X, higher blur reduces noise, lower blur differentiates marker from background)
  img_matrix <- img %>%
    image_convert("gray") %>%
    image_blur(radius = 1, sigma = 0.35) %>%
    image_data(channels = "gray") %>%
    as.numeric() %>%
    matrix(nrow = dims$height, ncol = dims$width) %>%
    standardize_image()
  
  # Verbose parameter for testing process
  if(verbose) {
    print("Image matrix summary:")
    print(summary(as.vector(img_matrix)))
  }
  
  brightness_threshold <- find_brightness_threshold(img_matrix, verbose = verbose)
  markers <- find_markers(img_matrix, brightness_threshold, 
                          marker_size = marker_size,
                          verbose = verbose)
  
  # Verbose parameter for testing process
  if(verbose) {
    print(paste("Number of markers found:", length(markers)))
  }
  
  return(markers)
}

# Function to test a single frame
test_single_frame <- function(frame_path, target_size = 300, marker_size = 10, verbose = FALSE) {
  # Process the frame
  markers <- process_frame(frame_path, target_size, marker_size = marker_size, verbose = verbose)
  
  # Load and resize the image maintaining aspect ratio
  img <- resize_image(frame_path, target_size)
  dims <- image_info(img)
  
  # Convert image to a matrix
  img_matrix <- img %>%
    image_convert(colorspace = "gray") %>%
    image_data() %>%
    as.numeric() %>%
    matrix(nrow = dims$height, ncol = dims$width) %>%
    standardize_image()
  
  # Create data frame for plotting
  img_df <- data.frame(
    x = rep(seq_len(dims$width), each = dims$height),
    y = rep(seq_len(dims$height), times = dims$width)
  )
  img_df$value <- as.vector(img_matrix)
  
  # Verbose parameter for testing process
  if(verbose) {
    print("Detected markers:")
    for(i in seq_along(markers)) {
      print(paste("Marker", i, "at position:", 
                  "x =", markers[[i]][2], 
                  "y =", markers[[i]][1]))
    }
  }
  
  # Create the plot
  p <- ggplot(img_df, aes(x = x, y = y)) +
    geom_raster(aes(fill = value)) +
    scale_fill_gradient(low = "black", high = "white") +
    coord_fixed() +
    scale_y_reverse() +
    theme_minimal() +
    labs(title = "Single Frame Analysis", x = "X coordinate", y = "Y coordinate")
  
  # Add markers to the plot if any were found
  if (length(markers) > 0) {
    marker_df <- data.frame(
      x = sapply(markers, function(m) m[2]),
      y = sapply(markers, function(m) m[1])
    )
    p <- p + geom_point(data = marker_df, aes(x = x, y = y), 
                        color = 'red', size = 3, shape = 4)
  }
  
  # Return both the plot and the markers
  return(list(plot = p, markers = markers, img_matrix = img_matrix))
}

# Nearest neighbour tracking function
track_markers <- function(current_markers, previous_markers, prev_velocities = NULL) {
  
  if(length(previous_markers) == 0) {
    # Keep first frame initialization exactly as is
    cat("Initializing first frame markers\n")
    
    markers_df <- data.frame(
      x = sapply(current_markers, function(m) m[2]),
      y = sapply(current_markers, function(m) m[1])
    )
    
    initial_markers <- list(
      list(name = "hip", num = 1, y = max(markers_df$y)),
      list(name = "knee", num = 2, y = sort(markers_df$y, decreasing=TRUE)[2]),
      list(name = "ankle", num = 3, y = sort(markers_df$y, decreasing=TRUE)[3]),
      list(name = "foot", num = 4, y = min(markers_df$y))
    )
    
    marker_assignments <- numeric(nrow(markers_df))
    for(i in 1:nrow(markers_df)) {
      distances <- sapply(initial_markers, function(m) abs(markers_df$y[i] - m$y))
      marker_assignments[i] <- initial_markers[[which.min(distances)]]$num
      cat(sprintf("Marker at y=%d, x=%d assigned number %d\n", 
                  markers_df$y[i], markers_df$x[i], marker_assignments[i]))
    }
    
    markers_df$marker_num <- marker_assignments
    return(markers_df)
  }
  
  # For subsequent frames
  max_distance <- 40
  
  current_df <- data.frame(
    x = sapply(current_markers, function(m) m[2]),
    y = sapply(current_markers, function(m) m[1])
  )
  
  # Calculate distances
  distances <- matrix(NA, nrow = nrow(current_df), ncol = nrow(previous_markers))
  for(i in 1:nrow(current_df)) {
    for(j in 1:nrow(previous_markers)) {
      distances[i,j] <- sqrt(
        (current_df$x[i] - previous_markers$x[j])^2 +
          (current_df$y[i] - previous_markers$y[j])^2
      )
    }
  }
  
  # Initialize assignments
  assignments <- rep(NA, nrow(current_df))
  used_prev_markers <- rep(FALSE, nrow(previous_markers))
  
  # Simple distance-based assignment first
  for(i in 1:nrow(current_df)) {
    if(sum(!used_prev_markers) == 0) break
    if(sum(!is.na(distances[i,!used_prev_markers])) == 0) next
    
    min_dist <- min(distances[i,!used_prev_markers], na.rm = TRUE)
    if(min_dist <= max_distance) {
      prev_idx <- which.min(distances[i,!used_prev_markers])
      prev_idx <- which(!used_prev_markers)[prev_idx]
      assignments[i] <- previous_markers$marker_num[prev_idx]
      used_prev_markers[prev_idx] <- TRUE
    }
  }
  
  # Handle any unassigned markers based on Y position
  if(any(is.na(assignments))) {
    unassigned <- which(is.na(assignments))
    available_nums <- setdiff(1:4, assignments[!is.na(assignments)])
    
    # Sort by Y position to maintain relative positions
    unassigned_y <- current_df$y[unassigned]
    sorted_indices <- unassigned[order(unassigned_y, decreasing = TRUE)]
    
    for(i in seq_along(sorted_indices)) {
      if(i <= length(available_nums)) {
        assignments[sorted_indices[i]] <- available_nums[i]
      }
    }
  }
  
  # Create final data frame
  markers_df <- data.frame(
    marker_num = assignments,
    x = current_df$x,
    y = current_df$y
  )
  
  return(markers_df)
}

# Function to process and visualize multiple frames
analyze_video_frames <- function(frames, target_size = 300, verbose = FALSE) {
  results <- list()
  
  # Process each frame
  for(i in seq_along(frames)) {
    if(verbose) {
      cat(sprintf("\nProcessing frame %d of %d\n", i, length(frames)))
    }
    
    # Process the frame - pass verbose parameter
    result <- test_single_frame(frames[i], 
                                target_size = target_size, 
                                marker_size = 10,
                                verbose = verbose)
    
    # Store the results
    results[[i]] <- list(
      frame_number = i,
      markers = result$markers,
      plot = result$plot,
      num_markers = length(result$markers)
    )
    
    if(verbose) {
      # Print number of markers found
      cat(sprintf("Found %d markers in frame %d\n", length(result$markers), i))
      
      # Display the plot
      print(result$plot)
      
      # Small delay to make visualization manageable
      Sys.sleep(0.5)
    }
  }
  
  # Summary statistics - always show regardless of verbose setting
  num_markers_per_frame <- sapply(results, function(x) x$num_markers)
  cat("\nSummary Statistics:\n")
  cat(sprintf("Average markers per frame: %.2f\n", mean(num_markers_per_frame)))
  cat(sprintf("Frames with 4 markers: %d of %d\n", 
              sum(num_markers_per_frame == 4), 
              length(frames)))
  
  return(results)
}

# Function to analyze marker detection across frames
analyze_marker_detection <- function(frames) {
  # Store results in a list
  results <- vector("list", length(frames))  # Pre-allocate the list
  
  # Data frame to store all marker positions
  all_markers <- data.frame(
    frame = integer(),
    marker_num = integer(),
    x = numeric(),
    y = numeric()
  )
  
  # Initialize previous_markers as NULL for tracking
  previous_markers <- NULL
  
  # Process each frame
  for(i in seq_along(frames)) {
    cat(sprintf("\rProcessing frame %d of %d", i, length(frames)))
    
    result <- test_single_frame(frames[i], target_size = 300)
    
    # Store the result, including number of markers found
    results[[i]] <- list(
      frame_number = i,
      markers = result$markers,
      num_markers = length(result$markers)
    )
    
    if(length(result$markers) > 0) {
      # Use tracking function with velocity information
      current_markers_df <- track_markers(result$markers, previous_markers)
      
      # Store for next iteration
      previous_markers <- current_markers_df
      
      # Add to all_markers data frame
      frame_markers <- data.frame(
        frame = i,
        marker_num = current_markers_df$marker_num,
        x = current_markers_df$x,
        y = current_markers_df$y
      )
      all_markers <- rbind(all_markers, frame_markers)
    }
  }
  
  # Create summary visualizations
  plots <- list()
  
  # Calculate markers per frame from results list
  markers_per_frame <- sapply(results, function(x) x$num_markers)
  
  # 1. Histogram of markers per frame
  plots$histogram <- ggplot(data.frame(markers = markers_per_frame), aes(x = markers)) +
    geom_histogram(binwidth = 1, fill = "blue", alpha = 0.7) +
    labs(title = "Distribution of Markers Detected per Frame",
         x = "Number of Markers",
         y = "Frequency") +
    theme_minimal()
  
  # 2. Marker positions over time
  plots$positions <- ggplot(all_markers, aes(x = frame, y = y, color = factor(marker_num))) +
    geom_point(alpha = 0.5) +
    geom_line(alpha = 0.3) +
    labs(title = "Marker Y-Positions Over Time",
         x = "Frame Number",
         y = "Y Position",
         color = "Marker Number") +
    scale_y_reverse() +
    theme_minimal()
  
  # 3. Heatmap of marker positions
  plots$heatmap <- ggplot(all_markers, aes(x = x, y = y)) +
    geom_bin2d(bins = 30) +
    scale_fill_viridis_c() +
    labs(title = "Heatmap of Marker Positions",
         x = "X Position",
         y = "Y Position") +
    theme_minimal() +
    coord_fixed()
  
  # 4. Success rate over time
  success_rate <- data.frame(
    frame = seq_along(frames),
    markers = markers_per_frame
  )
  plots$success_rate <- ggplot(success_rate, aes(x = frame, y = markers)) +
    geom_line() +
    geom_hline(yintercept = 4, color = "red", linetype = "dashed") +
    labs(title = "Number of Markers Detected Over Time",
         x = "Frame Number",
         y = "Number of Markers") +
    theme_minimal()
  
  # Calculate summary statistics
  stats <- list(
    total_frames = length(frames),
    mean_markers = mean(markers_per_frame),
    sd_markers = sd(markers_per_frame),
    perfect_detection = sum(markers_per_frame == 4),
    perfect_detection_rate = sum(markers_per_frame == 4) / length(frames) * 100,
    missed_frames = sum(markers_per_frame < 4),
    extra_detections = sum(markers_per_frame > 4)
  )
  
  # Return all results
  return(list(
    plots = plots,
    stats = stats,
    marker_data = all_markers,
    raw_results = results
  ))
}


#-------------------------------------------------------------------------------
#### MAIN WORKFLOW
#-------------------------------------------------------------------------------
# Create directory for the video
output_dir <- "output_directory_S2_4_normal"
dir.create(output_dir, showWarnings = FALSE)

# Set up parameters
video_path <- "S2_4_normal.MOV"

# Extract frames
av_video_images(video_path, output_dir, format = "jpg")

# Get all frames
frames <- list.files(output_dir, pattern = "*.jpg", full.names = TRUE)
print(paste("Number of frames extracted:", length(frames)))


#### FIRST FRAME CHECK (verbose mode)
first_frame <- frames[423]
test_result <- test_single_frame(first_frame, 
                                 target_size = 300, 
                                 marker_size = 10, 
                                 verbose = TRUE)



# Display results
print(test_result$plot)
cat("Image matrix summary:\n")
print(summary(as.vector(test_result$img_matrix)))

cat("\nDetected marker positions:\n")
if (length(test_result$markers) > 0) {
  for (i in seq_along(test_result$markers)) {
    cat(sprintf("Marker %d: row %d, column %d\n", 
                i, test_result$markers[[i]][1], test_result$markers[[i]][2]))
  }
} else {
  cat("No markers detected.\n")
}



#### MULTIPLE FRAME CHECK (verbose mode)

# Process first 10 frames with detailed output
test_frames <- frames[55:65]
test_results <- analyze_video_frames(test_frames, verbose = TRUE)



#### FULL ANALYSIS OF MARKERS (quiet/non-verbose mode)

# Clear workspace of previous results
rm(analysis_results)

# Run the analysis on all frames
analysis_results <- analyze_marker_detection(frames)

# Print summary statistics
cat("\nAnalysis Summary:\n")
cat("================\n")
cat(sprintf("Total frames analyzed: %d\n", 
            analysis_results$stats$total_frames))
cat(sprintf("Average markers per frame: %.2f (SD: %.2f)\n", 
            analysis_results$stats$mean_markers,
            analysis_results$stats$sd_markers))
cat(sprintf("Frames with exactly 4 markers: %d (%.1f%%)\n", 
            analysis_results$stats$perfect_detection,
            analysis_results$stats$perfect_detection_rate))
cat(sprintf("Frames with missing markers: %d\n", 
            analysis_results$stats$missed_frames))
cat(sprintf("Frames with extra detections: %d\n", 
            analysis_results$stats$extra_detections))

# Identify problematic frames (uncomment with Cmd+Shift+C )
{
  # Add to analysis_results to identify problematic frames
  problem_frames <- analysis_results$marker_data %>%
    group_by(frame) %>%
    filter(n() < 4) %>%
    select(frame, marker_num, x, y)
  
  print("Frames with excess markers:")
  print(n = 740, problem_frames)
  
  # Create plots list
  plots <- list()
  
  # Create visualization
  plots$problem_frames <- ggplot(problem_frames, aes(x = x, y = y)) +
    geom_point(alpha = 0.5) +
    labs(title = "Marker Positions in Frames with Missing Detections",
         x = "X Position",
         y = "Y Position") +
    scale_y_reverse() +
    theme_minimal()
  
  # Display the plot
  print(plots$problem_frames)
}

# Display analysis plots
print(analysis_results$plots$histogram)
print(analysis_results$plots$positions)
print(analysis_results$plots$heatmap)
print(analysis_results$plots$success_rate)


#-------------------------------------------------------------------------------
#### SAVE ANALYSIS RESULTS
#-------------------------------------------------------------------------------
# Video metadata
metadata <- list(
  subject_id = "S2",
  speed = 4,
  condition = "normal",
  date = "2024-12-20"
)

# Create the complete results object
results <- list(
  metadata = metadata,
  marker_data = analysis_results$marker_data,
  statistics = analysis_results$stats
)

# Function to save the results
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
  
  # Save complete results as RDS
  saveRDS(results, 
          file.path(base_path, "processed_data", 
                    paste0(results$metadata$subject_id, "_", 
                           filename_speed, "_",
                           results$metadata$condition, "_processed.rds")))
}

# Save the results
save_analysis_results(results)

#-------------------------------------------------------------------------------
##### DEBUG: FUNCTIONS
#-------------------------------------------------------------------------------
# Function to evaluate a single frame's marker detection
evaluate_frame_detection <- function(detected_markers, expected_positions, max_distance = 15, verbose = TRUE) {
  if(verbose) {
    cat("\nDetected markers:\n")
    for(m in detected_markers) {
      cat(sprintf("  y=%d, x=%d\n", m[1], m[2]))
    }
    cat("\nExpected markers:\n")
    for(m in expected_positions) {
      cat(sprintf("  y=%d, x=%d\n", m[1], m[2]))
    }
  }
  
  true_positives <- 0
  false_positives <- length(detected_markers)
  false_negatives <- length(expected_positions)
  matched_pairs <- list()
  
  # For each detected marker, find closest expected position
  for(detected in detected_markers) {
    min_dist <- Inf
    matched_idx <- NULL
    
    for(i in seq_along(expected_positions)) {
      expected <- expected_positions[[i]]
      if(all(expected == c(-999, -999))) next  # Skip already matched
      
      dist <- sqrt((detected[1] - expected[1])^2 + (detected[2] - expected[2])^2)
      if(dist < min_dist && dist < max_distance) {
        min_dist <- dist
        matched_idx <- i
      }
    }
    
    if(!is.null(matched_idx)) {
      true_positives <- true_positives + 1
      false_positives <- false_positives - 1
      false_negatives <- false_negatives - 1
      matched_pairs[[length(matched_pairs) + 1]] <- list(
        detected = detected,
        expected = expected_positions[[matched_idx]],
        distance = min_dist
      )
      expected_positions[[matched_idx]] <- c(-999, -999)  # Mark as matched
    }
  }
  
  if(verbose) {
    cat("\nMatched pairs:\n")
    for(pair in matched_pairs) {
      cat(sprintf("  Detected(y=%d,x=%d) -> Expected(y=%d,x=%d) [dist=%.2f]\n",
                  pair$detected[1], pair$detected[2],
                  pair$expected[1], pair$expected[2],
                  pair$distance))
    }
    
    cat("\nUnmatched expected positions:\n")
    for(exp in expected_positions) {
      if(!all(exp == c(-999, -999))) {
        cat(sprintf("  y=%d, x=%d\n", exp[1], exp[2]))
      }
    }
  }
  
  precision <- if(true_positives + false_positives > 0) {
    true_positives / (true_positives + false_positives)
  } else {
    0
  }
  
  recall <- if(true_positives + false_negatives > 0) {
    true_positives / (true_positives + false_negatives)
  } else {
    0
  }
  
  f1 <- if(precision + recall > 0) {
    2 * (precision * recall) / (precision + recall)
  } else {
    0
  }
  
  return(list(
    true_positives = true_positives,
    false_positives = false_positives,
    false_negatives = false_negatives,
    precision = precision,
    recall = recall,
    f1_score = f1,
    matched_pairs = matched_pairs
  ))
}


# Function to test a set of parameters on multiple frames
test_parameters <- function(frames, expected_markers, params, verbose = FALSE) {
  # Modify validate_marker function with current parameters
  validate_marker <- function(img_matrix, position, min_radius = 3, max_radius = 8, 
                              verbose = FALSE) {
    if(verbose) print("\n=== Marker Validation Details ===")
    
    # False detection zones
    if((position[2] >= 102 && position[2] <= 108 && 
        position[1] >= 94 && position[1] <= 121) ||
       (position[2] >= 105 && position[2] <= 130 && 
        position[1] >= 140 && position[1] <= 150)) {
      if(verbose) print("Rejected false detection zone")
      return(FALSE)
    }
    
    # Position constraints
    if(position[2] < 95 || position[2] > 185 || 
       position[1] < 45 || position[1] > 170) {
      if(verbose) print("Failed position check")
      return(FALSE)
    }
    
    # Size check
    size_multiplier <- if(position[1] > params$foot_height) {
      params$foot_size_mult
    } else {
      params$standard_size_mult
    }
    
    size_valid <- check_marker_size(img_matrix, position, 
                                    min_radius * size_multiplier, 
                                    max_radius * params$max_size_mult)
    
    if(!size_valid) return(FALSE)
    
    # Brightness checks
    center_brightness <- measure_block_brightness(img_matrix, position[2], position[1], 
                                                  block_size = 4)
    
    min_brightness <- if(position[1] > params$foot_height) {
      params$foot_brightness
    } else {
      params$standard_brightness
    }
    
    if(center_brightness < min_brightness) return(FALSE)
    
    # Circularity check
    circularity <- check_circularity(img_matrix, position, max_radius)
    
    circularity_threshold <- if(position[1] > params$foot_height) {
      params$foot_circularity
    } else {
      params$standard_circularity
    }
    
    if(circularity < circularity_threshold) return(FALSE)
    
    # Contrast check
    if(!check_contrast(img_matrix, position)) return(FALSE)
    
    return(TRUE)
  }
  
  # Test each frame
  results <- lapply(seq_along(frames), function(i) {
    if(verbose) cat(sprintf("\rTesting frame %d/%d", i, length(frames)))
    
    # Process frame
    frame_result <- test_single_frame(frames[i], target_size = 300, 
                                      marker_size = 10, verbose = FALSE)
    
    # Evaluate detection
    evaluate_frame_detection(frame_result$markers, expected_markers[[i]])
  })
  
  # Average results
  avg_results <- Reduce(`+`, lapply(results, unlist)) / length(results)
  
  return(avg_results)
}

# Main parameter tuning function
tune_marker_detection <- function(test_frames, expected_markers) {
  # Create parameter grid
  foot_heights <- c(145, 150, 155)
  foot_brightnesses <- c(0.15, 0.18, 0.20)
  foot_circularities <- c(0.25, 0.28, 0.30)
  foot_size_mults <- c(0.65, 0.68, 0.70)
  
  # Generate all combinations
  param_sets <- list()
  for(h in foot_heights) {
    for(b in foot_brightnesses) {
      for(c in foot_circularities) {
        for(s in foot_size_mults) {
          param_sets[[length(param_sets) + 1]] <- list(
            foot_height = h,
            foot_brightness = b,
            standard_brightness = 0.3,
            foot_circularity = c,
            standard_circularity = 0.45,
            foot_size_mult = s,
            standard_size_mult = 0.85,
            max_size_mult = 1.3
          )
        }
      }
    }
  }
  
  # Test each parameter set
  results <- lapply(param_sets, function(params) {
    cat("\nTesting parameters:\n")
    print(params)
    c(params, test_parameters(test_frames, expected_markers, params, verbose = TRUE))
  })
  
  # Sort results by F1 score
  results <- results[order(sapply(results, function(x) x$f1_score), decreasing = TRUE)]
  
  return(results)
}


#-------------------------------------------------------------------------------
#### DEBUG: WORKFLOW
#-------------------------------------------------------------------------------
# Example: Select frames 8, 50, 100, 150, 200
test_frame_indices <- c(8, 50, 100, 150, 200)
test_frames <- frames[8]

# Record expected marker positions for these frames
expected_markers <- list(
  # Frame 8
  list(
    c(48, 155),   # Hip
    c(98, 138),   # Knee
    c(157, 141),  # Ankle
    c(160, 126)   # Foot
  )
)

# Run tuning
results <- tune_marker_detection(test_frames, expected_markers)

# Look at results
print(results)



#-------------------------------------------------------------------------------
#### SINGLE IMAGE ANALYSIS (currently not needed)
#-------------------------------------------------------------------------------
# Load and resize the image
image_path <- "whitespot50.jpg"
resized_img <- resize_image(image_path, width = 100, height = 100)

# Extract numeric vector and standardize
img_vector <- extract_numeric_vector(resized_img)
std_img_vector <- standardize_image(img_vector)

# Convert to matrix
img_matrix <- matrix(std_img_vector, nrow = 100, ncol = 100)

# Print the matrix
print(img_matrix)

# Convert matrix to long format for plot
img_df <- expand.grid(x = 1:ncol(img_matrix), y = 1:nrow(img_matrix))
img_df$value <- as.vector(img_matrix)

# Find the coordinates of the white dot
white_dot <- find_brightest_spot(img_matrix)
print(paste("Brightest spot found at: row", white_dot[1], "column", white_dot[2]))

# Measure average brightness around the brightest spot
block_brightness <- measure_block_brightness(img_matrix, white_dot[2], white_dot[1], block_size = 5)
print(paste("Average brightness around the brightest spot:", round(block_brightness, 3)))

# Create a more informative brightness map
brightness_map <- img_matrix

# Convert brightness_map to a data frame suitable for ggplot
brightness_df <- expand.grid(x = 1:ncol(brightness_map), y = 1:nrow(brightness_map))
brightness_df$value <- as.vector(brightness_map)

# Plot the brightness map
p_brightness <- ggplot(brightness_df, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = 'black', high = 'white',
                      breaks = c(min(brightness_map), max(brightness_map)),
                      labels = c("Min", "Max")) +
  theme_minimal() +
  labs(title = "Brightness Map", x = "X coordinate", y = "Y coordinate") +
  coord_equal() +
  scale_y_reverse()

print(p_brightness)

# Main plot
p_main <- ggplot(img_df, aes(x = y, y = x)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'black', high = 'white') +
  geom_point(data = data.frame(x = white_dot[2], y = white_dot[1]),
             aes(x = x, y = y), color = 'red', size = 3, shape = 4) +
  coord_fixed() +
  scale_y_reverse() +
  theme_minimal() +
  labs(title = "Gait Analysis Image", x = "X coordinate", y = "Y coordinate") +
  theme(legend.position = "none")

# Display plot
print(p_main)

