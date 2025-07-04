#### Combined Complete Gait Analysis
#### All subjects compared
#### Author: Michael Leighton
#### Date: 31 Dec 2024

# Packages
library(tidyverse)
library(ggpubr)
library(patchwork)
library(zoo)
library(here)
library(viridis)
library(ez) # For ANOVA
library(rstatix) # For post-hoc tests
library(effectsize) # for effect sizes
library(fmsb) # For Radar Plots (not needed now)
library(scales) # For alpha function
library(knitr)
library(kableExtra)
library(webshot)
library(officer)
library(MASS)
library(caret)
library(pROC)
library(dplyr)

#-------------------------------------------------------------------------------
#### LOAD AND COMBINE SUBJECT DATA
#-------------------------------------------------------------------------------

# Load all kinematic data files
kinematic_files <- list.files(path = "gait_analysis/complete_analysis", 
                              pattern = "*kinematic_data.csv", 
                              full.names = TRUE)

# Read and combine all kinematic files
all_kinematic_data <- kinematic_files %>%
  map_df(read_csv) %>%
  mutate(
    condition = factor(condition, levels = c("normal","limp","PD")),
    subject_id = as.factor(subject_id),
    speed = as.factor(speed)
  ) %>%
  arrange(subject_id, condition, speed) %>% 
  select(subject_id, everything(),
         -cv_knee_angle, -cv_ankle_angle, -cv_knee_ROM, -cv_ankle_ROM)

# Load all temporal data files
temporal_files <- list.files(path = "gait_analysis/complete_analysis", 
                             pattern = "*temporal_data.csv", 
                             full.names = TRUE)

# Read and combine all temporal files
all_temporal_data <- temporal_files %>%
  map_df(read_csv) %>%  
  mutate(
    condition = factor(condition, levels = c("normal","limp","PD")),
    subject_id = as.factor(subject_id),
    speed = as.factor(speed)
  ) %>%
  arrange(subject_id, condition, speed) %>% 
  select(subject_id, everything(),
         -cv_stride_time, -cv_stride_length)


#-------------------------------------------------------------------------------
#### VALIDATION ANALYSIS: ALL DATA
#-------------------------------------------------------------------------------

# SD values for key parameters
validation_sd <- all_kinematic_data %>%
  group_by(condition, speed) %>%
  summarise(
    # Knee parameters
    sd_knee_ROM = sd(knee_ROM),
    sd_knee_angle = sd(mean_knee_angle),
    
    # Ankle parameters
    sd_ankle_ROM = sd(ankle_ROM),
    sd_ankle_angle = sd(mean_ankle_angle),
    .groups = 'drop'
  )

print(validation_sd)

# Calculate CV for key parameters across subjects
validation_cv <- all_kinematic_data %>%
  group_by(condition, speed) %>%
  summarise(
    # Knee parameters
    cv_knee_ROM_between = (sd(knee_ROM) / mean(knee_ROM)) * 100,
    cv_knee_angle_between = (sd(mean_knee_angle) / mean(mean_knee_angle)) * 100,
    
    # Ankle parameters
    cv_ankle_ROM_between = (sd(ankle_ROM) / mean(ankle_ROM)) * 100,
    cv_ankle_angle_between = (sd(mean_ankle_angle) / mean(mean_ankle_angle)) * 100,
    .groups = 'drop'
  )

print(validation_cv)

# Add temporal parameters SD
temporal_sd <- all_temporal_data %>%
  group_by(condition, speed) %>%
  summarise(
    sd_stride_length_between = sd(mean_stride_length),
    sd_stride_time_between = sd(mean_stride_time),
    sd_cadence = sd(cadence),
    .groups = 'drop'
  )

print(temporal_sd)

# Add temporal parameters CV
temporal_cv <- all_temporal_data %>%
  group_by(condition, speed) %>%
  summarise(
    cv_stride_length_between = (sd(mean_stride_length) / mean(mean_stride_length)) * 100,
    cv_stride_time_between = (sd(mean_stride_time) / mean(mean_stride_time)) * 100,
    cv_cadence_between = (sd(cadence) / mean(cadence)) * 100,
    .groups = 'drop'
  )

print(temporal_cv)


#-------------------------------------------------------------------------------
#### KINEMATIC STATISTICAL ANALYSIS
#-------------------------------------------------------------------------------

# Calculate confidence intervals
kinematic_ci <- all_kinematic_data %>%
  group_by(condition, speed) %>%
  summarise(
    knee_rom_mean = mean(knee_ROM),
    knee_rom_ci_lower = t.test(knee_ROM)$conf.int[1],
    knee_rom_ci_upper = t.test(knee_ROM)$conf.int[2],
    ankle_rom_mean = mean(ankle_ROM),
    ankle_rom_ci_lower = t.test(ankle_ROM)$conf.int[1],
    ankle_rom_ci_upper = t.test(ankle_ROM)$conf.int[2]
  )

print(kinematic_ci)

# Colour scheme
condition_colours <- c("normal" = "#377EB8", 
                      "limp" = "#4DAF4A", 
                      "PD" = "#E41A1C")
speed_fills <- c("2.5" = "#E5E5E5", 
                 "4" = "#B3B3B3", 
                 "5.5" = "#808080")

# Original colours (changed just for presentation)
# condition_colours <- c("normal" = "#382A54FF", 
#                       "limp" = "#357BA2FF", 
#                       "PD" = "#60CEACFF")
# speed_fills <- c("2.5" = "#E5E5E5", 
#                  "4" = "#B3B3B3", 
#                  "5.5" = "#808080")

# 1. Repeated Measures ANOVA for kinematic parameters
# Knee ROM
knee_rom_anova <- ezANOVA(
  data = all_kinematic_data,
  dv = knee_ROM,
  wid = subject_id,
  within = .(condition, speed),
  type = 3,
  detailed = TRUE
)

# Ankle ROM
ankle_rom_anova <- ezANOVA(
  data = all_kinematic_data,
  dv = ankle_ROM,
  wid = subject_id,
  within = .(condition, speed),
  type = 3,
  detailed = TRUE
)

# Mean knee angle
mean_knee_anova <- ezANOVA(
  data = all_kinematic_data,
  dv = mean_knee_angle,
  wid = subject_id,
  within = .(condition, speed),
  type = 3,
  detailed = TRUE
)

# Mean ankle angle
mean_ankle_anova <- ezANOVA(
  data = all_kinematic_data,
  dv = mean_ankle_angle,
  wid = subject_id,
  within = .(condition, speed),
  type = 3,
  detailed = TRUE
)


# 2. Post-hoc tests
# Knee: For condition effects
knee_posthoc_condition <- all_kinematic_data %>%
  group_by(speed) %>%
  pairwise_t_test(
    knee_ROM ~ condition,
    paired = TRUE,
    p.adjust.method = "bonferroni"
  )

# Knee: For speed effects
knee_posthoc_speed <- all_kinematic_data %>%
  group_by(condition) %>%
  pairwise_t_test(
    knee_ROM ~ speed,
    paired = TRUE,
    p.adjust.method = "bonferroni"
  )

# Ankle: For condition effects
ankle_posthoc_condition <- all_kinematic_data %>%
  group_by(speed) %>%
  pairwise_t_test(
    ankle_ROM ~ condition,
    paired = TRUE,
    p.adjust.method = "bonferroni"
  )

# Ankle: For speeed effects
ankle_posthoc_speed <- all_kinematic_data %>%
  group_by(condition) %>%
  pairwise_t_test(
    ankle_ROM ~ speed,
    paired = TRUE,
    p.adjust.method = "bonferroni"
  )


# 3. Effect Sizes (Cohen's d for paired comparisons)

# Knee effect sizes
knee_effect_sizes <- all_kinematic_data %>%
  group_by(speed) %>%
  rstatix::cohens_d(knee_ROM ~ condition, paired = TRUE) %>%
  ungroup()

# Ankle effect sizes
ankle_effect_sizes <- all_kinematic_data %>%
  group_by(speed) %>%
  rstatix::cohens_d(ankle_ROM ~ condition, paired = TRUE) %>%
  ungroup()

# Print results
print(knee_rom_anova)
print(ankle_rom_anova)
print(knee_posthoc_condition)
print(knee_posthoc_speed)
print(ankle_posthoc_condition)
print(ankle_posthoc_speed)
print(knee_effect_sizes)
print(ankle_effect_sizes)

#-------------------------------------------------------------------------------
#### KINEMATIC PCA
#-------------------------------------------------------------------------------

# Prepare data for PCA (focusing on key kinematic variables)
pca_data <- all_kinematic_data %>%
  select(subject_id, condition, speed, 
         knee_ROM, peak_flexion, peak_extension, 
         ankle_ROM, peak_ankle_flexion, peak_ankle_extension)

# Perform PCA
pca_vars <- pca_data %>%
  select(-subject_id, -condition, -speed) %>%
  scale() # Standardise variables

pca_result <- prcomp(pca_vars)

print(pca_result)

# Get standard deviations from pca_result
sdev <- pca_result$sdev

# Calculate variances
variances <- sdev^2

# Calculate proportion of variance explained by each component
prop_var <- variances/sum(variances)

# Calculate cumulative proportion
cum_prop <- cumsum(prop_var)

# Convert to percentages
percentages <- prop_var * 100

summary(pca_result)

# Extract PC scores and add back grouping variables
pca_scores <- as.data.frame(pca_result$x) %>%
  bind_cols(pca_data %>% select(subject_id, condition, speed))

print(pca_scores)

# Visualize PCA results
pca_kinematic_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, colour = condition, shape = factor(speed))) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Kinematic Parameters",
       colour = "Condition",
       shape = "Speed (km/h)") +
  scale_colour_manual(values = c("#382A54FF", "#357BA2FF", "#60CEACFF")) +
  theme(
    plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

print(pca_kinematic_plot)

ggsave("gait_analysis/complete_analysis/kinematic_PCA.png", 
       plot = pca_kinematic_plot,
       width = 8,
       height = 8,
       dpi = 300)


#-------------------------------------------------------------------------------
#### KINEMATIC VISUALIZATIONS 
#-------------------------------------------------------------------------------

# Create colour scheme
condition_colours <- c("normal" = "#377EB8", 
                      "limp" = "#4DAF4A", 
                      "PD" = "#E41A1C")

# Original colour scheme (changed just for presentation)
# condition_colours <- c("normal" = "#382A54FF", 
#                       "limp" = "#357BA2FF", 
#                       "PD" = "#60CEACFF")


# Plot theme
plot_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# Modify plots to remove titles
knee_rom_plot <- ggplot(all_kinematic_data, 
                        aes(x = speed, y = knee_ROM, fill = condition)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Speed (km/h)", 
       y = "Knee ROM (degrees)") +
  scale_fill_manual(values = condition_colours,
                    name = "Gait Condition") +
  plot_theme

ankle_rom_plot <- ggplot(all_kinematic_data, 
                         aes(x = speed, y = ankle_ROM, fill = condition)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Speed (km/h)",
       y = "Ankle ROM (degrees)") +
  scale_fill_manual(values = condition_colours,
                    name = "Gait Condition") +
  plot_theme

knee_interaction <- ggplot(all_kinematic_data, 
                           aes(x = speed, y = knee_ROM, 
                               colour = condition, group = condition)) +
  stat_summary(fun = mean, geom = "point", size = 2.5, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "line", size = 0.8, show.legend = FALSE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, show.legend = FALSE) +
  labs(x = "Speed (km/h)",
       y = "Knee ROM (degrees)") +
  scale_colour_manual(values = condition_colours,
                     name = "Gait Condition") +
  plot_theme

ankle_interaction <- ggplot(all_kinematic_data, 
                            aes(x = speed, y = ankle_ROM, 
                                colour = condition, group = condition)) +
  stat_summary(fun = mean, geom = "point", size = 2.5) +
  stat_summary(fun = mean, geom = "line", size = 0.8, show.legend = FALSE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, show.legend = FALSE) +
  labs(x = "Speed (km/h)",
       y = "Ankle ROM (degrees)") +
  scale_colour_manual(values = condition_colours,
                     name = "Gait Condition") +
  plot_theme

# Combine plots with lettered labels
combined_plot <- (knee_rom_plot + knee_interaction) /
  (ankle_rom_plot + ankle_interaction) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Kinematic Changes Across Gait Conditions",
    tag_levels = 'A',
    tag_prefix = '(',
    tag_suffix = ')',
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 10),

      plot.tag = element_text(size = 12, face = "bold")
    )
  )

print(combined_plot)

# Save at appropriate size for A4
ggsave("gait_analysis/complete_analysis/kinematic_analysis_plots.png", 
       plot = combined_plot,
       width = 8,
       height = 8,
       dpi = 300)

# Effect size visualization with condition focus
{
  # effect_size_plot <- ggplot(knee_effect_sizes,
  #                          aes(x = speed, y = effsize, fill = group2)) +
  # geom_bar(stat = "identity", position = "dodge") +
  # facet_wrap(~group1) +
  # theme_minimal() +
  # labs(title = "Effect Sizes of Gait Conditions on Knee ROM",
  #      x = "Speed (km/h)",
  #      y = "Cohen's d",
  #      fill = "Comparison Condition") +
  # scale_fill_manual(values = condition_colours) +
  # theme(legend.position = "bottom")
  }


#-------------------------------------------------------------------------------
#### TEMPORAL STATISTICAL ANALYSIS
#-------------------------------------------------------------------------------

# Calculate gait cycles per condition
gait_cycles <- all_temporal_data %>%
  group_by(condition, speed) %>%
  summarise(
    cycles_per_condition = mean(stride_frequency) * (10/60), # 10 seconds per trial
    total_cycles = sum(stride_frequency * (10/60))
  )

print(gait_cycles)

# Calculate CIs for temporal parameters
temporal_ci <- all_temporal_data %>%
  group_by(condition, speed) %>%
  summarise(
    stance_mean = mean(mean_stance),
    stance_ci_lower = t.test(mean_stance)$conf.int[1],
    stance_ci_upper = t.test(mean_stance)$conf.int[2],
    stride_length_mean = mean(mean_stride_length),
    stride_length_ci_lower = t.test(mean_stride_length)$conf.int[1],
    stride_length_ci_upper = t.test(mean_stride_length)$conf.int[2],
    cadence_mean = mean(cadence),
    cadence_ci_lower = t.test(cadence)$conf.int[1],
    cadence_ci_upper = t.test(cadence)$conf.int[2]
  )

print(temporal_ci)

#-------------------------------------------------------------------------------
#### TEMPORAL PCA
#-------------------------------------------------------------------------------

# Prepare data for temporal PCA
temporal_pca_data <- all_temporal_data %>%
  select(subject_id, condition, speed, 
         mean_stance, mean_stride_time, mean_stride_length, 
         cadence, stride_frequency)

# Perform PCA on temporal variables
temporal_pca_vars <- temporal_pca_data %>%
  select(-subject_id, -condition, -speed) %>%
  scale() # Standardise variables

temporal_pca_result <- prcomp(temporal_pca_vars)

# Get standard deviations from pca_result
temporal_sdev <- temporal_pca_result$sdev

# Calculate variances
temporal_variances <- temporal_sdev^2

# Calculate proportion of variance explained by each component
temporal_prop_var <- temporal_variances/sum(temporal_variances)

# Calculate cumulative proportion
temporal_cum_prop <- cumsum(temporal_prop_var)

# Convert to percentages
temporal_percentages <- temporal_prop_var * 100

summary(temporal_pca_result)

# Extract PC scores and add back grouping variables
temporal_pca_scores <- as.data.frame(temporal_pca_result$x) %>%
  bind_cols(temporal_pca_data %>% select(subject_id, condition, speed))

# Create PCA plot for temporal parameters
temporal_pca_plot <- ggplot(temporal_pca_scores, 
                            aes(x = PC1, y = PC2, colour = condition, shape = factor(speed))) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Temporal Parameters",
       colour = "Condition",
       shape = "Speed (km/h)") +
  scale_colour_manual(values = c("#382A54FF", "#357BA2FF", "#60CEACFF")) +
  theme(
    plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

print(temporal_pca_result)
print(temporal_pca_plot)

ggsave("gait_analysis/complete_analysis/temporal_PCA.png", 
       plot = temporal_pca_plot,
       width = 8,
       height = 8,
       dpi = 300)


#-------------------------------------------------------------------------------
#### TEMPORAL ANOVA
#-------------------------------------------------------------------------------

# 1. Repeated Measures ANOVA for temporal parameters
# Stance percentage
stance_anova <- ezANOVA(
  data = all_temporal_data,
  dv = mean_stance,
  wid = subject_id,
  within = .(condition, speed),
  type = 3,
  detailed = TRUE
)

# Stride length
stride_length_anova <- ezANOVA(
  data = all_temporal_data,
  dv = mean_stride_length,
  wid = subject_id,
  within = .(condition, speed),
  type = 3,
  detailed = TRUE
)

# Cadence
cadence_anova <- ezANOVA(
  data = all_temporal_data,
  dv = cadence,
  wid = subject_id,
  within = .(condition, speed),
  type = 3,
  detailed = TRUE
)

# 2. Post-hoc tests
# For condition effects
stance_posthoc_condition <- all_temporal_data %>%
  group_by(speed) %>%
  pairwise_t_test(
    mean_stance ~ condition,
    paired = TRUE,
    p.adjust.method = "bonferroni"
  )

# For speed effects
stance_posthoc_speed <- all_temporal_data %>%
  group_by(condition) %>%
  pairwise_t_test(
    mean_stance ~ speed,
    paired = TRUE,
    p.adjust.method = "bonferroni"
  )

# 3. Effect Sizes
stance_effect_sizes <- all_temporal_data %>%
  group_by(speed) %>%
  rstatix::cohens_d(mean_stance ~ condition, paired = TRUE) %>%
  ungroup()

# Print results
print(stance_anova)
print(stride_length_anova)
print(cadence_anova)
print(stance_posthoc_condition)
print(stance_effect_sizes)


# Post-hoc tests for Stride Length by Condition
stride_length_posthoc_condition <- all_temporal_data %>%
  group_by(speed) %>%
  pairwise_t_test(
    mean_stride_length ~ condition,
    paired = TRUE,
    p.adjust.method = "bonferroni"
  )

# Post-hoc tests for Stride Length by Speed
stride_length_posthoc_speed <- all_temporal_data %>%
  group_by(condition) %>%
  pairwise_t_test(
    mean_stride_length ~ speed,
    paired = TRUE,
    p.adjust.method = "bonferroni"
  )

# Effect Sizes for Stride Length
stride_length_effect_sizes <- all_temporal_data %>%
  group_by(speed) %>%
  rstatix::cohens_d(mean_stride_length ~ condition, paired = TRUE) %>%
  ungroup()

# Post-hoc tests for Cadence by Condition
cadence_posthoc_condition <- all_temporal_data %>%
  group_by(speed) %>%
  pairwise_t_test(
    cadence ~ condition,
    paired = TRUE,
    p.adjust.method = "bonferroni"
  )

# Post-hoc tests for Cadence by Speed
cadence_posthoc_speed <- all_temporal_data %>%
  group_by(condition) %>%
  pairwise_t_test(
    cadence ~ speed,
    paired = TRUE,
    p.adjust.method = "bonferroni"
  )

# Effect Sizes for Cadence
cadence_effect_sizes <- all_temporal_data %>%
  group_by(speed) %>%
  rstatix::cohens_d(cadence ~ condition, paired = TRUE) %>%
  ungroup()

print(stride_length_posthoc_condition)
print(stride_length_posthoc_speed)
print(stride_length_effect_sizes)

print(cadence_posthoc_condition)
print(cadence_posthoc_speed)
print(cadence_effect_sizes)

#-------------------------------------------------------------------------------
#### STANCE PHASE ANALYSIS
#-------------------------------------------------------------------------------

# Calculate stance-swing ratio
phase_analysis <- all_temporal_data %>%
  group_by(subject_id, condition, speed) %>%
  summarise(
    stance_percentage = mean_stance,
    swing_percentage = 100 - mean_stance,
    stance_swing_ratio = mean_stance / (100 - mean_stance),
    .groups = 'drop'
  )

print(phase_analysis)

phase_analysis %>% 
  arrange(subject_id, speed, condition) %>%
  print(n = Inf, width = Inf)

stance_summary <- phase_analysis %>%
  group_by(condition, speed) %>%
  summarise(
    mean_stance = mean(stance_percentage),
    sd_stance = sd(stance_percentage),
    mean_ratio = mean(stance_swing_ratio),
    sd_ratio = sd(stance_swing_ratio),
    .groups = 'drop'
  )

print(stance_summary)

#-------------------------------------------------------------------------------
#### TEMPORAL VISUALIZATIONS
#-------------------------------------------------------------------------------

# Create comprehensive temporal parameter plots with mixed visualisation types
temporal_plots <- list()

# 1. Stance Phase as line plot
temporal_plots$stance <- ggplot(all_temporal_data, 
                                aes(x = speed, y = mean_stance, colour = condition, group = condition)) +
  stat_summary(fun = mean, geom = "point", size = 3, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "line", size = 1, show.legend = FALSE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, show.legend = FALSE) +
  theme_minimal() +
  labs(x = "Speed (km/h)",
       y = "Stance Phase (%)") +
  scale_colour_manual(values = condition_colours)

# 2. Stacked bar plot showing stance-swing distribution
temporal_plots$stance_swing <- ggplot(all_temporal_data, 
                                      aes(x = speed, y = mean_stance, fill = condition)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(0.9), show.legend = FALSE) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", 
                position = position_dodge(0.9), width = 0.2, show.legend = FALSE) +
  geom_hline(yintercept = 60, linetype = "dashed", colour = "gray50") +
  scale_y_continuous(limits = c(0, 100),
                     sec.axis = sec_axis(~(100-.), name = "Swing Phase (%)")) +
  theme_minimal() +
  labs(x = "Speed (km/h)",
       y = "Stance Phase (%)") +
  scale_fill_manual(values = condition_colours)

# 3. Stride Length box plots
temporal_plots$stride_length <- ggplot(all_temporal_data, 
                                       aes(x = speed, y = mean_stride_length, fill = condition)) +
  geom_boxplot(show.legend = FALSE) +
  theme_minimal() +
  labs(x = "Speed (km/h)",
       y = "Stride Length (m)") +
  scale_fill_manual(values = condition_colours)

# 4. Stride Time line plot
temporal_plots$stride_time <- ggplot(all_temporal_data, 
                                     aes(x = speed, y = mean_stride_time, colour = condition, group = condition)) +
  stat_summary(fun = mean, geom = "point", size = 3, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "line", show.legend = FALSE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, show.legend = FALSE) +
  theme_minimal() +
  labs(x = "Speed (km/h)",
       y = "Stride Time (s)") +
  scale_colour_manual(values = condition_colours)

# 5. Cadence bar plot
temporal_plots$cadence <- ggplot(all_temporal_data, 
                                 aes(x = speed, y = cadence, fill = condition)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(0.9)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(0.9), width = 0.2, show.legend = FALSE) +
  theme_minimal() +
  labs(x = "Speed (km/h)",
       y = "Steps/min") +
  scale_fill_manual(values = condition_colours) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )


# 6. Gait Speed
temporal_plots$gait_speed <- ggplot(all_temporal_data, 
                                    aes(x = speed, y = gait_speed, fill = condition)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(0.9), show.legend = FALSE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(0.9), width = 0.2, show.legend = FALSE) +
  theme_minimal() +
  labs(x = "Speed (km/h)",
       y = "Gait Speed (m/s)") +
  scale_fill_manual(values = condition_colours)

# Common theme settings for all plots
plot_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# Apply theme to each plot and adjust point/line sizes
temporal_plots$stance <- temporal_plots$stance +
  scale_colour_manual(values = condition_colours, name = "Gait Condition")

temporal_plots$stride_time <- temporal_plots$stride_time +
  scale_colour_manual(values = condition_colours, name = "Gait Condition")

# Set the legend title for plots using fill
temporal_plots$stride_length <- temporal_plots$stride_length +
  scale_fill_manual(values = condition_colours, name = "Gait Condition")

temporal_plots$cadence <- temporal_plots$cadence +
  scale_fill_manual(values = condition_colours, name = "Gait Condition")

temporal_plots$gait_speed <- temporal_plots$gait_speed +
  scale_fill_manual(values = condition_colours, name = "Gait Condition")

temporal_plots$stance_swing <- temporal_plots$stance_swing +
  scale_fill_manual(values = condition_colours, name = "Gait Condition")

# Combine plots with adjusted overall theme
combined_plot <- (temporal_plots$stance + temporal_plots$stance_swing) /
  (temporal_plots$stride_length + temporal_plots$stride_time) /
  (temporal_plots$cadence) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Spatiotemporal Parameters Across Gait Conditions",
    tag_levels = 'A',
    tag_prefix = '(',
    tag_suffix = ')',
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 10),
      plot.tag = element_text(size = 12, face = "bold")
    )
  )

print(combined_plot)

# Save at A4 size
ggsave("gait_analysis/complete_analysis/temporal_analysis_plots_lettered_no_gait_speed_presentation.png", 
       plot = combined_plot,
       width = 8,
       height = 9,
       dpi = 300)



#-------------------------------------------------------------------------------
#### TEMPORAL NORMALIZATION PLOT
#-------------------------------------------------------------------------------

#### MULTIPLE PARAMETER NORMALIZATION PLOT

# Normalize and plot multiple parameters together
temporal_plots$normalized <- all_temporal_data %>%
  mutate(across(c(mean_stride_length, cadence, mean_stance), 
                ~scale(.), 
                .names = "{.col}_scaled")) %>%
  pivot_longer(cols = ends_with("_scaled"),
               names_to = "parameter",
               values_to = "scaled_value") %>%
  mutate(parameter = factor(parameter, 
                            levels = c("mean_stance_scaled", 
                                       "cadence_scaled", 
                                       "mean_stride_length_scaled"),
                            labels = c("Stance Phase", 
                                       "Cadence", 
                                       "Stride Length"))) %>%
  ggplot(aes(x = speed, y = scaled_value, colour = condition, group = condition)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  facet_wrap(~parameter) +
  theme_minimal() +
  labs(title = "Normalized Parameter Comparison",
       x = "Speed (km/h)",
       y = "Standardized Value") +
  scale_colour_manual(values = condition_colours, name = "Gait Condition") +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 13),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )


print(temporal_plots$normalized)

# Save for A4 size
ggsave("gait_analysis/complete_analysis/normalized_temporal_plots.png", 
       plot = temporal_plots$normalized,
       width = 8,
       height = 4,
       dpi = 300,)


#-------------------------------------------------------------------------------
#### INDIVIDUAL SUBJECT DATA PLOTS
#-------------------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(ggplot2)
library(patchwork)

# Function to load and process raw marker data
load_marker_data <- function(path) {
  # All CSV files
  files <- list.files(path = path, 
                      pattern = "*_raw.csv", 
                      full.names = TRUE)
  
  # Extract condition and speed from filenames
  file_info <- data.frame(
    file = files,
    subject = str_extract(files, "S[2-5]"),
    speed = str_extract(files, "(2_5|4|5_5)"),
    condition = str_extract(files, "(normal|limp|PD)")
  )
  
  # Load and combine all files
  all_data <- map_df(1:nrow(file_info), function(i) {
    read_csv(file_info$file[i]) %>%
      mutate(
        subject = file_info$subject[i],
        speed = file_info$speed[i],
        condition = file_info$condition[i]
      )
  })
  
  return(all_data)
}

# Load data
marker_data <- load_marker_data("gait_analysis/raw_data") %>%
  mutate(condition = factor(condition, levels = c("normal", "limp", "PD")))

# Custom labels for the speeds
speed_labels <- c("2_5" = "2.5 km/h", "4" = "4 km/h", "5_5" = "5.5 km/h")

# Custom labels for trajectory plot
marker_labels <- c("1" = "foot", "2" = "ankle", "3" = "knee", "4" = "hip")

# Updated heatmap plot
heatmap_plot <- ggplot(marker_data, aes(x = x, y = y)) +
  geom_bin2d(bins = 50) +
  facet_grid(condition ~ speed, labeller = labeller(speed = speed_labels)) +
  scale_fill_viridis_c(trans = "log", breaks = c(1, 10, 100, 300)) +
  theme_minimal(base_size = 12) +
  labs(
    x = "X Position",
    y = "Y Position"
  ) +
  coord_fixed(ratio = 0.5) +
  scale_y_reverse(breaks = seq(160, 40, by = -40),
                  labels = seq(40, 160, by = 40)) +
  theme(
    plot.title = element_text(size = 13, hjust = 0.5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# 2. Trajectory Plot
trajectory_plot <- ggplot(marker_data, aes(x = x, y = y, colour = factor(marker_num))) +
  geom_path(alpha = 0.5) +
  facet_grid(condition ~ speed, labeller = labeller(speed = speed_labels)) +
  scale_colour_manual(name = "Marker",
                        breaks = c('4','3','2','1'),
                        labels = c('hip','knee','ankle','foot'),
                     values = c(
                       "1" = "#440154FF",
                       "2" = "#55C667FF",
                       "3" = "#287D8EFF",
                       "4" = "#FDE725FF"
                     )
                        ) +
  theme_minimal() +
  labs(
       x = "X Position",
       y = "Y Position") +
  coord_fixed(ratio = 0.5) +
  scale_y_reverse(breaks = seq(160, 40, by = -40),
                  labels = seq(40, 160, by = 40)) +
  theme(
    plot.title = element_text(size = 13, hjust = 0.5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# 3. Combined Average Path with Variance
{
  # average_path <- marker_data %>%
  #   group_by(condition, speed, marker_num) %>%
  #   summarise(
  #     mean_x = mean(x),
  #     mean_y = mean(y),
  #     sd_x = sd(x),
  #     sd_y = sd(y),
  #     .groups = 'drop'
  #   )
  # 
  # combined_plot <- ggplot() +
  #   # Add variance regions
  #   geom_rect(data = average_path,
  #             aes(xmin = mean_x - sd_x,
  #                 xmax = mean_x + sd_x,
  #                 ymin = mean_y - sd_y,
  #                 ymax = mean_y + sd_y,
  #                 fill = factor(marker_num)),
  #             alpha = 0.2) +
  #   # Add mean paths
  #   geom_path(data = average_path,
  #             aes(x = mean_x, y = mean_y, 
  #                 colour = factor(marker_num),
  #                 group = interaction(condition, marker_num))) +
  #   facet_grid(condition ~ speed, labeller = labeller(speed = speed_labels)) +
  #   scale_colour_viridis_d(name = "Marker") +
  #   scale_fill_viridis_d(name = "Marker") +
  #   theme_minimal() +
  #   labs(title = "Average Marker Paths with Variance",
  #        x = "X Position",
  #        y = "Y Position") +
  #   coord_fixed(ratio = 0.615) +
  #   scale_y_reverse(breaks = seq(160, 40, by = -40),
  #                   labels = seq(40, 160, by = 40)) +
  #   theme(
  #     plot.title = element_text(size = 13, hjust = 0.5),
  #     axis.title = element_text(size = 11),
  #     axis.text = element_text(size = 10),
  #     strip.text = element_text(size = 10),
  #     legend.title = element_text(size = 11),
  #     legend.text = element_text(size = 10),
  #     plot.background = element_rect(fill = "white", colour = NA),
  #     panel.background = element_rect(fill = "white", colour = NA)
  #   )
  # 
  # print(combined_plot)
}

# Combine plots
all_plots <-  (trajectory_plot) / (heatmap_plot) +
  plot_annotation(
    title = NULL,
    tag_levels = 'A',
    tag_prefix = '(',
    tag_suffix = ')',
    theme = theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"))
  )

print(all_plots)

# Save the plot
ggsave("gait_analysis/complete_analysis/marker_position_analysis.png", 
       all_plots, 
       width = 8, 
       height = 8, 
       dpi = 300)


#-------------------------------------------------------------------------------
#### TABLES FOR KINEMATIC AND TEMPORAL DATA
#-------------------------------------------------------------------------------

# For Kinematic Parameters
kinematic_summary <- all_kinematic_data %>%
  group_by(condition, speed) %>%
  summarise(
    knee_rom_mean = mean(knee_ROM, na.rm = TRUE),
    knee_rom_sd = sd(knee_ROM, na.rm = TRUE),
    ankle_rom_mean = mean(ankle_ROM, na.rm = TRUE),
    ankle_rom_sd = sd(ankle_ROM, na.rm = TRUE)
  )

kinematic_table <- kinematic_summary %>%
  kbl(caption = "Table 1. Kinematic Parameters Across Walking Speeds and Conditions",
      col.names = c("Condition", "Speed (km/h)", "Mean", "SD", "Mean", "SD"),
      format = "html",
      digits = 1) %>%
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE) %>%
  add_header_above(c(" " = 2, 
                     "Knee ROM (degrees)" = 2,
                     "Ankle ROM (degrees)" = 2))

print(kinematic_table)

# Save the table as an HTML file with wider styling
save_kable(kinematic_table, "gait_analysis/complete_analysis/table_kinematic.html", self_contained = TRUE)

# Open the HTML file and edit the <table> style
html_content <- readLines("gait_analysis/complete_analysis/table_kinematic.html")
html_content <- gsub("<table", "<table style='width:100%; border-collapse:collapse; font-family:Arial; font-size:11pt;'", html_content)
writeLines(html_content, "gait_analysis/complete_analysis/table_kinematic.html")

# Convert the updated HTML to JPG
webshot("gait_analysis/complete_analysis/table_kinematic.html", 
        "gait_analysis/complete_analysis/table_kinematic.jpg", 
        vwidth = 1200, vheight = 400, zoom = 2)

# For Spatiotemporal Parameters
spatiotemporal_summary <- all_temporal_data %>%
  group_by(condition, speed) %>%
  summarise(
    stance_mean = mean(mean_stance, na.rm = TRUE),
    stance_sd = sd(mean_stance, na.rm = TRUE),
    stride_length_mean = mean(mean_stride_length, na.rm = TRUE),
    stride_length_sd = sd(mean_stride_length, na.rm = TRUE),
    cadence_mean = mean(cadence, na.rm = TRUE),
    cadence_sd = sd(cadence, na.rm = TRUE)
  )

print(spatiotemporal_summary)

spatiotemporal_table <- spatiotemporal_summary %>%
  kbl(caption = "Table 2. Spatiotemporal Parameters Across Walking Speeds and Conditions",
      col.names = c("Condition", "Speed (km/h)", "Mean", "SD", "Mean", "SD", "Mean", "SD"),
      format = "html",
      digits = c(NA, 1, 1, 1, 3, 3, 1, 1)) %>%
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE) %>%
  add_header_above(c(" " = 2, 
                     "Stance (%)" = 2,
                     "Stride Length (m)" = 2,
                     "Cadence (steps/min)" = 2))

print(spatiotemporal_table)

# Save the table as an HTML file with wider styling
save_kable(spatiotemporal_table, "gait_analysis/complete_analysis/table_spatiotemporal.html", self_contained = TRUE)

# Open the HTML file and edit the <table> style
html_content <- readLines("gait_analysis/complete_analysis/table_spatiotemporal.html")
html_content <- gsub("<table", "<table style='width:100%; border-collapse:collapse; font-family:Arial; font-size:11pt;'", html_content)
writeLines(html_content, "gait_analysis/complete_analysis/table_spatiotemporal.html")

# Convert the HTML to JPG and save it in the complete_analysis folder
webshot("gait_analysis/complete_analysis/table_spatiotemporal.html", 
        "gait_analysis/complete_analysis/table_spatiotemporal.jpg", 
        vwidth = 1200, vheight = 400, zoom = 2)



#-------------------------------------------------------------------------------
#### CLASSIFICATION ANALYSIS
#-------------------------------------------------------------------------------

library(randomForest)
library(caret)

# Prepare combined dataset for classification
classification_data <- all_kinematic_data %>%
  left_join(all_temporal_data, by = c("subject_id", "condition", "speed")) %>%
  select(
    subject_id, condition, speed,
    knee_ROM, ankle_ROM, mean_knee_angle, mean_ankle_angle,
    mean_stance, mean_stride_length, cadence
  )

# Prepare features and target
features <- classification_data %>%
  select(-subject_id, -condition, -speed) %>%
  scale()

# Convert to matrix/dataframe format (RandomForest expects this!!)
features <- as.data.frame(features)
target <- factor(classification_data$condition)

# Train basic random forest first
rf_basic <- randomForest(x = features, 
                         y = target, 
                         ntree = 500,
                         importance = TRUE)

print(rf_basic)

# Get feature importance directly from randomForest
importance_scores <- importance(rf_basic)
importance_df <- data.frame(
  Feature = rownames(importance_scores),
  Importance = importance_scores[, "MeanDecreaseGini"]
) %>%
  arrange(desc(Importance))

print(importance_df)

# Calculate confusion matrix
predictions <- predict(rf_basic, features)
conf_matrix <- confusionMatrix(predictions, target)

print(conf_matrix)

# Create feature importance plot
feature_importance_plot <- ggplot(importance_df, 
                                  aes(x = reorder(Feature, Importance), 
                                      y = Importance)) +
  geom_bar(stat = "identity", fill = "#4682B4") +
  coord_flip() +
  theme_minimal() +
  labs(x = "Feature", 
       y = "Importance Score",
       title = "Feature Importance in Gait Classification")

print(feature_importance_plot)
ggsave("feature_importance.png", feature_importance_plot, width = 8, height = 6)

# Cross-validation
if(exists("rf_basic")) {
  # Set up cross-validation
  set.seed(42)
  train_control <- trainControl(
    method = "cv",
    number = 4,  # 4-fold CV for the 4 subjects
    classProbs = TRUE
  )
  
  # Train model with cross-validation
  rf_model <- train(
    x = features,
    y = target,
    method = "rf",
    trControl = train_control,
    importance = TRUE
  )
  
  print("\nCross-validated Model Performance:")
  print(rf_model)
}


#### PLOTS

# A. Feature Importance Plot (refined version)
importance_plot <- ggplot(importance_df, 
                          aes(x = reorder(Feature, Importance), 
                              y = Importance)) +
  geom_bar(stat = "identity", fill = "#377EB8") +
  coord_flip() +
  theme_minimal() +
  labs(x = "Gait Parameter",  
       y = "Relative Importance Score (%)") +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# B. Confusion Matrix Heatmap
# Get the raw confusion matrix
conf_raw <- conf_matrix$table
# Convert to df
conf_data <- as.data.frame(conf_raw)
names(conf_data) <- c("Actual", "Predicted", "Count")

confusion_plot <- ggplot(conf_data, aes(x = Actual, y = Predicted, fill = Count)) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = "white", high = "#377EB8") +
  geom_text(aes(label = Count), colour = "black", size = 4) +
  theme_minimal() +
  labs(x = "Actual Condition",
       y = "Predicted Condition") +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# Combine plots with lettered labels
combined_plot <- (importance_plot + confusion_plot) +
  plot_annotation(
    title = "Gait Classification Analysis",
    tag_levels = 'A',
    tag_prefix = '(',
    tag_suffix = ')',
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.tag = element_text(size = 12, face = "bold")
    )
  )

print(combined_plot)

ggsave("gait_analysis/complete_analysis/classification_analysis.png", 
       plot = combined_plot,
       width = 8,
       height = 5,
       dpi = 300)

print(conf_matrix$overall)
print(conf_matrix$byClass)


#### LDA Classification Analysis

# Prepare data more explicitly
classification_data <- all_kinematic_data %>%
  dplyr::left_join(all_temporal_data, by = c("subject_id", "condition", "speed")) %>%
  dplyr::select(dplyr::all_of(c(
    "subject_id", "condition", "speed",
    "knee_ROM", "ankle_ROM", "mean_knee_angle", "mean_ankle_angle",
    "mean_stance", "mean_stride_length", "cadence"
  )))

str(classification_data)

# Prepare features and scale
features <- classification_data %>%
  dplyr::select(-subject_id, -condition, -speed) %>%
  scale() %>%
  as.data.frame()

# Fit LDA model
lda_model <- lda(features, classification_data$condition)

# Get predictions
predictions <- predict(lda_model, features)
predicted_classes <- predictions$class

# Create confusion matrix
conf_matrix <- table(Actual = classification_data$condition, 
                     Predicted = predicted_classes)

# Convert to df for plotting
conf_data <- as.data.frame(conf_matrix)
names(conf_data) <- c("Actual", "Predicted", "Count")

# Calculate feature importance from LDA coefficients
importance_scores <- abs(lda_model$scaling) %>%
  as.data.frame() %>%
  mutate(Feature = rownames(.)) %>%
  gather(key = "LD", value = "Score", -Feature) %>%
  group_by(Feature) %>%
  summarise(Importance = mean(abs(Score))) %>%
  arrange(desc(Importance))

# Plots
# A. Feature Importance
importance_plot <- ggplot(importance_scores, 
                          aes(x = reorder(Feature, Importance), 
                              y = Importance)) +
  geom_bar(stat = "identity", fill = "#377EB8") +
  coord_flip() +
  theme_minimal() +
  labs(x = "Gait Parameter",
       y = "LDA Coefficient Magnitude") +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# B. Confusion Matrix Plot
confusion_plot <- ggplot(conf_data, aes(x = Actual, y = Predicted, fill = Count)) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = "white", high = "#377EB8") +
  geom_text(aes(label = Count), colour = "black", size = 4) +
  theme_minimal() +
  labs(x = "Actual Condition",
       y = "Predicted Condition") +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# C. LDA Plot
lda_scores <- as.data.frame(predictions$x)
lda_scores$condition <- classification_data$condition
lda_scores$speed <- classification_data$speed

lda_plot <- ggplot(lda_scores, 
                   aes(x = LD1, y = LD2, 
                       colour = condition, 
                       shape = speed)) +
  geom_point(size = 3) +
  stat_ellipse(aes(colour = condition), level = 0.95) +
  theme_minimal() +
  labs(x = "First Linear Discriminant",
       y = "Second Linear Discriminant",
       colour = "Condition",
       shape = "Speed (km/h)") +
  scale_colour_manual(values = c("normal" = "#377EB8", 
                                "limp" = "#4DAF4A", 
                                "PD" = "#E41A1C")) +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

# Combine plots
combined_plot <- (importance_plot + confusion_plot) / lda_plot +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title = "Gait Classification Analysis Using LDA",
    tag_levels = 'A',
    tag_prefix = '(',
    tag_suffix = ')',
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.tag = element_text(size = 12, face = "bold")
    )
  )

print(combined_plot)

# Save plot
ggsave("gait_analysis/complete_analysis/lda_classification_analysis.png", 
       plot = combined_plot,
       width = 10,
       height = 8,
       dpi = 300)

# Print accuracy metrics
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
print(paste("Overall accuracy:", round(accuracy * 100, 1), "%"))

print(conf_matrix)


##### LDA ANALYSIS - RESULTS

# Print overall accuracy
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
print(paste("Overall accuracy:", round(accuracy * 100, 1), "%"))

print(conf_matrix)

# Calculate per-class metrics
class_metrics <- data.frame(
  Condition = rownames(conf_matrix),
  Sensitivity = diag(prop.table(conf_matrix, 1)),
  Precision = diag(prop.table(conf_matrix, 2))
) %>%
  mutate(
    Sensitivity = round(Sensitivity * 100, 1),
    Precision = round(Precision * 100, 1)
  )

print(class_metrics)

# Print LDA coefficients for each discriminant
print(round(lda_model$scaling, 3))

# Print proportion of trace (importance of each discriminant)
print(round(prop.table(lda_model$svd^2), 3))


#### LDA Cross Validation

# Prepare df for training more explicitly
training_data <- classification_data %>%
  dplyr::select(dplyr::all_of(c("knee_ROM", "ankle_ROM", "mean_knee_angle", "mean_ankle_angle",
                                "mean_stance", "mean_stride_length", "cadence", "condition"))) %>%
  mutate(across(-condition, scale))

str(training_data)

# Set up cross-validation
set.seed(42)
train_control <- trainControl(
  method = "cv",
  number = 4,  # 4-fold CV because we have 4 subjects
  savePredictions = TRUE,
  classProbs = TRUE
)

# Train LDA with cross-validation
lda_cv <- train(
  condition ~ .,
  data = training_data,
  method = "lda",
  trControl = train_control
)

# Get cross-validated predictions
cv_preds <- lda_cv$pred

# Calculate ROC curves
roc_curves <- list()
auc_values <- numeric(3)
conditions <- levels(classification_data$condition)

# Calculate ROC and AUC for each condition
for(i in seq_along(conditions)) {
  # Binary outcomes for each class
  binary_outcome <- ifelse(classification_data$condition == conditions[i], 1, 0)
  
  # Probabilities for this class
  class_probs <- cv_preds[, conditions[i]]
  
  # ROC curve
  roc_curves[[i]] <- roc(binary_outcome, class_probs)
  auc_values[i] <- auc(roc_curves[[i]])
}

print(lda_cv)

# AUC values for each class
print("\nAUC Values:")
for(i in seq_along(conditions)) {
  print(paste(conditions[i], "AUC:", round(auc_values[i], 3)))
}

# Create ROC plot
roc_plot <- ggroc(roc_curves, legacy.axes = TRUE) +
  geom_line(size = 1) +
  scale_colour_manual(
    values = c("#377EB8", "#4DAF4A", "#E41A1C"),
    labels = paste(conditions, "AUC:", round(auc_values, 3))
  ) +
  theme_minimal() +
  labs(
    title = "ROC Curves for Gait Classification",
    x = "False Positive Rate",
    y = "True Positive Rate",
    colour = "Condition"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

print(roc_plot)

# Calculate confusion matrix from cross-validated predictions
cv_conf_matrix <- confusionMatrix(cv_preds$pred, cv_preds$obs)

print(cv_conf_matrix)

# Save the ROC plot
ggsave("gait_analysis/complete_analysis/roc_curves.png", 
       plot = roc_plot,
       width = 8,
       height = 6,
       dpi = 300)