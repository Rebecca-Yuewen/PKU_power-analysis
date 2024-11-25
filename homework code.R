library(dplyr)
require(nlme)

data <- read.csv("overview_herbicide_biodiv_effects.csv")

df_pol<-data %>% 
  subset(taxon %in% "pollinators") %>% 
  select(Treatment, pred_sr, lower_sr, upper_sr, taxon) %>%
  mutate(sd_sr = ((upper_sr-lower_sr)/(2 * 1.96)))
df_pol$Treatment <- factor(df_pol$Treatment, levels = c("CONTROL", "LIGHT", "MODERATE", "INTENSIVE"))

df_moth<-data %>% 
  subset(taxon %in% "moths") %>% 
  select(Treatment, pred_sr, lower_sr, upper_sr, taxon) %>%
  mutate(sd_sr = ((upper_sr-lower_sr)/(2 * 1.96)))
df_moth$Treatment <- factor(df_moth$Treatment, levels = c("CONTROL", "LIGHT", "MODERATE", "INTENSIVE"))

##############simulation 1000 times######################
power_pol <- function(sample_size, block_sd, alpha){
  set.seed(110)
#construct data frame 
  df1 <- data.frame(block = rep(1:sample_size, 4),
                    Treatment = rep(c("CONTROL", "LIGHT", "MODERATE", "INTENSIVE"), each = sample_size))
  df1$Treatment <- factor(df1$Treatment, levels = c("CONTROL", "LIGHT", "MODERATE", "INTENSIVE"))
  df2 <- left_join(df1, df_pol[, c("Treatment", "pred_sr", "sd_sr")], by = "Treatment")
  light_results <- c()
  moderate_results <- c()
  intensive_results <- c() 
for (i in 1:1000) {
  # simulation
  df2_sim <- df2 %>%
    group_by(Treatment) %>%
    mutate(richness = rnorm(n(), mean = pred_sr, sd = sd_sr)) %>%
    ungroup() %>%
    mutate(block_effect = rnorm(n = length(unique(df2$block)), mean = 0, sd = block_sd)[df2$block]) %>%
    mutate(richness_new = richness + block_effect)
  # construct the model
  model <- lme(richness_new ~ Treatment, random = ~1 | block, data = df2_sim)
  # get the results
  p_values <- summary(model)$tTable[, "p-value"][-1]  # overlook the intercept p value
  light_results[i] <- p_values[1]<= alpha
  moderate_results[i] <- p_values[2]<= alpha
  intensive_results[i] <- p_values[3]<= alpha
}
  list(
    LIGHT = mean(light_results),
    MODERATE = mean(moderate_results),
    INTENSIVE = mean(intensive_results)
  )
}
power_pol(sample_size = 7, block_sd = 3.54, alpha = 0.05)

power_moth <- function(sample_size, block_sd, alpha){
  set.seed(120)
  #construct data frame 
  df1 <- data.frame(block = rep(1:sample_size, 4),
                    Treatment = rep(c("CONTROL", "LIGHT", "MODERATE", "INTENSIVE"), each = sample_size))
  df1$Treatment <- factor(df1$Treatment, levels = c("CONTROL", "LIGHT", "MODERATE", "INTENSIVE"))
  df2 <- left_join(df1, df_moth[, c("Treatment", "pred_sr", "sd_sr")], by = "Treatment")
  light_results <- c()
  moderate_results <- c()
  intensive_results <- c()  
  for (i in 1:1000) {
    # simulation
    df2_sim <- df2 %>%
      group_by(Treatment) %>%
      mutate(richness = rnorm(n(), mean = pred_sr, sd = sd_sr)) %>%
      ungroup() %>%
      mutate(block_effect = rnorm(n = length(unique(df2$block)), mean = 0, sd = block_sd)[df2$block]) %>%
      mutate(richness_new = richness + block_effect)
    # construct the model
    model <- lme(richness_new ~ Treatment, random = ~1 | block, data = df2_sim)
    # get the results
    p_values <- summary(model)$tTable[, "p-value"][-1]  # overlook the intercept p value
    light_results[i] <- p_values[1]<= alpha
    moderate_results[i] <- p_values[2]<= alpha
    intensive_results[i] <- p_values[3]<= alpha
  }
  list(
    LIGHT = mean(light_results),
    MODERATE = mean(moderate_results),
    INTENSIVE = mean(intensive_results)
  )
}
power_moth(sample_size = 8, block_sd = 5.12, alpha = 0.05)

######################test when the power reach 0.8###################
###pollinator###
power_levels_pol <- list(LIGHT = c(), MODERATE = c(), INTENSIVE = c())
sample_sizes <- c(2,3,4,6,8,10,15,20,30,40,50,100,200,300,400,500,600)
for (i in 1:length(sample_sizes)) {
  power_levels_pol$LIGHT[i] <- power_pol(sample_size = sample_sizes[i], block_sd = 3.54, alpha = 0.05)$LIGHT
  power_levels_pol$MODERATE[i] <- power_pol(sample_size = sample_sizes[i], block_sd = 3.54, alpha = 0.05)$MODERATE
  power_levels_pol$INTENSIVE[i] <- power_pol(sample_size = sample_sizes[i], block_sd = 3.54, alpha = 0.05)$INTENSIVE
}
power_levels_pol
power_df_pol <- data.frame(sample_size = sample_sizes,
  Light = power_levels_pol$LIGHT, Moderate = power_levels_pol$MODERATE, Intensive = power_levels_pol$INTENSIVE)
write.csv(power_df_pol,"pollinator power test.csv")
###moth###
power_levels_moth <- list(LIGHT = c(), MODERATE = c(), INTENSIVE = c())
for (i in 1:length(sample_sizes)) {
  power_levels_moth$LIGHT[i] <- power_moth(sample_size = sample_sizes[i], block_sd = 5.12, alpha = 0.05)$LIGHT
  power_levels_moth$MODERATE[i] <- power_moth(sample_size = sample_sizes[i], block_sd = 5.12, alpha = 0.05)$MODERATE
  power_levels_moth$INTENSIVE[i] <- power_moth(sample_size = sample_sizes[i], block_sd = 5.12, alpha = 0.05)$INTENSIVE
}
power_levels_moth
power_df_moth <- data.frame(sample_size = sample_sizes,
                            Light = power_levels_moth$LIGHT, Moderate = power_levels_moth$MODERATE, Intensive = power_levels_moth$INTENSIVE)
write.csv(power_df_moth,"moth power test.csv")

# plot curve
library(ggplot2)
library(tidyr)
power_long_pol <- pivot_longer(power_df_pol, cols = c("Light", "Moderate", "Intensive"), 
                               names_to = "Treatment", values_to = "Power")
power_long_moth <- pivot_longer(power_df_moth, cols = c("Light", "Moderate", "Intensive"), 
                                names_to = "Treatment", values_to = "Power")

plot_pol<-ggplot(power_long_pol, aes(x = sample_size, y = Power, color = Treatment, group = Treatment)) +
  geom_line(linewidth = 1.2, alpha=0.7) +
  geom_point(size = 1.7) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgrey", size = 1) +
  labs(x = "Sample Size",y = "Power") +
  theme_bw() + 
  scale_color_manual(values = c("Light" = "blue4", "Moderate" = "goldenrod1", "Intensive" = "red3")) +
  scale_x_continuous(limits = c(0, 50)) +
  theme(legend.title = element_blank(), legend.position = c(0.9, 0.82),
        axis.title = element_text(size = 14), axis.text = element_text(size = 14), 
        legend.text = element_text(size = 14), legend.background = element_blank())
ggsave("power_pol_plot.png", plot = plot_pol, width = 8, height = 7, dpi = 300)

plot_moth<-ggplot(power_long_moth, aes(x = sample_size, y = Power, color = Treatment, group = Treatment)) +
  geom_line(linewidth = 1.2, alpha=0.7) +
  geom_point(size = 1.7) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgrey", size = 1) +
  labs(x = "Sample Size",y = "Power") +
  theme_bw() + 
  scale_color_manual(values = c("Light" = "blue4", "Moderate" = "goldenrod1", "Intensive" = "red3")) +
  scale_x_continuous(limits = c(0, 700)) +
  theme(legend.title = element_blank(), legend.position = c(0.91, 0.86),
        axis.title = element_text(size = 14), axis.text = element_text(size = 14), 
        legend.text = element_text(size = 14), legend.background = element_blank())
ggsave("power_moth_plot.png", plot = plot_moth, width = 8, height = 7, dpi = 300)
