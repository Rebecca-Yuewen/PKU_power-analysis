library(dplyr)
require(nlme)

setwd("D:/PKU/生态保护科学研究设计/power analysis file/paper practice")
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

###############demo#########################
df1<-data.frame(block=rep(1:7,4),Treatment=rep(c("CONTROL", "LIGHT", "MODERATE", "INTENSIVE"),each=7))
df1$Treatment <- factor(df1$Treatment, levels = c("CONTROL", "LIGHT", "MODERATE", "INTENSIVE"))
df2 <- left_join(df1, df_pol[, c("Treatment", "pred_sr", "sd_sr")], by = "Treatment")

df2 <- df2 %>%
  group_by(Treatment) %>%
  mutate(richness = rnorm(n(), mean = pred_sr, sd = sd_sr)) %>%  #n() gives the current group size
  ungroup() %>%
  mutate(block_effect = rnorm(n = length(unique(df2$block)), mean = 0, sd = 3.54)[df2$block]) %>%
  mutate(richness_new = richness + block_effect)

model_pol<-lme(richness_new ~ Treatment, random = ~1 | block, data = df2)
summary(model_pol)$tTable[, "p-value"]

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
power_pol(sample_size = 5, block_sd = 3.54, alpha = 0.05)

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
power_moth(sample_size = 471, block_sd = 5.12, alpha = 0.05)

######################test when the power reach 0.8###################
###pollinator###
power_levels_pol <- list(LIGHT = c(), MODERATE = c(), INTENSIVE = c())
sample_sizes <- c(2,3,4,6,8,10,15,20,30,40)
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
sample_sizes1 <- c(2,3,4,6,8,10,15,20,30,40,50,100,200,300,400,500,600)
power_levels_moth2 <- list(LIGHT = c(), MODERATE = c(), INTENSIVE = c())
sample_sizes2 <- c(610,620,630)
for (i in 1:length(sample_sizes2)) {
  power_levels_moth2$LIGHT[i] <- power_moth(sample_size = sample_sizes2[i], block_sd = 5.12, alpha = 0.05)$LIGHT
  power_levels_moth2$MODERATE[i] <- power_moth(sample_size = sample_sizes2[i], block_sd = 5.12, alpha = 0.05)$MODERATE
  power_levels_moth2$INTENSIVE[i] <- power_moth(sample_size = sample_sizes2[i], block_sd = 5.12, alpha = 0.05)$INTENSIVE
}
power_levels_moth2
power_df_moth <- data.frame(sample_size = sample_sizes,
                            Light = power_levels_moth$LIGHT, Moderate = power_levels_moth$MODERATE, Intensive = power_levels_moth$INTENSIVE)
power_df_moth2 <- data.frame(sample_size = sample_sizes2,
                            Light = power_levels_moth2$LIGHT, Moderate = power_levels_moth2$MODERATE, Intensive = power_levels_moth2$INTENSIVE)

write.csv(power_df_moth,"moth power test.csv")


# plot curve
library(ggplot2)
library(tidyr)
power_df_pol<-read.csv("pollinator power test.csv")
power_df_moth<-read.csv("moth power test.csv")
# 将数据框转换为长格式，便于 ggplot 使用
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

l_pol<-power_df_pol[,2:3]
m_pol<-power_df_pol[,c(2,4)]
i_pol<-power_df_pol[,c(2,5)]

l_moth<-power_df_moth[,2:3]
m_moth<-power_df_moth[,c(2,4)]
i_moth<-power_df_moth[,c(2,5)]

install.packages("patchwork")
library(patchwork)

p1<-ggplot(l_pol, aes(x = sample_size, y = Light)) +
  geom_line(linewidth = 1.2, alpha=0.7, color = "blue4") +
  geom_point(size = 1.7) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgrey", size = 1) +
  labs(x = "Sample Size",y = "Power") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 20)) 
p2<-ggplot(m_pol, aes(x = sample_size, y = Moderate)) +
  geom_line(linewidth = 1.2, alpha=0.7, color = "goldenrod1") +
  geom_point(size = 1.7,color = "goldenrod1") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgrey", size = 1) +
  labs(x = "Sample Size",y = "Power") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 8)) 
p3<-ggplot(i_pol, aes(x = sample_size, y = Intensive)) +
  geom_line(linewidth = 1.2, alpha=0.7, color = "red3") +
  geom_point(size = 1.7, color = "red3") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgrey", size = 1) +
  labs(x = "Sample Size",y = "Power") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 8)) 
p4<-ggplot(l_moth, aes(x = sample_size, y = Light)) +
  geom_line(linewidth = 1.2, alpha=0.7, color = "blue4") +
  geom_point(size = 1.7) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgrey", size = 1) +
  labs(x = "Sample Size",y = "Power") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 600)) 
p5<-ggplot(m_moth, aes(x = sample_size, y = Moderate)) +
  geom_line(linewidth = 1.2, alpha=0.7, color = "goldenrod1") +
  geom_point(size = 1.7,color = "goldenrod1") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgrey", size = 1) +
  labs(x = "Sample Size",y = "Power") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 8)) 
p6<-ggplot(i_moth, aes(x = sample_size, y = Intensive)) +
  geom_line(linewidth = 1.2, alpha=0.7, color = "red3") +
  geom_point(size = 1.7, color = "red3") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgrey", size = 1) +
  labs(x = "Sample Size",y = "Power") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 20)) 

combined_plot <- (p1 + p2 + p3) / (p4 + p5 + p6)

# 导出为PNG图片
ggsave("combined_plot.png", combined_plot, width = 10, height = 6, dpi = 300)
