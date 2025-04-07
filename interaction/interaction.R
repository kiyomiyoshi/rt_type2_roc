#+ message = F, warning = F
library(tidyverse)
library(gghalves)
library(GGally)
library(car)
library(doParallel)
library(sjPlot)
library(cowplot)
theme_set(theme_publication()) 

script <- c("Maniscalco_2017_expt1.R",
            "Maniscalco_2017_expt3.R",
            "Maniscalco_2017_expt4.R",
            "Hainguerlot_2018.R", 
            "Hainguerlot_unpub.R",
            "Reyes_2015.R")

exp <- c("Maniscalco_2017_expt1",
         "Maniscalco_2017_expt3",
         "Maniscalco_2017_expt4",
         "Hainguerlot_2018", 
         "Hainguerlot_unpub",
         "Reyes_2015")

df_grand <- c()
df_accuracy <- c()
df_logit <- c()
sl_accuracy <- c()
sl_logit <- c()

for (j in 1:length(script)) {
  source(script[j]) 
  df$exp <- exp[j]
  df_grand <- rbind(df_grand, df)
  
  df1$exp <- exp[j]
  df2$exp <- exp[j]
  df_accuracy <- rbind(df_accuracy, df1) 
  df_logit <- rbind(df_logit, df2) 
  
  slope_1$exp <- exp[j]
  slope_2$exp <- exp[j]
  sl_accuracy <- rbind(sl_accuracy, slope_1) 
  sl_logit <- rbind(sl_logit, slope_2) 
}  

# visualization
levels(df_accuracy$RT_bin) <- c("q1", "q2", "q3", "q4")
df_accuracy$exp <- factor(df_accuracy$exp, levels = exp)
ggplot(df_accuracy) + 
  geom_line(aes(x = Confidence, y = mean_accuracy, color = RT_bin, group = RT_bin), size = 0.3, 
            position = position_dodge(width = 0.2)) +
  geom_point(aes(x = Confidence, y = mean_accuracy, color = RT_bin), size = 0.7, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(x = Confidence, ymin = mean_accuracy - se, ymax = mean_accuracy + se, color = RT_bin), 
                size = 0.3, width = 0.1, position = position_dodge(width = 0.2)) +
  scale_color_manual(values = c("#3300CC", "#7900FF", "#FF0077", "#CC0033")) +
  theme(legend.position = "top") + labs(color = "RT bin") +
  xlab("Confidence") + ylab("Accuracy") + coord_cartesian(ylim = c(0.4, 1)) +
  facet_wrap(.~ exp) +
  theme(axis.title = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(0.5, 0.5, 0.5, 0.5, "pt")),
        plot.margin = unit(c(1, 1, 1, 1), "pt")) -> p1
p1
sjPlot::save_plot("figure_s8.jpg", p1, width = 12, height = 8, dpi = 300)

levels(df_logit$RT_bin) <- c("q1", "q2", "q3", "q4")
df_logit$exp <- factor(df_logit$exp, levels = exp)
ggplot(df_logit) + 
  geom_line(aes(x = Confidence, y = mean_logit, color = RT_bin, group = RT_bin), size = 0.3, 
            position = position_dodge(width = 0.2)) +
  geom_point(aes(x = Confidence, y = mean_logit, color = RT_bin), size = 0.7, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(x = Confidence, ymin = mean_logit - se, ymax = mean_logit + se, color = RT_bin),
                size = 0.3, width = 0.1, position = position_dodge(width = 0.2)) +
  scale_color_manual(values = c("#3300CC", "#7900FF", "#FF0077", "#CC0033")) +
  theme(legend.position = "top") + labs(color = "RT bin") +
  xlab("Confidence") + ylab("Logit-transformed accuracy") + coord_cartesian(ylim = c(-0.5, 3.5)) +
  facet_wrap(.~ exp) +
  theme(axis.title = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(0.5, 0.5, 0.5, 0.5, "pt")),
        plot.margin = unit(c(1, 1, 1, 1), "pt")) -> p2
p2
sjPlot::save_plot("figure_s9.jpg", p2, width = 12, height = 8, dpi = 300)

plot_grid(p1, p2, nrow = 2)


df_grand %>%                   
  group_by(RT_bin, exp, sub) %>%
  summarise(sd = sd(accuracy)) %>%
  ungroup(sub) %>%
  summarise(mean_sd_expwise = sd(sd))



df_accuracy %>%                              # U シェイプ
  group_by(RT_bin, exp) %>%                  # 確信度1も「全くわからんという確信が強い（即決）」ケースと
  summarise(sd_conf = sd(mean_accuracy)) %>% #「悩みまくって結局何の手掛かりもつかめなかった（遅い）」ケースがあるのでは
  ungroup(exp) %>%
  summarise(mean_sd_conf = sd(sd_conf))

df_accuracy %>%                            　# 確信度が高くて既にaccuracy帯がすでにかなり高いにもかかわらず、    
  group_by(Confidence, exp) %>%              # そういうときほどRTによる正答率の違いは大きい（RTが反応正誤について持っている情報は大きい）
  summarise(sd_RT = sd(mean_accuracy)) %>%   # じっくり考えた場合も即決だった場合も主観的な確信度は同じだが、即決のほうが正答率はずっととかい
  ungroup(exp) %>%                           # 確信度4でも反応遅い(q4)やつは他3つから思いっきり離されてる
  summarise(mean_sd_RT = sd(sd_RT))          # 自覚的でない成分をとらえている

df_logit %>%
  group_by(RT_bin, exp) %>%
  summarise(sd_conf = sd(mean_logit)) %>%
  ungroup(exp) %>%
  summarise(mean_sd_conf = sd(sd_conf))

df_logit %>%
  group_by(Confidence, exp) %>%
  summarise(sd_RT = sd(mean_logit)) %>%
  ungroup(exp) %>%
  summarise(mean_sd_RT = sd(sd_RT))

# pairwise slope
sl_accuracy %>%
  na.omit() %>%
  group_by(RT_bin, pair) %>%
  summarise(mean_slope = mean(slope),
            sd = sd(slope),
            se = sd / sqrt(length(exp))) %>%
  ggplot(aes(x = pair, y = mean_slope, fill = RT_bin)) + 
  geom_bar(aes(x = pair, y = mean_slope, fill = RT_bin), 
           size = 0.7, stat = "identity", position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = mean_slope - se, ymax = mean_slope + se),
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#3300CC", "#7900FF", "#FF0077", "#CC0033")) +
  theme(legend.position = "top") + labs(color = "RT bin") +
  xlab("Pair") + ylab("Slope of accuracy") + coord_cartesian(ylim = c(0, 0.2)) +
  theme(axis.title = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(-0.5, -0.5, -0.5, -0.5, "pt")),
        plot.margin = unit(c(1, 1, 1, 1), "pt")) -> p3

sl_logit %>%
  na.omit() %>%
  group_by(RT_bin, pair) %>%
  summarise(mean_slope = mean(slope),
            sd = sd(slope),
            se = sd / sqrt(length(exp))) %>%
  ggplot(aes(x = pair, y = mean_slope, fill = RT_bin)) + 
  geom_bar(aes(x = pair, y = mean_slope, fill = RT_bin), 
           size = 0.7, stat = "identity", position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = mean_slope - se, ymax = mean_slope + se),
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#3300CC", "#7900FF", "#FF0077", "#CC0033")) +
  theme(legend.position = "top") + labs(color = "RT bin") +
  xlab("Pair") + ylab("Slope of logit") + coord_cartesian(ylim = c(0, 1.5)) +
  theme(axis.title = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(-0.5, -0.5, -0.5, -0.5, "pt")),
        plot.margin = unit(c(1, 1, 1, 1), "pt")) -> p4

plot_grid(p3, p4)
save_plot("accuracy.jpg", p1)
save_plot("logit.jpg", p2)

df_accuracy$text1 <- c("β[interaction] == -0.139", "β[interaction] == -0.317", "β[interaction] == -0.448",
                       "β[interaction] == -0.234", "β[interaction] == -0.138", "β[interaction] == -0.356")[df_accuracy$exp]

df_accuracy$text2 <- c("discrete: 1, 2.67, 4.34, 6", "discrete: 1, 2.67, 4.34, 6", "discrete: 1, 2.67, 4.34, 6",
                       "discrete: 1, 2, 3, 4, 5, 6", "discrete: 1, 2, 3, 4, 5, 6", "continuous: 1-6")[df_accuracy$exp]


ggplot(df_accuracy) + 
  geom_line(aes(x = Confidence, y = mean_accuracy, color = RT_bin, group = RT_bin), size = 0.3, 
            position = position_dodge(width = 0.2)) +
  geom_point(aes(x = Confidence, y = mean_accuracy, color = RT_bin), size = 0.7, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(x = Confidence, ymin = mean_accuracy - se, ymax = mean_accuracy + se, color = RT_bin), 
                size = 0.3, width = 0.1, position = position_dodge(width = 0.2)) +
  scale_color_manual(values = c("#3300CC", "#7900FF", "#FF0077", "#CC0033")) +
  theme(legend.position = "top") + labs(color = "RT bin") +
  xlab("Confidence") + ylab("Accuracy") + coord_cartesian(ylim = c(0.4, 1)) +
  facet_wrap(.~ exp) +
  geom_text(aes(x = 2.5, y = 0.5, label = text1), size = 2, parse = TRUE) +
  geom_text(aes(x = 2.5, y = 0.45, label = text2), size = 2) +
  theme(axis.title = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(0.5, 0.5, 0.5, 0.5, "pt")),
        plot.margin = unit(c(1, 1, 1, 1), "pt")) -> p1
p1