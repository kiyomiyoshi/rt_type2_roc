library(tidyverse)
library(data.table)
library(doParallel)

cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

dat <- fread("data_Hainguerlot_unpub.csv", header = T)

df <- foreach(i = 1:length(unique(dat$Subj_idx)), .combine = "rbind", .packages = c("dplyr", "metaSDT")) %dopar% {
  
  d <- subset(dat, dat$Subj_idx == unique(dat$Subj_idx)[i])
  d <- na.omit(d)
  d <- mutate(d, Correct = ifelse(Stimulus == Response, 1, 0))
  d <- mutate(d, Confidence = ifelse(Confidence < 0.7, 1,
                                     ifelse(Confidence < 0.9, 2, 3)))
  
  # RT bin
  q <- quantile(d$RT_dec, probs = seq(0, 1, by = 0.25), na.rm = TRUE)
  d$RT_bin <- cut(d$RT_dec, breaks = q, labels = FALSE, include.lowest = TRUE)
  d$RT_bin <- as.factor(d$RT_bin)
  d %>%
    group_by(Confidence, RT_bin) %>%
    summarise(accuracy = mean(Correct)) %>%
    mutate(logit = ifelse(accuracy == 1, qlogis(accuracy - 0.02),
                          ifelse(accuracy == 0, qlogis(accuracy + 0.02), qlogis(accuracy)))) %>%
    mutate(sub = i) %>%
    print()
  
}

df %>%
  group_by(Confidence, RT_bin) %>%
  summarise(mean_accuracy = mean(accuracy),
            sd = sd(accuracy),
            se = sd/sqrt(length(unique(dat$Subj_idx)))) -> df1
ggplot(df1) + 
  geom_line(aes(x = Confidence, y = mean_accuracy, color = RT_bin, group = RT_bin), size = 0.3) +
  geom_point(aes(x = Confidence, y = mean_accuracy, color = RT_bin), size = 0.7) +
  geom_errorbar(aes(x = Confidence, ymin = mean_accuracy - se, ymax = mean_accuracy + se, color = RT_bin), width = 0.2) +
  scale_color_manual(values = c("#3300CC", "#7900FF", "#FF0077", "#CC0033")) +
  theme(legend.position = "top") + labs(color = "RT bin") +
  xlab("Confidence") + ylab("Accuracy") + coord_cartesian(ylim = c(0.5, 1)) +
  theme(axis.title = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(-0.5, -0.5, -0.5, -0.5, "pt")),
        plot.margin = unit(c(1, 1, 1, 1), "pt")) + ggtitle("Hainguerlot_unpub") -> p1

df %>%
  group_by(Confidence, RT_bin) %>%
  summarise(mean_logit = mean(logit),
            sd = sd(logit),
            se = sd/sqrt(length(unique(dat$Subj_idx)))) -> df2
ggplot(df2) + 
  geom_line(aes(x = Confidence, y = mean_logit, color = RT_bin, group = RT_bin), size = 0.3) +
  geom_point(aes(x = Confidence, y = mean_logit, color = RT_bin), size = 0.7) +
  geom_errorbar(aes(x = Confidence, ymin = mean_logit - se, ymax = mean_logit + se, color = RT_bin), width = 0.2) +
  scale_color_manual(values = c("#3300CC", "#7900FF", "#FF0077", "#CC0033")) +
  theme(legend.position = "top") + labs(color = "RT bin") +
  xlab("Confidence") + ylab("Log odds ratio of correct response") + coord_cartesian(ylim = c(-0.5, 3.5)) +
  theme(axis.title = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(-0.5, -0.5, -0.5, -0.5, "pt")),
        plot.margin = unit(c(1, 1, 1, 1), "pt")) + ggtitle("Hainguerlot_unpub") -> p2

df1 %>%
  arrange(RT_bin, Confidence) %>%
  group_by(RT_bin) %>%
  mutate(slope = lead(mean_accuracy, 1) - mean_accuracy,
         pair = rep(c("2-1", "3-2", "NA"))) -> slope_1

df2 %>%
  arrange(RT_bin, Confidence) %>%
  group_by(RT_bin) %>%
  mutate(slope = lead(mean_logit, 1) - mean_logit,
         pair = rep(c("2-1", "3-2", "NA"))) -> slope_2