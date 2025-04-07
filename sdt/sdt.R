#+ message = F, warning = F
library(tidyverse)
library(gghalves)
library(GGally)
library(car)
library(doParallel)
library(sjPlot)
theme_set(theme_publication()) 

script <- c("Hainguerlot_2018.R",        # too fast or slow response discouraged
            "Hainguerlot_unpub.R",       # too fast or slow response discouraged
            "Maniscalco_2017_expt1.R",   # trial pacing depended on subject RT, 5s time limit
            "Maniscalco_2017_expt2_1.R", # trial pacing depended on subject RT, 5s time limit
            "Maniscalco_2017_expt2_2.R", # trial pacing depended on subject RT, 5s time limit
            "Maniscalco_2017_expt2_3.R", # trial pacing depended on subject RT, 5s time limit
            "Maniscalco_2017_expt3.R",   # trial pacing depended on subject RT, 5s time limit
            "Maniscalco_2017_expt4.R",   # trial pacing depended on subject RT, 5s time limit
            "Massoni_2017_1.R",          # no time limit in response
            "Massoni_2017_2.R",          # no time limit in response
            "Massoni_2017_3.R",          # no time limit in response
            "Massoni_unpub_1_1.R",       # no information
            "Massoni_unpub_1_2.R",       # no information
            "Massoni_unpub_2_1.R",       # no information
            "Massoni_unpub_2_2.R",       # no information
            "Reyes_2015.R")              # instructed to respond in less than 2000 ms

data <- c()

for (j in 1:length(script)) {
  source(script[j]) 
  df$exp <- script[j] 
  data <- rbind(data, df)
}  

data <- as.data.frame(data)
data <- mutate(data, mratio_conf =  mdp_conf / dp,
                     mratio_rt =    mdp_rt / dp,
                     mratio_logit = mdp_logit / dp,
                     mdiff_conf =   mdp_conf - dp,
                     mdiff_rt =     mdp_rt - dp,
                     mdiff_logit =  mdp_logit - dp)


#'# d prime
data %>%
  select(dp, mdp_conf, mdp_rt, mdp_logit) %>%
  summary()

g1 <- GGally::ggpairs(data[, c("dp", "mdp_conf", "mdp_rt", "mdp_logit")], 
                      columnLabels = c('"d\'"', '"meta-d\'"[confidence]', '"meta-d\'"[RT]', '"meta-d\'"[confidence+RT]'),
                      labeller = label_parsed,
                      diag = list(continuous = wrap("densityDiag", size = 0.4)),
                      upper = list(continuous = wrap("cor", size = 3.5, color = "#333333")),
                      lower = list(continuous = wrap("points", alpha = 0.07))) + 
  xlim(-1, 4) + ylim(-1, 4)
g1[1, 1] <- g1[1, 1] + geom_vline(xintercept = mean(data$dp), linetype = "dashed", size = 0.3) + 
  annotate("text", x = 2.2, y = 0.6, label = "1.24", size = 3.2) + ylim(0, 0.7)
g1[1, 3] <- g1[1, 3] + xlim(-1.324082, 2.683723)
g1[2, 2] <- g1[2, 2] + geom_vline(xintercept = mean(data$mdp_conf), linetype = "dashed", size = 0.3) + 
  annotate("text", x = 2.0, y = 0.6, label = "0.94", size = 3.2) + ylim(0, 0.7)
g1[2, 3] <- g1[2, 3] + xlim(-1.324082, 2.683723)
g1[3, 3] <- g1[3, 3] + geom_vline(xintercept = mean(data$mdp_rt), linetype = "dashed", size = 0.3) + 
  annotate("text", x = 1.7, y = 0.6, label = "0.62", size = 3.2) + ylim(0, 0.7)
g1[3, 4] <- g1[3, 4] + xlim(-0.3627197, 3.8788264) + ylim(-1.324082, 2.683723)
g1[4, 4] <- g1[4, 4] + geom_vline(xintercept = mean(data$mdp_logit), linetype = "dashed", size = 0.3) + 
  annotate("text", x = 2.0, y = 0.6, label = "1.08", size = 3.2) + ylim(0, 0.7)
g1 <- g1 + theme(strip.text = element_text(size = 7))
g1
ggsave("figure_1.jpg", g1, height = 4.72, width = 4.72)
# ggplot_build(g1[2, 3])$layout$panel_scales_x[[1]]$limits


#'# stats
t.test(data$dp, data$mdp_logit, paired = T)
t.test(data$mdp_conf, data$mdp_rt, paired = T)
t.test(data$mdp_conf, data$mdp_logit, paired = T)

nrow(subset(data, data$dp < 0)) / nrow(data)
nrow(subset(data, data$mdp_conf < 0)) / nrow(data)
nrow(subset(data, data$mdp_rt < 0)) / nrow(data)
nrow(subset(data, data$mdp_logit < 0)) / nrow(data)


#'# m-ratio
data %>%
  select(dp, mratio_conf, mratio_rt, mratio_logit) %>%
  filter(dp > 0 &
         mratio_conf  > 0 & mratio_conf  < 4 &
         mratio_rt    > 0 & mratio_rt    < 4 & 
         mratio_logit > 0 & mratio_logit < 4) %>%
  summary()

data %>%
  select(dp, mratio_conf, mratio_rt, mratio_logit) %>%
  filter(dp > 0 &
         mratio_conf  > 0 & mratio_conf  < 4 &
         mratio_rt    > 0 & mratio_rt    < 4 & 
         mratio_logit > 0 & mratio_logit < 4) %>%
  GGally::ggpairs(columnLabels = c('"d\'"', '"m-ratio"[confidence]', '"m-ratio"[RT]', '"m-ratio"[confidence+RT]'),
                  labeller = label_parsed,
                  diag = list(continuous = wrap("densityDiag", size = 0.4)),
                  upper = list(continuous = wrap("cor", size = 3.5, color = "#333333")),
                  lower = list(continuous = wrap("points", alpha = 0.07))) -> g2
g2[1, 1] <- g2[1, 1] + geom_vline(xintercept = 1.38, linetype = "dashed", size = 0.3) + 
  annotate("text", x = 2.1, y = 0.5, label = "1.38", size = 3.2) + xlim(0, 4)
g2[2, 1] <- g2[2, 1] + xlim(0, 4) + ylim(0, 4)
g2[2, 2] <- g2[2, 2] + geom_vline(xintercept = 0.89, linetype = "dashed", size = 0.3) + 
  annotate("text", x = 1.8, y = 0.8, label = "0.89", size = 3.2) + xlim(0, 4)
g2[3, 1] <- g2[3, 1] + xlim(0, 4) + ylim(0, 4)
g2[3, 2] <- g2[3, 2] + xlim(0, 4) + ylim(0, 4)
g2[3, 3] <- g2[3, 3] + geom_vline(xintercept = 0.64, linetype = "dashed", size = 0.3) + 
  annotate("text", x = 1.6, y = 0.8, label = "0.64", size = 3.2) + xlim(0, 4)
g2[4, 1] <- g2[4, 1] + xlim(0, 4) + ylim(0, 4)
g2[4, 2] <- g2[4, 2] + xlim(0, 4) + ylim(0, 4)
g2[4, 3] <- g2[4, 3] + xlim(0, 4) + ylim(0, 4)
g2[4, 4] <- g2[4, 4] + geom_vline(xintercept = 0.98, linetype = "dashed", size = 0.3) + 
  annotate("text", x = 1.9, y = 0.8, label = "0.98", size = 3.2) + xlim(0, 4)
g2 <- g2 + theme(strip.text = element_text(size = 7))
g2
ggsave("figure_A3.jpg", g2, height = 4.72, width = 4.72)


#'# m-ratio (dp > 1)
data %>%
  select(dp, mratio_conf, mratio_rt, mratio_logit) %>%
  filter(dp > 1 &
           mratio_conf  > 0 & mratio_conf  < 4 &
           mratio_rt    > 0 & mratio_rt    < 4 & 
           mratio_logit > 0 & mratio_logit < 4) %>%
  summary()

data %>%
  select(dp, mratio_conf, mratio_rt, mratio_logit) %>%
  filter(dp > 1 &
           mratio_conf  > 0 & mratio_conf  < 4 &
           mratio_rt    > 0 & mratio_rt    < 4 & 
           mratio_logit > 0 & mratio_logit < 4) %>%
  GGally::ggpairs(columnLabels = c('"d\'"', '"m-ratio"[confidence]', '"m-ratio"[RT]', '"m-ratio"[confidence+RT]'),
                  labeller = label_parsed,
                  diag = list(continuous = wrap("densityDiag", size = 0.4)),
                  upper = list(continuous = wrap("cor", size = 3.5, color = "#333333")),
                  lower = list(continuous = wrap("points", alpha = 0.07))) -> g3
g3[1, 1] <- g3[1, 1] + geom_vline(xintercept = 1.74, linetype = "dashed", size = 0.3) + 
  annotate("text", x = 2.6, y = 0.5, label = "1.74", size = 3.2) + xlim(0, 4)
g3[2, 1] <- g3[2, 1] + xlim(0, 4) + ylim(0, 2.5)
g3[2, 2] <- g3[2, 2] + geom_vline(xintercept = 0.80, linetype = "dashed", size = 0.3) + 
  annotate("text", x = 1.7, y = 0.8, label = "0.80", size = 3.2) + xlim(0, 2.5)
g3[3, 1] <- g3[3, 1] + xlim(0, 4) + ylim(0, 2.5)
g3[3, 2] <- g3[3, 2] + xlim(0, 2.5) + ylim(0, 2.5)
g3[3, 3] <- g3[3, 3] + geom_vline(xintercept = 0.55, linetype = "dashed", size = 0.3) + 
  annotate("text", x = 1.4, y = 0.8, label = "0.55", size = 3.2) + xlim(0, 2.5)
g3[4, 1] <- g3[4, 1] + xlim(0, 4) + ylim(0, 2.5)
g3[4, 2] <- g3[4, 2] + xlim(0, 2.5) + ylim(0, 2.5)
g3[4, 3] <- g3[4, 3] + xlim(0, 2.5) + ylim(0, 2.5)
g3[4, 4] <- g3[4, 4] + geom_vline(xintercept = 0.85, linetype = "dashed", size = 0.3) + 
  annotate("text", x = 1.7, y = 0.8, label = "0.85", size = 3.2) + xlim(0, 2.5)
g3 <- g3 + theme(strip.text = element_text(size = 7),
                 axis.text.x = element_text(angle = 45, hjust = 1))
g3
ggsave("figure_A4.jpg", g3, height = 4.72, width = 4.72)


#'# type-1 auc
data %>%
  select(auc1_conf, auc1_rt, auc1_logit) %>%
  summary()

data %>%
  select(dp, auc1_conf, auc1_rt, auc1_logit) %>%
  GGally::ggpairs(columnLabels = c("d'", "AUC1 confidence", "AUC1 RT", "AUC1 logit"),
                  diag = list(continuous = wrap("densityDiag", size = 0.4)),
                  upper = list(continuous = wrap("cor", size = 3.5, color = "#333333")),
                  lower = list(continuous = wrap("points", alpha = 0.07))) -> g4
g4[1, 1] <- g4[1, 1] + geom_vline(xintercept = mean(data$dp),        linetype = "dashed", size = 0.3) + 
  annotate("text", x = 2, y = 0.07, label = "1.24", size = 3.2) + xlim(-1, 4)
g4[2, 1] <- g4[2, 1] + xlim(-1, 4)
g4[2, 2] <- g4[2, 2] + geom_vline(xintercept = mean(data$auc1_conf), linetype = "dashed", size = 0.3) + 
  annotate("text", x = 0.87, y = 0.5, label = "0.76", size = 3.2)
g4[3, 1] <- g4[3, 1] + xlim(-1, 4)
g4[3, 3] <- g4[3, 3] + geom_vline(xintercept = mean(data$auc1_rt), linetype = "dashed", size = 0.3) + 
  annotate("text", x = 0.86, y = 0.5, label = "0.75", size = 3.2)
g4[4, 1] <- g4[4, 1] + xlim(-1, 4) + ylim(0.4, 1)
g4[4, 2] <- g4[4, 2] + ylim(0.4, 1)
g4[4, 3] <- g4[4, 3] + ylim(0.4, 1)
g4[4, 4] <- g4[4, 4] + geom_vline(xintercept = mean(data$auc1_logit), linetype = "dashed", size = 0.3) + 
  annotate("text", x = 0.88, y = 0.5, label = "0.77", size = 3.2) + xlim(0.4, 1)
g4 <- g4 + theme(strip.text = element_text(size = 7),
                 axis.text.x = element_text(angle = 45, hjust = 1))
g4
ggsave("figure_A1.jpg", g4, height = 4.72, width = 4.72)


#'# type-2 auc
data %>%
  select(auc2_conf, auc2_rt, auc2_logit) %>%
  summary()

data %>%
  select(dp, auc2_conf, auc2_rt, auc2_logit) %>%
  GGally::ggpairs(columnLabels = c("d'", "AUC2 confidence", "AUC2 RT", "AUC2 logit"),
                  diag = list(continuous = wrap("densityDiag", size = 0.4)),
                  upper = list(continuous = wrap("cor", size = 3.5, color = "#333333")),
                  lower = list(continuous = wrap("points", alpha = 0.07))) -> g5
g5[1, 1] <- g5[1, 1] + geom_vline(xintercept = mean(data$dp), linetype = "dashed", size = 0.3) + 
  annotate("text", x = 2, y = 0.07, label = "1.24", size = 3.2) + xlim(-1, 4)
g5[2, 1] <- g5[2, 1] + xlim(-1, 4) + ylim(0.4, 1)
g5[2, 2] <- g5[2, 2] + geom_vline(xintercept = mean(data$auc2_conf), linetype = "dashed", size = 0.3) + 
  annotate("text", x = 0.73, y = 0.48, label = "0.63", size = 3.2) + xlim(0.4, 1)
g5[3, 1] <- g5[3, 1] + xlim(-1, 4) + ylim(0.4, 1)
g5[3, 2] <- g5[3, 2] + xlim(0.4, 1) + ylim(0.4, 1)
g5[3, 3] <- g5[3, 3] + geom_vline(xintercept = mean(data$auc2_rt), linetype = "dashed", size = 0.3) + 
  annotate("text", x = 0.69, y = 0.48, label = "0.59", size = 3.2) + xlim(0.4, 1)
g5[4, 1] <- g5[4, 1] + xlim(-1, 4) + ylim(0.4, 1)
g5[4, 2] <- g5[4, 2] + xlim(0.4, 1) + ylim(0.4, 1)
g5[4, 3] <- g5[4, 3] + xlim(0.4, 1) + ylim(0.4, 1)
g5[4, 4] <- g5[4, 4] + geom_vline(xintercept = mean(data$auc2_logit), linetype = "dashed", size = 0.3) + 
  annotate("text", x = 0.76, y = 0.48, label = "0.66", size = 3.2) + xlim(0.4, 1)
g5 <- g5 + theme(strip.text = element_text(size = 7),
                 axis.text.x = element_text(angle = 45, hjust = 1))
g5
ggsave("figure_A2.jpg", g5, height = 4.72, width = 4.72)


#'# meta-d' boost
data %>%
  select(dp, mratio_conf, mratio_rt, mratio_logit) %>%
  filter(dp > 0 &
           mratio_conf  > 0 & mratio_conf  < 4 &
           mratio_rt    > 0 & mratio_rt    < 4 & 
           mratio_logit > 0 & mratio_logit < 4) %>% 
  mutate(ratio =mratio_logit / mratio_conf) %>%
  filter(ratio < 4) %>% # 406 cases left
  summary()

data %>%
  select(dp, mratio_conf, mratio_rt, mratio_logit) %>%
  filter(dp > 0 &
           mratio_conf  > 0 & mratio_conf  < 4 &
           mratio_rt    > 0 & mratio_rt    < 4 & 
           mratio_logit > 0 & mratio_logit < 4) %>% # 412 cases left
  mutate(ratio =mratio_logit / mratio_conf) %>%
  filter(ratio < 4) %>%
  ggplot() + geom_density(aes(x = ratio)) +
  geom_vline(xintercept = 1.1757, linetype = "dashed")

# Across-subject median split on mean RT
p1 <- ggplot(data, aes(x = dp, y = mdp_conf, color = group)) +
  geom_point(alpha = 0.4) + xlim(0, 4) + ylim(0, 4) +
  geom_smooth(method = "loess", span = 0.75, se = FALSE)


data %>%
  group_by(group) %>%
  select(dp, mdp_conf, mratio_conf) %>%
  filter(dp > 0 &
           mratio_conf  > 0 & mratio_conf  < 4) %>%
  summarise(mean_dp = mean(dp), mean_mdp_conf = mean(mdp_conf), mean_mratio_conf = mean(mratio_conf))


p2 <- data %>%
  filter(dp > 0 & mratio_conf > 0 & mratio_conf < 4) %>%
  ggplot() + 
  geom_violin(aes(x = group, y = dp, color = group)) +
  stat_summary(aes(x = group, y = dp), 
               fun = mean, geom = "point", 
               shape = 2, size = 3, color = "black") +
  stat_summary(aes(x = group, y = dp, label = round(after_stat(y), 2)), 
               fun = mean, geom = "text", 
               vjust = -1, size = 4, color = "black")


p3 <- data %>%
  filter(dp > 0 & mratio_conf > 0 & mratio_conf < 4) %>%
  ggplot() + 
  geom_violin(aes(x = group, y = mdp_conf, color = group)) +
  stat_summary(aes(x = group, y = mdp_conf), 
               fun = mean, geom = "point", 
               shape = 2, size = 3, color = "black") +
  stat_summary(aes(x = group, y = mdp_conf, label = round(after_stat(y), 2)), 
               fun = mean, geom = "text", 
               vjust = -1, size = 4, color = "black")


p4 <- data %>%
  filter(dp > 0 & mratio_conf > 0 & mratio_conf < 4) %>%
  ggplot() + 
  geom_violin(aes(x = group, y = mratio_conf, color = group)) +
  stat_summary(aes(x = group, y = mratio_conf), 
               fun = mean, geom = "point", 
               shape = 2, size = 3, color = "black") +
  stat_summary(aes(x = group, y = mratio_conf, label = round(after_stat(y), 2)), 
               fun = mean, geom = "text", 
               vjust = -1, size = 4, color = "black") 

p <- plot_grid(p1, p2, p3, p4, nrow = 2, labels = c("a", "b", "c", "d"))

p <- ggdraw(p) + ggtitle("Across-subject median split on mean RT") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        plot.title.position = "plot")
p