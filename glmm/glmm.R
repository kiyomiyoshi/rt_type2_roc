#+ message = F, warning = F
library(tidyverse)
library(gghalves)
library(ggeffects)
library(GGally)
library(car)
library(doParallel)
library(sjPlot)
library(ggstance) 

source("theme_publication.R")
theme_set(theme_publication()) 

cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

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

exp <- sub("\\.R$", "", script)


#'# model fits
results <- foreach(j = 1:length(script), .combine = "rbind", .multicombine = T) %dopar% {
  source(script[j]) 
  list(Anova(m3), summary(m3), Anova(m2), summary(m2), roc_1, roc_2, roc_1_sdt, roc_2_sdt, sdt, g1, g2, g3, script[j])
}  


#'# GLMM significance
glmm_list <- results[, 1]
Pr <- foreach(k = 1:length(script), .combine = "rbind") %dopar% {
  glmm_list[[k]]$`Pr(>Chisq)`
}

apply(Pr, 2, function(x) sum(x < 0.05))
matrix(as.numeric(Pr < 0.05), ncol = ncol(Pr))
Pr <- as.data.frame(Pr < 0.05)
rownames(Pr) <- gsub("^result\\.", "", exp)
colnames(Pr) <- c("Confidence", "RT", "Confidence * RT")
Pr


#'# Visualization
summary_list <- results[, 4]
rt <- foreach(k = 1:length(script), .combine = "rbind") %dopar% {
  summary_list[[k]]$coefficients[3, 1:2]
}
rt <- as.data.frame(rt)
rt$Upper <- rt[, 1] + 1.96 * rt[, 2]
rt$Lower <- rt[, 1] - 1.96 * rt[, 2]
rt$Variable <- "RT"

summary_list <- results[, 4]
confidence <- foreach(k = 1:length(script), .combine = "rbind") %dopar% {
  summary_list[[k]]$coefficients[2, 1:2]
}
confidence <- as.data.frame(confidence)
confidence$Upper <- confidence[, 1] + 1.96 * confidence[, 2]
confidence$Lower <- confidence[, 1] - 1.96 * confidence[, 2]
confidence$Variable <- "confidence"

summary_list <- results[, 2]
interaction <- foreach(k = 1:length(script), .combine = "rbind") %dopar% {
  summary_list[[k]]$coefficients[4, 1:2]
}
interaction <- as.data.frame(interaction)
interaction$Upper <- interaction[, 1] + 1.96 * interaction[, 2]
interaction$Lower <- interaction[, 1] - 1.96 * interaction[, 2]
interaction$Variable <- "interaction"

df <- rbind(rt, confidence, interaction)
df$Variable <- factor(df$Variable, levels = c("interaction", "confidence", "RT"))
Dataset <- c("Hainguerlot_2018", "Hainguerlot_unpub", "Maniscalco_2017_expt1", "Maniscalco_2017_expt2_1",
              "Maniscalco_2017_expt2_2", "Maniscalco_2017_expt2_3", "Maniscalco_2017_expt3", "Maniscalco_2017_expt4",
              "Massoni_2017_1", "Massoni_2017_2", "Massoni_2017_3", "Massoni_unpub_1_1", 
              "Massoni_unpub_1_2", "Massoni_unpub_2_1", "Massoni_unpub_2_2", "Reyes_2015")
df$Dataset <- rep(Dataset, 3)
df$Dataset <- factor(df$Dataset, levels = rev(Dataset))

df$Variable <- factor(df$Variable, levels = c("interaction", "RT", "confidence"))
ggplot(df, aes(x = Estimate, y = Dataset, 
               color = Variable, shape = Variable)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.35)) +
  geom_point(position = ggstance::position_dodgev(height = 1.0), size = 2.5) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper),
                position = ggstance::position_dodgev(height = 1.0), width = 0.6) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  theme_light(base_size = 12) +
  labs(x = "Estimate", y = NULL, color = "Variable", shape = "Variable") +
  theme(
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "top",
    legend.margin = margin(t = 0, b = 0)
  ) +
  scale_color_manual(
    limits = c("confidence", "RT", "interaction"),
    values = c("confidence" = "#E69F00", "RT" = "#009E73", "interaction" = "#0072B2")
  ) +
  scale_shape_manual(
    limits = c("confidence", "RT", "interaction"),
    values = c("confidence" = 16, "RT" = 17, "interaction" = 15)
  ) -> p1
p1
ggsave("figure_2.jpg", p1, width = 14, height = 13, units = "cm", dpi = 300)

#'# type-1 roc from glmm logit
roc_1_list <- results[, 5]
roc_1 <- foreach(k = 1:length(script), .combine = "rbind") %dopar% {
  roc_1_list[[k]]
}

level_mapping <- c("1" = "Model 1", "2" = "Model 2", "3" = "Model 3")
roc_1 %>%
  group_by(variable, order) %>%
  summarize(mean_hr_1 = mean(hr_1), mean_far_1 = mean(far_1)) %>%
  mutate(model = factor(variable, levels = names(level_mapping), labels = level_mapping)) %>%
  ggplot() + geom_point(aes(x = mean_far_1, y = mean_hr_1, color = model), size = 0.8, alpha = 0.8) +
  theme(legend.position = c(0.75, 0.3), legend.direction = "vertical") + 
  guides(color = guide_legend(title = NULL)) + 
  xlab("Type-1 hit rate") + ylab("Type-1 FA rate") -> p2
p2


#'# type-2 roc from glmm logit
roc_2_list <- results[, 6]
roc_2 <- foreach(k = 1:length(script), .combine = "rbind") %dopar% {
  roc_2_list[[k]]
}

level_mapping <- c("1" = "Model 1", "2" = "Model 2", "3" = "Model 3")
roc_2 %>%
  group_by(variable, order) %>%
  summarize(mean_hr_2 = mean(hr_2), mean_far_2 = mean(far_2)) %>%
  mutate(model = factor(variable, levels = names(level_mapping), labels = level_mapping)) %>%
  ggplot() + geom_point(aes(x = mean_far_2, y = mean_hr_2, color = model), size = 0.8, alpha = 0.8) +
  theme(legend.position = c(0.75, 0.3), legend.direction = "vertical") + 
  guides(color = guide_legend(title = NULL)) + 
  xlab("Type-2 hit rate") + ylab("Type-2 FA rate") -> p3
p3


#'# sdt stats from glmm logit
sdt_list <- results[, 9]
sdt <- foreach(k = 1:length(script), .combine = "rbind") %dopar% {
  sdt_list[[k]]
}

# summary table
sdt %>%
  group_by(variable) %>%
  summarise(across(c(da, meta_da, M_diff, M_ratio), mean))

# m-ratio
sdt %>%
  mutate(index = factor(variable, levels = c("conf_logit", "RT_logit", "interaction_logit"))) %>%
  ggplot(aes(x = variable, y = M_ratio, color = variable)) + geom_half_violin() + geom_point() +
  stat_summary(fun = "mean", size = 0.7, geom = "crossbar", width = 0.2, color = "black") +
  theme(axis.title.x = element_blank(), legend.position = "none") + 
  scale_x_discrete(labels = c("Model 1", "Model 2", "Model 3")) + ylab("M-ratio") -> p4
p4

# d prime
sdt_list <- results[, 9]
dp <- foreach(k = 1:length(script), .combine = "rbind") %dopar% {
  df <- as.data.frame(sdt_list[[k]])
  cbind(c(df[1, 1], df[, 2]), c("dp", "mdp_conf", "mdp_add", "mdp_int"))
}

colors = c("dp" = "black", "mdp_conf" = "#F8766D", "mdp_add" = "#00BA38", "mdp_int" = "#619CFF")
as.data.frame(dp) %>% 
  mutate(V1 = as.numeric(V1), V2 = factor(V2, levels = c("dp", "mdp_conf", "mdp_add", "mdp_int"))) %>%
  ggplot(aes(x = V2, y = V1, color = V2)) + geom_half_violin() + geom_point() +
  scale_color_manual(values = colors) +
  stat_summary(fun = "mean", size = 0.7, geom = "crossbar", width = 0.2, color = "black") +
  theme(axis.title.x = element_blank(), legend.position = "none")  + 
  scale_x_discrete(labels = c("d'", "meta-d'\nModel 1", "meta-d'\nModel 2", "meta-d'\nModel 3")) + ylab("value") -> p5
p5


#'# t tests on glmm m-ratio
pairwise.t.test(sdt[, 4], sdt[, 5], p.adjust.method = "holm", paired = T)