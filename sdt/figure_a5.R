#+ message = F, warning = F
library(tidyverse)
library(gghalves)
library(GGally)
library(car)
library(doParallel)
library(sjPlot)
library(psych)
library(scales)

source("theme_publication.R")
theme_set(theme_publication()) 

#'# summary stats
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
data <- mutate(data, 
               mratio_logit = mdp_logit / dp,
               mdpr_mdpc    = mdp_rt / mdp_conf,
               mdpl_mdpc    = mdp_logit / mdp_conf)

data %>%
  filter(dp > 0 &
  mratio_logit > 0 & mratio_logit  < 4 &
  mdpr_mdpc    > 0 & mdpr_mdpc    < 4 & 
  mdpl_mdpc    > 0 & mdpl_mdpc < 4) -> data

stats <- describeBy(data[, -c(11)], group = data$exp, digits = 3)
print(stats)


#'# visualization
dataset <- c("Hainguerlot_2018", "Hainguerlot_unpub", "Maniscalco_2017_expt1", "Maniscalco_2017_expt2_1",
             "Maniscalco_2017_expt2_2", "Maniscalco_2017_expt2_3", "Maniscalco_2017_expt3", "Maniscalco_2017_expt4",
             "Massoni_2017_1", "Massoni_2017_2", "Massoni_2017_3", "Massoni_unpub_1_1", 
             "Massoni_unpub_1_2", "Massoni_unpub_2_1", "Massoni_unpub_2_2", "Reyes_2015")

mrl <- foreach(i = 1:length(script), .combine = "rbind") %dopar% {
  c(stats[[i]][14, ]$mean, stats[[i]][14, ]$se, "mrl", dataset[i])
}
mdpr_c <- foreach(i = 1:length(script), .combine = "rbind") %dopar% {
  c(stats[[i]][15, ]$mean, stats[[i]][15, ]$se, "mdpr_c", dataset[i])
}
mdpl_c <- foreach(i = 1:length(script), .combine = "rbind") %dopar% {
  c(stats[[i]][16, ]$mean, stats[[i]][16, ]$se, "mdpl_c", dataset[i])
}

df1 <- as.data.frame(rbind(mrl, mdpr_c, mdpl_c))
colnames(df1) <- c("value", "se", "index", "dataset")
df1$value <- as.numeric(df1$value)
df1$se <- as.numeric(df1$se)

df1 %>%
  group_by(index) %>%
  summarize(mean = mean(value), se <- sd(value) / sqrt(16)) -> dsav
dsav$dataset <- "Across datasets"
colnames(dsav) <- c("index", "value", "se", "dataset")
df1 <- rbind(df1, dsav)
df1$dataset <- as.factor(df1$dataset)
df1$dataset <- factor(df1$dataset, levels = c("Across datasets", rev(dataset)))

df1 %>%
  mutate(index = fct_relevel(index,  "mdpl_c", "mdpr_c", "mrl")) %>%
  ggplot(aes(x = value, y = dataset, color = index)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = gray(0.35)) +
  geom_point(aes(shape = ifelse(dataset == "Across datasets", "triangle", "circle")), size = 3,
             position = position_dodge(width = 0.7), show.legend = F) +
  geom_errorbar(aes(xmin = value - se, xmax = value + se), 
                position = position_dodge(width = 0.7), width = 0.4) +
  theme_light(base_size = 12) +  # Use a clean theme
  labs(x = "value", y = NULL, color = "index") +
  theme(axis.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        legend.position = "top",
        legend.justification = c(1, 1.2)) +
  scale_color_manual(values = c(hue_pal()(3)[3], "#FFC425", "#440154FF"), 
                     limits = c("mrl", "mdpr_c", "mdpl_c"),
                     labels = c(expression(paste(m, "-", ratio[confidence+RT])),
                                expression(paste(meta, "-", d, "'"[RT]/meta, "-", d, "'"[confidence])), 
                                expression(paste(meta, "-", d, "'"[confidence+RT]/meta, "-", d, "'"[confidence])))) +
  guides(color = guide_legend(nrow = 2)) -> p1
p1
ggsave("figure_a5.jpg", p1, dpi = 300)


#'# t tests
data %>%
  group_by(exp) %>%
  summarise(dp_mdpl_t =   t.test(dp, mdp_logit, paired = TRUE)$statistic,
            dp_mdpl_p =   t.test(dp, mdp_logit, paired = TRUE)$p.value,
            mdpc_mdpr_t = t.test(mdp_conf, mdp_rt, paired = TRUE)$statistic,
            mdpc_mdpr_p = t.test(mdp_conf, mdp_rt, paired = TRUE)$p.value)