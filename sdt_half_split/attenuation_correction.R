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

data <- data %>%
  group_by(exp, sub) %>% 
  filter(all(c("odd", "even") %in% rn)) %>%
  ungroup()

wide_data <- data %>%
  select(sub, exp, rn, mdp_conf, mdp_rt) %>%
  pivot_wider(names_from = rn, values_from = c(mdp_conf, mdp_rt), names_sep = "_")


#'# d prime
summary(wide_data)

g1 <- GGally::ggpairs(wide_data,
                      columns = c("mdp_conf_odd", "mdp_conf_even", "mdp_rt_odd", "mdp_rt_even"),
                      columnLabels = c('"meta-d\'"[confidence]~(odd)', '"meta-d\'"[confidence]~(even)', '"meta-d\'"[RT]~(odd)', '"meta-d\'"[RT]~(even)'),
                      labeller = label_parsed,
                      diag = list(continuous = wrap("densityDiag", size = 0.4)),
                      upper = list(continuous = wrap("cor", size = 3.5, color = "#333333")),
                      lower = list(continuous = wrap("points", alpha = 0.07)))

g1[1, 1] <- g1[1, 1] + geom_vline(xintercept = mean(wide_data$mdp_conf_odd), linetype = "dashed", size = 0.3) + 
  annotate("text", x = 2.2, y = 0.6, label = "0.89", size = 3.2) + xlim(-1, 4) + ylim(0, 0.7)
g1[2, 1] <- g1[2, 1] + xlim(-1, 4) + ylim(-1, 4)
g1[3, 1] <- g1[3, 1] + xlim(-1, 4) + ylim(-1, 4)
g1[4, 1] <- g1[4, 1] + xlim(-1, 4) + ylim(-1, 4)
g1[2, 2] <- g1[2, 2] + geom_vline(xintercept = mean(wide_data$mdp_conf_even), linetype = "dashed", size = 0.3) + 
  annotate("text", x = 2.2, y = 0.6, label = "0.89", size = 3.2) + xlim(-1, 4) + ylim(0, 0.7)
g1[3, 2] <- g1[3, 2] + xlim(-1, 4) + ylim(-1, 4)
g1[4, 2] <- g1[4, 2] + xlim(-1, 4) + ylim(-1, 4)
g1[3, 3] <- g1[3, 3] + geom_vline(xintercept = mean(wide_data$mdp_rt_odd), linetype = "dashed", size = 0.3) + 
  annotate("text", x = 1.9, y = 0.6, label = "0.59", size = 3.2) + xlim(-1, 4) + ylim(0, 0.7)
g1[4, 3] <- g1[4, 3] + xlim(-1, 4) + ylim(-1, 4)
g1[4, 4] <- g1[4, 4] + geom_vline(xintercept = mean(wide_data$mdp_rt_even), linetype = "dashed", size = 0.3) + 
  annotate("text", x = 1.9, y = 0.6, label = "0.57", size = 3.2) + xlim(-1, 4) + ylim(0, 0.7)
g1 <- g1 + theme(strip.text = element_text(size = 7))
g1
ggsave("figure_A.jpg", g1, height = 5.2, width = 5.2)

# attenuation correction
reliability_rt   <- (2 * 0.431) / (1 + 0.431)
reliability_conf <- (2 * 0.588) / (1 + 0.588)

0.518 / sqrt(reliability_rt * reliability_conf)