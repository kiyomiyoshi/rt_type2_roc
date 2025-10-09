#+ message = F, warning = F
library(tidyverse)
library(metafor)
library(gridExtra)
library(sjPlot)
library(ggtext)

source("theme_publication.R")
theme_set(theme_publication()) 

script <- c("Hainguerlot_2018.R",       
            "Hainguerlot_unpub.R",      
            "Maniscalco_2017_expt1.R",  
            "Maniscalco_2017_expt2_1.R",
            "Maniscalco_2017_expt2_2.R",
            "Maniscalco_2017_expt2_3.R",
            "Maniscalco_2017_expt3.R",
            "Maniscalco_2017_expt4.R",
            "Massoni_2017_1.R",
            "Massoni_2017_2.R",
            "Massoni_2017_3.R",
            "Massoni_unpub_1_1.R",
            "Massoni_unpub_1_2.R",
            "Massoni_unpub_2_1.R",
            "Massoni_unpub_2_2.R",
            "Reyes_2015.R") 

dataset <- c("Hainguerlot_2018",       
             "Hainguerlot_unpub",      
             "Maniscalco_2017_expt1",  
             "Maniscalco_2017_expt2_1",
             "Maniscalco_2017_expt2_2",
             "Maniscalco_2017_expt2_3",
             "Maniscalco_2017_expt3",
             "Maniscalco_2017_expt4",
             "Massoni_2017_1",
             "Massoni_2017_2",
             "Massoni_2017_3",
             "Massoni_unpub_1_1",
             "Massoni_unpub_1_2",
             "Massoni_unpub_2_1",
             "Massoni_unpub_2_2",
             "Reyes_2015") 

data <- c()

for (j in 1:length(script)) {
  source(script[j]) 
  cor_df_mratio$exp <- dataset[j] 
  data <- rbind(data, cor_df_mratio)
}  

# Fisher's Z Transformation
data_z <- data %>%
  mutate(across(
    .cols = 1:6,
    .fns = ~ atanh(.x),
    .names = "{.col}_z"
  )) %>%
  mutate(across(
    .cols = 1:6,
    .fns = ~ 1 / (n - 3), # variance for Z
    .names = "{.col}_var"
  ))

# random effect model
results <- list()
pairs <- colnames(data)[1:6]

for (pair in pairs) {
  z_col <- paste0(pair, "_z")
  var_col <- paste0(pair, "_var")
  
  res <- rma(yi = data_z[[z_col]], vi = data_z[[var_col]], method = "REML")
  results[[pair]] <- res
}

# visualization by base functions
plots <- list()

par(mfrow = c(2,3)) # 2 x 3
for (pair in pairs) {
  res <- results[[pair]]
  forest(res, slab = data$exp, atransf = tanh, xlab = "r", main = pair, cex = 1)
}
par(mfrow = c(1,1))

# visualization by ggplot
pairs <- colnames(data)[1:6]
results <- list()

# pairwise meta analysis
for (p in pairs) {
  r_vals <- data[[p]]
  n_vals <- data$n
  z <- atanh(r_vals)
  var <- 1/(n_vals - 3)
  res <- rma(yi = z, vi = var, method = "REML")
  results[[p]] <- res
}

# data frames
df_forest <- map_dfr(pairs, function(p){
  df_pair <- data %>%
    select(exp, n, !!p) %>%
    rename(r = !!p) %>%
    mutate(
      pair = p,
      z = atanh(r),
      var = 1/(n - 3),
      ci_lb_r_study = tanh(z - 1.96 * sqrt(var)),
      ci_ub_r_study = tanh(z + 1.96 * sqrt(var))
    )
  df_pair
})

df_meta_rows <- map_dfr(pairs, function(p){
  res <- results[[p]]
  tibble(
    exp = "Meta-analysis",
    n = NA,
    r = NA,
    pair = p,
    z = NA,
    var = NA,
    ci_lb_r_study = NA,
    ci_ub_r_study = NA,
    est_r_meta = tanh(res$b),
    ci_lb_r_meta = tanh(res$ci.lb),
    ci_ub_r_meta = tanh(res$ci.ub)
  )
})

df_forest <- bind_rows(df_forest, df_meta_rows)

df_forest <- df_forest %>%
  group_by(pair) %>%
  mutate(
    exp = factor(exp, levels = rev(c(setdiff(unique(exp), "Meta-analysis"), "Meta-analysis")))
  ) %>%
  arrange(pair, exp) %>%
  ungroup()

df_forest_relabeled <- df_forest %>%
  mutate(pair = factor(pair,
                       levels = c("dp-mratio_conf", "dp-mratio_rt", "dp-mratio_logit",
                                  "mratio_conf-mratio_rt", "mratio_conf-mratio_logit", "mratio_rt-mratio_logit"),
                       labels = c(
                         "d*\"'\"*\" vs. \"*m*\"-\"*ratio*\"\"[confidence]",
                         "d*\"'\"*\" vs. \"*m*\"-\"*ratio*\"\"[RT]",
                         "d*\"'\"*\" vs. \"*m*\"-\"*ratio*\"\"[confidence+RT]",
                         "m*\"-\"*ratio*\"\"[confidence]*\" vs. \"*m*\"-\"*ratio*\"\"[RT]",
                         "m*\"-\"*ratio*\"\"[confidence]*\" vs. \"*m*\"-\"*ratio*\"\"[confidence+RT]",
                         "m*\"-\"*ratio*\"\"[RT]*\" vs. \"*m*\"-\"*ratio*\"\"[confidence+RT]"
                       )))


ggplot(df_forest_relabeled, aes(y = exp, x = r)) +
  geom_errorbarh(aes(xmin = ci_lb_r_study, xmax = ci_ub_r_study), alpha = 0.6) +
  geom_point(color = "darkgray") +
  geom_point(data = df_forest_relabeled %>% filter(exp == "Meta-analysis"),
             aes(x = est_r_meta), shape = 18, size = 3) +
  geom_errorbarh(data = df_forest_relabeled %>% filter(exp == "Meta-analysis"),
                 aes(y = exp, xmin = ci_lb_r_meta, xmax = ci_ub_r_meta),
                 height = 0.3) +
  facet_wrap(~pair, labeller = label_parsed) +
  labs(x = "Correlation coefficient", y = NULL) +
  theme_bw() +
  coord_cartesian(xlim = c(-1.03, 1.03)) +
  theme(axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        strip.text = element_text(size = 8)) -> g
g
ggsave("figure_a5.png", g, height = 4.72 * 1.1, width = 4.72 * 1.55)