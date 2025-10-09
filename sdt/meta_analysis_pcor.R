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
  pcor$exp <- dataset[j] 
  data <- rbind(data, pcor)
}  

# Fisher's Z Transformation
data_z <- data %>%
  mutate(across(
    .cols = 1,
    .fns = ~ atanh(.x),
    .names = "{.col}_z"
  )) %>%
  mutate(across(
    .cols = 1,
    .fns = ~ 1 / (n - 3), # variance for Z
    .names = "{.col}_var"
  ))

# random effect model
res <- rma(yi = estimate_z, vi = estimate_var, data = data_z, method = "REML")
summary(res)

z <- 0.3007
ci_lb_z <- 0.2119
ci_ub_z <- 0.3894

# z -> r
r <- (exp(2*z) - 1) / (exp(2*z) + 1)
ci_lb <- (exp(2*ci_lb_z) - 1) / (exp(2*ci_lb_z) + 1)
ci_ub <- (exp(2*ci_ub_z) - 1) / (exp(2*ci_ub_z) + 1)

r; ci_lb; ci_ub