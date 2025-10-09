#+ message = F, warning = F
library(tidyverse)
library(gghalves)
library(GGally)
library(car)
library(doParallel)
library(sjPlot)
library(ggdist)

source("theme_publication.R")
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

data_sh <- c()

for (j in 1:length(script)) {
  source(script[j]) 
  df$exp <- dataset[j] 
  data_sh <- rbind(data_sh, df)
}  

data_sh <- data_sh %>%
  group_by(exp, sub) %>% 
  filter(all(c("odd", "even") %in% rn)) %>%
  ungroup()

data_filtered <- data %>%
  inner_join(dplyr::select(data_sh, exp, sub), 
             by = c("Subj_idx" = "sub", "exp" = "exp"))
nrow(data_filtered) # 970
data_filtered <- unique(data_filtered)
nrow(data_filtered) # 485

data_sh_filtered <- data_sh %>%
  dplyr::inner_join(
    dplyr::select(data_filtered, exp, Subj_idx),
    by = c("sub" = "Subj_idx", "exp" = "exp"))
nrow(data_sh_filtered) # 970

# write.csv(data_filtered, "data_filtered.csv", row.names = FALSE)
# write.csv(data_sh_filtered, "data_sh_filtered.csv", row.names = FALSE)

