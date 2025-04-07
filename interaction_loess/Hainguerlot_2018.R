



library(tidyverse)
library(data.table)
library(doParallel)

cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

dat <- fread("data_Hainguerlot_2018.csv", header = T)

df1 <- foreach(i = 1:length(unique(dat$Subj_idx)), .combine = "rbind", .packages = c("dplyr", "metaSDT")) %dopar% {
  
  d <- subset(dat, dat$Subj_idx == unique(dat$Subj_idx)[i])
  d <- na.omit(d)
  d <- mutate(d, Correct = ifelse(Stimulus == Response, 1, 0))
  d <- mutate(d, Confidence = ifelse(Confidence < 0.7, 1,
                                     ifelse(Confidence < 0.9, 2, 3)))
  
  # RT bin
  q <- unique(quantile(d$RT_dec, probs = seq(0, 1, by = 0.05), na.rm = TRUE))
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

df2 <- foreach(i = 1:length(unique(dat$Subj_idx)), .combine = "rbind", .packages = c("dplyr", "metaSDT")) %dopar% {
  
  d <- subset(dat, dat$Subj_idx == unique(dat$Subj_idx)[i])
  d <- na.omit(d)
  d <- mutate(d, Correct = ifelse(Stimulus == Response, 1, 0))
  d <- mutate(d, Confidence = ifelse(Confidence < 0.7, 1,
                                     ifelse(Confidence < 0.9, 2, 3)))
  
  # RT bin
  q <- unique(quantile(d$RT_dec, probs = seq(0, 1, by = 0.05), na.rm = TRUE))
  d$RT_bin <- cut(d$RT_dec, breaks = q, labels = FALSE, include.lowest = TRUE)
  d$RT_bin <- as.factor(d$RT_bin)
  print(d)
  
}



