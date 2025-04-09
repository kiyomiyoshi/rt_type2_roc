library(metaSDT)
library(tidyverse)
library(GGally)
library(data.table)
library(car)
library(doParallel)

cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

dat <- fread("data_Massoni_unpub.csv", header = T)
dat <- subset(dat, dat$Study == "2")
dat <- subset(dat, dat$Difficulty == "2")
dat <- subset(dat, dat$Confidence >= 0.5)

df1 <- foreach(i = 1:length(unique(dat$Subj_idx)), .combine = "rbind", .packages = c("dplyr", "metaSDT")) %dopar% {
  
  d <- subset(dat, dat$Subj_idx == unique(dat$Subj_idx)[i])
  d <- filter(d, row_number() %% 2 == 1)
  d <- na.omit(d)
  d <- mutate(d, Correct = ifelse(Stimulus == Response, 1, 0))
  
  # confidence
  nR_S1 <- c(sum(d$Stimulus == "1" & d$Response == "1" & d$Confidence == 1.0, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$Confidence == 0.9, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$Confidence == 0.8, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$Confidence == 0.7, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$Confidence == 0.6, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$Confidence == 0.5, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$Confidence == 0.5, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$Confidence == 0.6, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$Confidence == 0.7, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$Confidence == 0.8, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$Confidence == 0.9, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$Confidence == 1.0, na.rm = TRUE))
  
  nR_S2 <- c(sum(d$Stimulus == "2" & d$Response == "1" & d$Confidence == 1.0, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$Confidence == 0.9, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$Confidence == 0.8, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$Confidence == 0.7, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$Confidence == 0.6, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$Confidence == 0.5, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$Confidence == 0.5, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$Confidence == 0.6, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$Confidence == 0.7, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$Confidence == 0.8, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$Confidence == 0.9, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$Confidence == 1.0, na.rm = TRUE))
  
  fit_c <- data.frame(NA, NA, NA)
  tryCatch({
    fit_c <- fit_meta_d_MLE(nR_S1, nR_S2)
  }, error = function(e) {
    cat("Error: Fitting the model failed with message:", e$message, "\n")
  })
  
  n_hit <- cumsum(nR_S1)
  n_fa  <- cumsum(nR_S2)
  p_hit <- n_hit/sum(nR_S1)
  p_fa  <- n_fa/sum(nR_S2)
  auc1_c <- p_hit[1] * p_fa[1] / 2
  for (n in 1:(length(nR_S1) - 1)) {
    auc1_c <- auc1_c + (p_hit[n] + p_hit[n + 1]) * (p_fa[n + 1] - p_fa[n]) / 2
  }
  
  cl <- length(nR_S1) / 2
  hit2 <- nR_S1[1:cl] + rev(nR_S2[(cl + 1):(cl * 2)])
  fa2  <- rev(nR_S1[(cl + 1):(cl * 2)]) + nR_S2[1:cl]
  n_hit2 <- cumsum(hit2)
  n_fa2  <- cumsum(fa2)
  p_hit2 <- n_hit2/sum(hit2)
  p_fa2  <- n_fa2/sum(fa2)
  auc2_c <- p_hit2[1] * p_fa2[1] / 2
  for (n in 1:(cl - 1)) {
    auc2_c <- auc2_c + (p_hit2[n] + p_hit2[n + 1]) * (p_fa2[n + 1] - p_fa2[n]) / 2
  }
  
  # RT
  q <- quantile(d$RT_dec, probs = seq(0, 1, by = 1/6), na.rm = TRUE)
  d$RT_bin <- cut(d$RT_dec, breaks = q, labels = FALSE, include.lowest = TRUE)
  
  nR_S1 <- c(sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 5, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 6, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$RT_bin == 6, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$RT_bin == 5, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$RT_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$RT_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$RT_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$RT_bin == 1, na.rm = TRUE))
  
  nR_S2 <- c(sum(d$Stimulus == "2" & d$Response == "1" & d$RT_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$RT_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$RT_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$RT_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$RT_bin == 5, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$RT_bin == 6, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$RT_bin == 6, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$RT_bin == 5, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$RT_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$RT_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$RT_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$RT_bin == 1, na.rm = TRUE))
  
  fit_r <- data.frame(NA, NA, NA)
  tryCatch({
    fit_r <- fit_meta_d_MLE(nR_S1, nR_S2)
  }, error = function(e) {
    cat("Error: Fitting the model failed with message:", e$message, "\n")
  })
  
  n_hit <- cumsum(nR_S1)
  n_fa  <- cumsum(nR_S2)
  p_hit <- n_hit/sum(nR_S1)
  p_fa  <- n_fa/sum(nR_S2)
  auc1_r <- p_hit[1] * p_fa[1] / 2
  for (n in 1:(length(nR_S1) - 1)) {
    auc1_r <- auc1_r + (p_hit[n] + p_hit[n + 1]) * (p_fa[n + 1] - p_fa[n]) / 2
  }
  
  cl <- length(nR_S1) / 2
  hit2 <- nR_S1[1:cl] + rev(nR_S2[(cl + 1):(cl * 2)])
  fa2  <- rev(nR_S1[(cl + 1):(cl * 2)]) + nR_S2[1:cl]
  n_hit2 <- cumsum(hit2)
  n_fa2  <- cumsum(fa2)
  p_hit2 <- n_hit2/sum(hit2)
  p_fa2  <- n_fa2/sum(fa2)
  auc2_r <- p_hit2[1] * p_fa2[1] / 2
  for (n in 1:(cl - 1)) {
    auc2_r <- auc2_r + (p_hit2[n] + p_hit2[n + 1]) * (p_fa2[n + 1] - p_fa2[n]) / 2
  }
  
  # logit
  m <- glm(Correct ~ Confidence * RT_dec, family = binomial, data = d)
  d$logit <- predict(m)
  q <- quantile(d$logit, probs = seq(0, 1, by = 1/6), na.rm = TRUE)
  d$logit_bin <- cut(d$logit, breaks = q, labels = FALSE, include.lowest = TRUE)
  
  nR_S1 <- c(sum(d$Stimulus == "1" & d$Response == "1" & d$logit_bin == 6, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$logit_bin == 5, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$logit_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$logit_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$logit_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$logit_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$logit_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$logit_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$logit_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$logit_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$logit_bin == 5, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$logit_bin == 6, na.rm = TRUE))
  
  nR_S2 <- c(sum(d$Stimulus == "2" & d$Response == "1" & d$logit_bin == 6, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$logit_bin == 5, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$logit_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$logit_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$logit_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$logit_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$logit_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$logit_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$logit_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$logit_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$logit_bin == 5, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$logit_bin == 6, na.rm = TRUE))
  
  fit_l <- data.frame(NA, NA, NA)
  tryCatch({
    fit_l <- fit_meta_d_MLE(nR_S1, nR_S2)
  }, error = function(e) {
    cat("Error: Fitting the model failed with message:", e$message, "\n")
  })
  
  n_hit <- cumsum(nR_S1)
  n_fa  <- cumsum(nR_S2)
  p_hit <- n_hit/sum(nR_S1)
  p_fa  <- n_fa/sum(nR_S2)
  auc1_l <- p_hit[1] * p_fa[1] / 2
  for (n in 1:(length(nR_S1) - 1)) {
    auc1_l <- auc1_l + (p_hit[n] + p_hit[n + 1]) * (p_fa[n + 1] - p_fa[n]) / 2
  }
  
  cl <- length(nR_S1) / 2
  hit2 <- nR_S1[1:cl] + rev(nR_S2[(cl + 1):(cl * 2)])
  fa2  <- rev(nR_S1[(cl + 1):(cl * 2)]) + nR_S2[1:cl]
  n_hit2 <- cumsum(hit2)
  n_fa2  <- cumsum(fa2)
  p_hit2 <- n_hit2/sum(hit2)
  p_fa2  <- n_fa2/sum(fa2)
  auc2_l <- p_hit2[1] * p_fa2[1] / 2
  for (n in 1:(cl - 1)) {
    auc2_l <- auc2_l + (p_hit2[n] + p_hit2[n + 1]) * (p_fa2[n + 1] - p_fa2[n]) / 2
  }
  
  print(as.numeric(c(fit_c[1, 1], fit_c[1, 3], fit_r[1, 3], fit_l[1, 3], 
                     auc1_c, auc1_r, auc1_l, auc2_c, auc2_r, auc2_l, d$Subj_idx[1])))
  
}

df2 <- foreach(i = 1:length(unique(dat$Subj_idx)), .combine = "rbind", .packages = c("dplyr", "metaSDT")) %dopar% {
  
  d <- subset(dat, dat$Subj_idx == unique(dat$Subj_idx)[i])
  d <- filter(d, row_number() %% 2 == 0)
  d <- na.omit(d)
  d <- mutate(d, Correct = ifelse(Stimulus == Response, 1, 0))
  
  # confidence
  nR_S1 <- c(sum(d$Stimulus == "1" & d$Response == "1" & d$Confidence == 1.0, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$Confidence == 0.9, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$Confidence == 0.8, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$Confidence == 0.7, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$Confidence == 0.6, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$Confidence == 0.5, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$Confidence == 0.5, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$Confidence == 0.6, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$Confidence == 0.7, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$Confidence == 0.8, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$Confidence == 0.9, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$Confidence == 1.0, na.rm = TRUE))
  
  nR_S2 <- c(sum(d$Stimulus == "2" & d$Response == "1" & d$Confidence == 1.0, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$Confidence == 0.9, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$Confidence == 0.8, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$Confidence == 0.7, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$Confidence == 0.6, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$Confidence == 0.5, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$Confidence == 0.5, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$Confidence == 0.6, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$Confidence == 0.7, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$Confidence == 0.8, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$Confidence == 0.9, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$Confidence == 1.0, na.rm = TRUE))
  
  fit_c <- data.frame(NA, NA, NA)
  tryCatch({
    fit_c <- fit_meta_d_MLE(nR_S1, nR_S2)
  }, error = function(e) {
    cat("Error: Fitting the model failed with message:", e$message, "\n")
  })
  
  n_hit <- cumsum(nR_S1)
  n_fa  <- cumsum(nR_S2)
  p_hit <- n_hit/sum(nR_S1)
  p_fa  <- n_fa/sum(nR_S2)
  auc1_c <- p_hit[1] * p_fa[1] / 2
  for (n in 1:(length(nR_S1) - 1)) {
    auc1_c <- auc1_c + (p_hit[n] + p_hit[n + 1]) * (p_fa[n + 1] - p_fa[n]) / 2
  }
  
  cl <- length(nR_S1) / 2
  hit2 <- nR_S1[1:cl] + rev(nR_S2[(cl + 1):(cl * 2)])
  fa2  <- rev(nR_S1[(cl + 1):(cl * 2)]) + nR_S2[1:cl]
  n_hit2 <- cumsum(hit2)
  n_fa2  <- cumsum(fa2)
  p_hit2 <- n_hit2/sum(hit2)
  p_fa2  <- n_fa2/sum(fa2)
  auc2_c <- p_hit2[1] * p_fa2[1] / 2
  for (n in 1:(cl - 1)) {
    auc2_c <- auc2_c + (p_hit2[n] + p_hit2[n + 1]) * (p_fa2[n + 1] - p_fa2[n]) / 2
  }
  
  # RT
  q <- quantile(d$RT_dec, probs = seq(0, 1, by = 1/6), na.rm = TRUE)
  d$RT_bin <- cut(d$RT_dec, breaks = q, labels = FALSE, include.lowest = TRUE)
  
  nR_S1 <- c(sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 5, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 6, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$RT_bin == 6, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$RT_bin == 5, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$RT_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$RT_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$RT_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$RT_bin == 1, na.rm = TRUE))
  
  nR_S2 <- c(sum(d$Stimulus == "2" & d$Response == "1" & d$RT_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$RT_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$RT_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$RT_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$RT_bin == 5, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$RT_bin == 6, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$RT_bin == 6, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$RT_bin == 5, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$RT_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$RT_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$RT_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$RT_bin == 1, na.rm = TRUE))
  
  fit_r <- data.frame(NA, NA, NA)
  tryCatch({
    fit_r <- fit_meta_d_MLE(nR_S1, nR_S2)
  }, error = function(e) {
    cat("Error: Fitting the model failed with message:", e$message, "\n")
  })
  
  n_hit <- cumsum(nR_S1)
  n_fa  <- cumsum(nR_S2)
  p_hit <- n_hit/sum(nR_S1)
  p_fa  <- n_fa/sum(nR_S2)
  auc1_r <- p_hit[1] * p_fa[1] / 2
  for (n in 1:(length(nR_S1) - 1)) {
    auc1_r <- auc1_r + (p_hit[n] + p_hit[n + 1]) * (p_fa[n + 1] - p_fa[n]) / 2
  }
  
  cl <- length(nR_S1) / 2
  hit2 <- nR_S1[1:cl] + rev(nR_S2[(cl + 1):(cl * 2)])
  fa2  <- rev(nR_S1[(cl + 1):(cl * 2)]) + nR_S2[1:cl]
  n_hit2 <- cumsum(hit2)
  n_fa2  <- cumsum(fa2)
  p_hit2 <- n_hit2/sum(hit2)
  p_fa2  <- n_fa2/sum(fa2)
  auc2_r <- p_hit2[1] * p_fa2[1] / 2
  for (n in 1:(cl - 1)) {
    auc2_r <- auc2_r + (p_hit2[n] + p_hit2[n + 1]) * (p_fa2[n + 1] - p_fa2[n]) / 2
  }
  
  # logit
  m <- glm(Correct ~ Confidence * RT_dec, family = binomial, data = d)
  d$logit <- predict(m)
  q <- quantile(d$logit, probs = seq(0, 1, by = 1/6), na.rm = TRUE)
  d$logit_bin <- cut(d$logit, breaks = q, labels = FALSE, include.lowest = TRUE)
  
  nR_S1 <- c(sum(d$Stimulus == "1" & d$Response == "1" & d$logit_bin == 6, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$logit_bin == 5, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$logit_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$logit_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$logit_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$logit_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$logit_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$logit_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$logit_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$logit_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$logit_bin == 5, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "2" & d$logit_bin == 6, na.rm = TRUE))
  
  nR_S2 <- c(sum(d$Stimulus == "2" & d$Response == "1" & d$logit_bin == 6, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$logit_bin == 5, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$logit_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$logit_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$logit_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "1" & d$logit_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$logit_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$logit_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$logit_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$logit_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$logit_bin == 5, na.rm = TRUE),
             sum(d$Stimulus == "2" & d$Response == "2" & d$logit_bin == 6, na.rm = TRUE))
  
  fit_l <- data.frame(NA, NA, NA)
  tryCatch({
    fit_l <- fit_meta_d_MLE(nR_S1, nR_S2)
  }, error = function(e) {
    cat("Error: Fitting the model failed with message:", e$message, "\n")
  })
  
  n_hit <- cumsum(nR_S1)
  n_fa  <- cumsum(nR_S2)
  p_hit <- n_hit/sum(nR_S1)
  p_fa  <- n_fa/sum(nR_S2)
  auc1_l <- p_hit[1] * p_fa[1] / 2
  for (n in 1:(length(nR_S1) - 1)) {
    auc1_l <- auc1_l + (p_hit[n] + p_hit[n + 1]) * (p_fa[n + 1] - p_fa[n]) / 2
  }
  
  cl <- length(nR_S1) / 2
  hit2 <- nR_S1[1:cl] + rev(nR_S2[(cl + 1):(cl * 2)])
  fa2  <- rev(nR_S1[(cl + 1):(cl * 2)]) + nR_S2[1:cl]
  n_hit2 <- cumsum(hit2)
  n_fa2  <- cumsum(fa2)
  p_hit2 <- n_hit2/sum(hit2)
  p_fa2  <- n_fa2/sum(fa2)
  auc2_l <- p_hit2[1] * p_fa2[1] / 2
  for (n in 1:(cl - 1)) {
    auc2_l <- auc2_l + (p_hit2[n] + p_hit2[n + 1]) * (p_fa2[n + 1] - p_fa2[n]) / 2
  }
  
  print(as.numeric(c(fit_c[1, 1], fit_c[1, 3], fit_r[1, 3], fit_l[1, 3], 
                     auc1_c, auc1_r, auc1_l, auc2_c, auc2_r, auc2_l, d$Subj_idx[1])))
  
}

df1 <- as.data.frame(df1)
df1$rn <- "odd"
df2 <- as.data.frame(df2)
df2$rn <- "even"
df <- rbind(df1, df2)

colnames(df) <- c("dp", "mdp_conf", "mdp_rt", "mdp_logit", 
                  "auc1_conf", "auc1_rt", "auc1_logit",
                  "auc2_conf", "auc2_rt", "auc2_logit", "sub", "rn")

df <- as.data.frame(df)
df <- na.omit(df)
summary(df)