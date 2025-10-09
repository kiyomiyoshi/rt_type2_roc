library(metaSDT)
library(tidyverse)
library(GGally)
library(data.table)
library(car)
library(doParallel)
library(ppcor)

cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)

dat <- fread("data_Maniscalco_2017_expt3.csv", header = T)
dat <- filter(dat, Subj_idx != 14) 

df <- foreach(i = 1:length(unique(dat$Subj_idx)), .combine = "rbind", .packages = c("dplyr", "metaSDT")) %dopar% {
  
  d <- subset(dat, dat$Subj_idx == unique(dat$Subj_idx)[i])
  d <- na.omit(d)
  d <- mutate(d, Correct = ifelse(Stimulus == Response, 1, 0))
  
  # confidence
  nR_S1 <- c(sum(d$Stimulus == "1" & d$Response == "1" & d$Confidence == 4, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$Confidence == 3, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$Confidence == 2, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$Confidence == 1, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "0" & d$Confidence == 1, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "0" & d$Confidence == 2, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "0" & d$Confidence == 3, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "0" & d$Confidence == 4, na.rm = TRUE))
  
  nR_S2 <- c(sum(d$Stimulus == "0" & d$Response == "1" & d$Confidence == 4, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "1" & d$Confidence == 3, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "1" & d$Confidence == 2, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "1" & d$Confidence == 1, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "0" & d$Confidence == 1, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "0" & d$Confidence == 2, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "0" & d$Confidence == 3, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "0" & d$Confidence == 4, na.rm = TRUE))
  
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
  q <- quantile(d$RT_dec, probs = seq(0, 1, by = 0.25), na.rm = TRUE)
  d$RT_bin <- cut(d$RT_dec, breaks = q, labels = FALSE, include.lowest = TRUE)
  
  nR_S1 <- c(sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "0" & d$RT_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "0" & d$RT_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "0" & d$RT_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "0" & d$RT_bin == 1, na.rm = TRUE))
  
  nR_S2 <- c(sum(d$Stimulus == "0" & d$Response == "1" & d$RT_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "1" & d$RT_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "1" & d$RT_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "1" & d$RT_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "0" & d$RT_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "0" & d$RT_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "0" & d$RT_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "0" & d$RT_bin == 1, na.rm = TRUE))
  
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
  q <- quantile(d$logit, probs = seq(0, 1, by = 0.25), na.rm = TRUE)
  d$logit_bin <- cut(d$logit, breaks = q, labels = FALSE, include.lowest = TRUE)
  
  nR_S1 <- c(sum(d$Stimulus == "1" & d$Response == "1" & d$logit_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$logit_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$logit_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "1" & d$logit_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "0" & d$logit_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "0" & d$logit_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "0" & d$logit_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "1" & d$Response == "0" & d$logit_bin == 4, na.rm = TRUE))
  
  nR_S2 <- c(sum(d$Stimulus == "0" & d$Response == "1" & d$logit_bin == 4, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "1" & d$logit_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "1" & d$logit_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "1" & d$logit_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "0" & d$logit_bin == 1, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "0" & d$logit_bin == 2, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "0" & d$logit_bin == 3, na.rm = TRUE),
             sum(d$Stimulus == "0" & d$Response == "0" & d$logit_bin == 4, na.rm = TRUE))
  
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

df <- as.data.frame(df)
colnames(df) <- c("dp", "mdp_conf", "mdp_rt", "mdp_logit", 
                  "auc1_conf", "auc1_rt", "auc1_logit",
                  "auc2_conf", "auc2_rt", "auc2_logit", "Subj_idx")
df %>%
  dplyr::filter(
      dp        > -1 & dp        < 4 &
      mdp_conf  > -1 & mdp_conf  < 4 &
      mdp_rt    > -1 & mdp_rt    < 4 & 
      mdp_logit > -1 & mdp_logit < 4) -> df
nrow(df)

rt_means <- dat %>%
  group_by(Subj_idx) %>%
  summarise(mean_RT_Dec = mean(RT_dec, na.rm = TRUE)) %>%
  mutate(group = if_else(mean_RT_Dec <= median(mean_RT_Dec), "Fast", "Slow"))

df <- df %>%
  left_join(rt_means, by = "Subj_idx")

df <- na.omit(df)
summary(df)
t.test(df$mdp_conf, df$mdp_logit, paired = T)
ggpairs(df[, c("dp", "mdp_conf", "mdp_rt", "mdp_logit")]) + xlim(-2, 4)

# correlation coefficients
vars <- c("dp", "mdp_conf", "mdp_rt", "mdp_logit")
cor_mat <- cor(df[, vars], use = "pairwise.complete.obs")
r_vec <- cor_mat[lower.tri(cor_mat)]
pair_names <- combn(vars, 2, FUN = function(x) paste(x, collapse = "-"))
cor_df <- as.data.frame(t(r_vec))
colnames(cor_df) <- pair_names
cor_df$n <- nrow(df)
pcor <- pcor.test(df$mdp_conf, df$mdp_rt, df$dp)

# correlation coefficients for m-ratio
df %>%
  dplyr::filter(dp > 0 & mdp_conf > 0 & mdp_rt > 0 & mdp_logit > 0) %>%
  dplyr::mutate(mratio_conf =  mdp_conf / dp,
                mratio_rt =    mdp_rt / dp,
                mratio_logit = mdp_logit / dp) -> df_filtered

vars <- c("dp", "mratio_conf", "mratio_rt", "mratio_logit")
cor_mat <- cor(df_filtered[, vars], use = "pairwise.complete.obs")
r_vec <- cor_mat[lower.tri(cor_mat)]
pair_names <- combn(vars, 2, FUN = function(x) paste(x, collapse = "-"))
cor_df_mratio <- as.data.frame(t(r_vec))
colnames(cor_df_mratio) <- pair_names
cor_df_mratio$n <- nrow(df_filtered)

# correlation coefficients for m-diff
df %>%
  dplyr::mutate(mdiff_conf =   mdp_conf - dp,
                mdiff_rt =     mdp_rt - dp,
                mdiff_logit =  mdp_logit - dp) -> df_mutated

vars <- c("dp", "mdiff_conf", "mdiff_rt", "mdiff_logit")
cor_mat <- cor(df_mutated[, vars], use = "pairwise.complete.obs")
r_vec <- cor_mat[lower.tri(cor_mat)]
pair_names <- combn(vars, 2, FUN = function(x) paste(x, collapse = "-"))
cor_df_mdiff <- as.data.frame(t(r_vec))
colnames(cor_df_mdiff) <- pair_names
cor_df_mdiff$n <- nrow(df_mutated)