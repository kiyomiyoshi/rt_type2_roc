library(doParallel)
library(forcats)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tictoc)
library(cowplot)
library(sjPlot)
library(ggsci)
library(matrixStats)
library(metaSDT)
library(car)
library(scales)
library(ggeffects)
library(ggh4x)

source("theme_publication.R")
theme_set(theme_publication()) 
options(scipen = 6)

cores <- getOption("mc.cores", detectCores())
cl <- makeCluster(cores)
registerDoParallel(cl)


#'# simulation
nu_tar <- c(0.0015, 0.0025, 0.0035)
nu_dis <- 0.001
eta <- c(0.0005, 0.005)
s <- 0.03
sz <- c(0.1, 0.6)
tnd <- 0.3
a <- 1.2
trial <- 100000
sample <- 3501

df1 <- c()
df2 <- c()

tic()
for (i in nu_tar) {
  for (j in eta) {
    for (k in sz) {
      dr <- cbind(rnorm(trial, i,    j),
                  rnorm(trial, nu_dis, j))
      a_1 <- runif(trial, -k, k)
      a_2 <- runif(trial, -k, k)
      d <- foreach(l = 1:trial, .combine = "rbind", .packages = c("dplyr", "matrixStats")) %dopar% {
        accum <- matrix(0, nrow = sample + 1, ncol = 2)
        accum[, 1] <- c(a_1[l], rep(dr[l, 1], sample))
        accum[, 2] <- c(a_2[l], rep(dr[l, 2], sample))
        accum[-1, 1] <- accum[-1, 1] + rnorm(sample, 0, s)
        accum[-1, 2] <- accum[-1, 2] + rnorm(sample, 0, s)
        accum <- matrixStats::colCumsums(accum)
        mvn <- matrixStats::rowOrderStats(-accum, which = 2) - matrixStats::rowOrderStats(-accum, which = 1)
        simrt  <- which(mvn > a)[1]
        response <- which.max(accum[simrt, ])
        second <- which.max(-accum[simrt, ])
        e_chosen_post <- ifelse(simrt <= 3000, accum[simrt + 150, response], NA)
        e_second_post <- ifelse(simrt <= 3000, accum[simrt + 150, second], NA)
        as.numeric(c(simrt, response, mvn[simrt], mvn[1], mvn[350], e_chosen_post, e_second_post))
      }
      d <- as.data.frame(d)
      colnames(d) <- c("rt", "response", "mvn_response", "mvn_0", "mvn_350", "e_chosen_post", "e_second_post")
      d$rt   <- d$rt + tnd * 1000
      d$nu_tar <- i
      d$eta <- j
      d$sz <- k
      df1 <- rbind(df1, filter(d, rt <= 3000))
    }
  }
}

for (i in nu_tar) {
  for (j in eta) {
    for (k in sz) {
      dr <- cbind(rnorm(trial, nu_dis, j),
                  rnorm(trial, i,    j)) 
      a_1 <- runif(trial, -k, k)
      a_2 <- runif(trial, -k, k)
      d <- foreach(l = 1:trial, .combine = "rbind", .packages = c("dplyr", "matrixStats")) %dopar% {
        accum <- matrix(0, nrow = sample + 1, ncol = 2)
        accum[, 1] <- c(a_1[l], rep(dr[l, 1], sample))
        accum[, 2] <- c(a_2[l], rep(dr[l, 2], sample))
        accum[-1, 1] <- accum[-1, 1] + rnorm(sample, 0, s)
        accum[-1, 2] <- accum[-1, 2] + rnorm(sample, 0, s)
        accum <- matrixStats::colCumsums(accum)
        mvn <- matrixStats::rowOrderStats(-accum, which = 2) - matrixStats::rowOrderStats(-accum, which = 1)
        simrt  <- which(mvn > a)[1]
        response <- which.max(accum[simrt, ])
        second <- which.max(-accum[simrt, ])
        e_chosen_post <- ifelse(simrt <= 3000, accum[simrt + 150, response], NA)
        e_second_post <- ifelse(simrt <= 3000, accum[simrt + 150, second], NA)
        as.numeric(c(simrt, response, mvn[simrt], mvn[1], mvn[350], e_chosen_post, e_second_post))  
      }
      d <- as.data.frame(d)
      colnames(d) <- c("rt", "response", "mvn_response", "mvn_0", "mvn_350", "e_chosen_post", "e_second_post")
      d$rt   <- d$rt + tnd * 1000
      d$nu_tar <- i
      d$eta <- j
      d$sz <- k
      df2 <- rbind(df2, filter(d, rt <= 3000))
    }
  }
}
toc()


#'# preprocess
df1$stim <- "s1"
df2$stim <- "s2"
df1$choice <- ifelse(df1$response == 1, "target", "distractor")
df2$choice <- ifelse(df2$response == 2, "target", "distractor")
simdat_ddm <- rbind(df1, df2)
simdat_ddm$choice <- factor(simdat_ddm$choice, levels = c("target", "distractor")) # this is required
simdat_ddm$correct = ifelse(simdat_ddm$choice == "target", 1, 0)
simdat_ddm <- na.omit(simdat_ddm)
simdat_ddm <- mutate(simdat_ddm, early_conf = mvn_350 - mvn_0)
simdat_ddm <- mutate(simdat_ddm, post_conf  = e_chosen_post - e_second_post - mvn_response)

# logistic regression
simdf <- c()
coef <- list()
test <- list()

for (i in nu_tar) {
  for (j in eta) {
    for (k in sz) {
      df <- subset(simdat_ddm, simdat_ddm$nu_tar == i & simdat_ddm$eta == j & simdat_ddm$sz == k)
      m1 <- glm(correct ~ mvn_response * rt, family = binomial, data = df)
      m2 <- glm(correct ~ post_conf * rt, family = binomial, data = df)
      df$logit_be <- predict(m1)
      df$logit_post <- predict(m2)
      simdf <- rbind(simdf, df)
      coef <- append(coef, summary(m1))
      test <- append(test, Anova(m1))
    }
  }
}
simdat_ddm <- simdf

# roc
sdt_ddm <- c()
roc_1_ddm <- c()
roc_2_ddm <- c()
for (i in unique(simdat_ddm$nu_tar)) {
  for (j in unique(simdat_ddm$eta)) {
    for (k in unique(simdat_ddm$sz)) {
      
      sd <- subset(simdat_ddm, simdat_ddm$nu_tar == i & simdat_ddm$eta == j & simdat_ddm$sz == k)
      
      cq <- quantile(sd$post_conf, probs = seq(0, 1, 0.1))
      all_levels <- seq(length(cq) - 1)
      sd$conf_bin <- cut(sd$post_conf, breaks = cq, labels = FALSE)
      sd$conf_bin <- factor(sd$conf_bin, levels = all_levels)
      
      h1c <- c(rev(table(sd[which(sd$response == 2 & sd$stim == "s2"), "conf_bin"])),
               table(sd[which(sd$response == 1 & sd$stim == "s2"), "conf_bin"]))
      f1c <- c(rev(table(sd[which(sd$response == 2 & sd$stim == "s1"), "conf_bin"])),
               table(sd[which(sd$response == 1 & sd$stim == "s1"), "conf_bin"]))
      dfc_1 <- data.frame(hr_1 = cumsum(h1c)/sum(h1c), far_1 = cumsum(f1c)/sum(f1c))
      dfc_1$variable <- "post"
      sdt_ddm_c <- fit_meta_d_MLE(h1c, f1c)
      
      h2c <- rev(table(sd[which(sd$choice == "target"), "conf_bin"]))
      f2c <- rev(table(sd[which(sd$choice == "distractor"), "conf_bin"]))
      dfc_2 <- data.frame(hr_2 = cumsum(h2c)/sum(h2c), far_2 = cumsum(f2c)/sum(f2c))
      dfc_2$variable <- "post"
      
      bq <- quantile(sd$mvn_response, probs = seq(0, 1, 0.1))
      all_levels <- seq(length(bq) - 1)
      sd$be_bin <- cut(sd$mvn_response, breaks = bq, labels = FALSE)
      sd$be_bin <- factor(sd$be_bin, levels = all_levels)
      
      h1b <- c(rev(table(sd[which(sd$response == 2 & sd$stim == "s2"), "be_bin"])),
               table(sd[which(sd$response == 1 & sd$stim == "s2"), "be_bin"]))
      f1b <- c(rev(table(sd[which(sd$response == 2 & sd$stim == "s1"), "be_bin"])),
               table(sd[which(sd$response == 1 & sd$stim == "s1"), "be_bin"]))
      dfb_1 <- data.frame(hr_1 = cumsum(h1b)/sum(h1b), far_1 = cumsum(f1b)/sum(f1b))
      dfb_1$variable <- "be"
      sdt_ddm_b <- fit_meta_d_MLE(h1b, f1b)
      
      h2b <- rev(table(sd[which(sd$choice == "target"), "be_bin"]))
      f2b <- rev(table(sd[which(sd$choice == "distractor"), "be_bin"]))
      dfb_2 <- data.frame(hr_2 = cumsum(h2b)/sum(h2b), far_2 = cumsum(f2b)/sum(f2b))
      dfb_2$variable <- "be"
      
      rq <- quantile(sd$rt, probs = seq(0, 1, 0.1))
      sd$rt_bin <- cut(sd$rt, breaks = rq, labels = FALSE)
      sd$rt_bin <- factor(sd$rt_bin, levels = all_levels)
      
      h1r <- c(table(sd[which(sd$response == 2 & sd$stim == "s2"), "rt_bin"]),
               rev(table(sd[which(sd$response == 1 & sd$stim == "s2"), "rt_bin"])))
      f1r <- c(table(sd[which(sd$response == 2 & sd$stim == "s1"), "rt_bin"]),
               rev(table(sd[which(sd$response == 1 & sd$stim == "s1"), "rt_bin"])))
      dfr_1 <- data.frame(hr_1 = cumsum(h1r)/sum(h1r), far_1 = cumsum(f1r)/sum(f1r))
      dfr_1$variable <- "RT"
      sdt_ddm_r <- fit_meta_d_MLE(h1r, f1r)
      
      h2r <- table(sd[which(sd$choice == "target"), "rt_bin"])
      f2r <- table(sd[which(sd$choice == "distractor"), "rt_bin"])
      dfr_2 <- data.frame(hr_2 = cumsum(h2r)/sum(h2r), far_2 = cumsum(f2r)/sum(f2r))
      dfr_2$variable <- "RT"
      
      lbq <- quantile(sd$logit_be, probs = seq(0, 1, 0.1))
      sd$logit_be_bin <- cut(sd$logit_be, breaks = lbq, labels = FALSE)
      sd$logit_be_bin <- factor(sd$logit_be_bin, levels = all_levels)
      
      h1lb <- c(rev(table(sd[which(sd$response == 2 & sd$stim == "s2"), "logit_be_bin"])),
                table(sd[which(sd$response == 1 & sd$stim == "s2"), "logit_be_bin"]))
      f1lb <- c(rev(table(sd[which(sd$response == 2 & sd$stim == "s1"), "logit_be_bin"])),
                table(sd[which(sd$response == 1 & sd$stim == "s1"), "logit_be_bin"]))
      dflb_1 <- data.frame(hr_1 = cumsum(h1lb)/sum(h1lb), far_1 = cumsum(f1lb)/sum(f1lb))
      dflb_1$variable <- "logit_be"
      sdt_ddm_lb <- fit_meta_d_MLE(h1lb, f1lb)
      
      h2lb <- rev(table(sd[which(sd$choice == "target"), "logit_be_bin"]))
      f2lb <- rev(table(sd[which(sd$choice == "distractor"), "logit_be_bin"]))
      dflb_2 <- data.frame(hr_2 = cumsum(h2lb)/sum(h2lb), far_2 = cumsum(f2lb)/sum(f2lb))
      dflb_2$variable <- "logit_be"
      
      lpq <- quantile(sd$logit_post, probs = seq(0, 1, 0.1))
      sd$logit_po_bin <- cut(sd$logit_post, breaks = lpq, labels = FALSE)
      sd$logit_po_bin <- factor(sd$logit_po_bin, levels = all_levels)
      
      h1lp <- c(rev(table(sd[which(sd$response == 2 & sd$stim == "s2"), "logit_po_bin"])),
                table(sd[which(sd$response == 1 & sd$stim == "s2"), "logit_po_bin"]))
      f1lp <- c(rev(table(sd[which(sd$response == 2 & sd$stim == "s1"), "logit_po_bin"])),
                table(sd[which(sd$response == 1 & sd$stim == "s1"), "logit_po_bin"]))
      dflp_1 <- data.frame(hr_1 = cumsum(h1lp)/sum(h1lp), far_1 = cumsum(f1lp)/sum(f1lp))
      dflp_1$variable <- "logit_post"
      sdt_ddm_lp <- fit_meta_d_MLE(h1lp, f1lp)
      
      h2lp <- rev(table(sd[which(sd$choice == "target"), "logit_po_bin"]))
      f2lp <- rev(table(sd[which(sd$choice == "distractor"), "logit_po_bin"]))
      dflp_2 <- data.frame(hr_2 = cumsum(h2lp)/sum(h2lp), far_2 = cumsum(f2lp)/sum(f2lp))
      dflp_2$variable <- "logit_post"
      
      df_1 <- rbind(dfc_1, dfb_1, dfr_1, dflb_1, dflp_1)
      df_2 <- rbind(dfc_2, dfb_2, dfr_2, dflb_2, dflp_2)
      df_1$nu_tar <- i
      df_2$nu_tar <- i
      df_1$eta <- j
      df_2$eta <- j
      df_1$sz <- k
      df_2$sz <- k
      
      sdt_ddm <- rbind(sdt_ddm, cbind(sdt_ddm_b[1, 1], sdt_ddm_b[1, 3], sdt_ddm_c[1, 3], sdt_ddm_r[1, 3], sdt_ddm_lb[1, 3], sdt_ddm_lp[1, 3], i, j, k))
      roc_1_ddm <- rbind(roc_1_ddm, df_1)
      roc_2_ddm <- rbind(roc_2_ddm, df_2)
    }
  }
}

# sdt
colnames(sdt_ddm) <- c("dp", "mdp_be", "mdp_post", "mdp_rt", "mdp_logit_be", "mdp_logit_post", "nu_tar", "eta", "sz")
sdt_ddm <- as.data.frame(sdt_ddm)
sdt_ddm <- mutate(sdt_ddm, 
              mratio_be = mdp_be / dp,
              mratio_post = mdp_post / dp,
              mratio_rt = mdp_rt / dp,
              mratio_logit_be = mdp_logit_be / dp,
              mratio_logit_post = mdp_logit_post / dp,
              be_logit_to_be = mdp_logit_be / mdp_be,
              po_logit_to_po = mdp_logit_post / mdp_post)

roc_1_sdt_ddm <- c()
roc_2_sdt_ddm <- c()
for (l in 1:nrow(sdt_ddm)) {
  dp <- sdt_ddm[l, 1]
  x2 <- rnorm(10000, dp, 1)
  x1 <- rnorm(10000, 0,  1)
  q <- quantile(c(x1, x2), probs = seq(0, 1, 0.02))
  all_levels <- seq(length(q) - 1)
  
  resp1 <- cut(x1, breaks = q, labels = FALSE)
  resp1 <- factor(resp1, levels = all_levels)
  resp1 <- rev(table(resp1))
  resp2 <- cut(x2, breaks = q, labels = FALSE)
  resp2 <- factor(resp2, levels = all_levels)
  resp2 <- rev(table(resp2))
  
  h2 <- resp2[1:(length(resp2) / 2)] + rev(resp1[(length(resp1) / 2 + 1):length(resp1)])
  f2 <- rev(resp2[(length(resp2) / 2 + 1):length(resp2)]) + resp1[1:(length(resp1) / 2)]
  
  roc_1_sdt_ddm <- rbind(roc_1_sdt_ddm, data.frame(hr_1 = cumsum(resp2)/sum(resp2), far_1 = cumsum(resp1)/sum(resp1), 
                                           nu_tar = sdt_ddm$nu_tar[l], eta = sdt_ddm$eta[l], sz = sdt_ddm$sz[l]))
  roc_2_sdt_ddm <- rbind(roc_2_sdt_ddm, data.frame(hr_2 = cumsum(h2)/sum(h2), far_2 = cumsum(f2)/sum(f2),
                                           nu_tar = sdt_ddm$nu_tar[l], eta = sdt_ddm$eta[l], sz = sdt_ddm$sz[l]))
}


#'# visualization
# choice proportion
simdat_ddm %>%
  group_by(choice, nu_tar, eta, sz) %>%
  mutate(n = n()) %>%
  group_by(nu_tar, eta, sz) %>%
  mutate(p = ifelse(choice == "target", n[choice == "target"] / n(),
                    n[choice == "distractor"] / n())) %>%
  distinct(nu_tar, eta, sz, choice, n, p) %>%
  mutate(condition = ifelse(eta == 5e-04 & sz == 0.1, "η == 0.0005 ~ phantom(',') ~ s[z] == 0.1",
                            ifelse(eta == 5e-04 & sz == 0.6, "η == 0.0005 ~ phantom(',') ~ s[z] == 0.6",
                                   ifelse(eta == 0.005 & sz == 0.1, "η == 0.005 ~ phantom(',') ~ s[z] == 0.1",
                                          "η == 0.005 ~ phantom(',') ~ s[z] == 0.6")))) %>%
  ggplot() + geom_line(aes(x = nu_tar, y = p, color = choice), size = 0.3) + 
  geom_point(aes(x = nu_tar, y = p, color = factor(choice)), size = 0.3) + 
  ylab("Choice proportion") +
  scale_color_discrete(labels = c("target", "distractor")) + labs(color = "choice") +
  scale_x_continuous(breaks = c(0.0015, 0.0025, 0.0035)) + xlab(expression(bold(paste({ν[target]})))) +
  theme(legend.text = element_text(size = 7),
        legend.position = "top", 
        legend.direction = "horizontal",
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = NA)) +
  facet_wrap(~ condition, labeller = label_parsed, nrow = 2)  +
  scale_color_npg() + coord_cartesian(ylim = c(0, 1)) -> p1_ddm
p1_ddm
p1_ddm <- p1_ddm + guides(color = F) + theme(legend.text = element_text(size = 3),
                                             legend.position = "top", 
                                             legend.direction = "horizontal",
                                             legend.background = element_rect(fill = NA, colour = NA),
                                             legend.key = element_rect(fill = NA),
                                             strip.text = element_text(size = 5.8, face = "bold", color = "black"),
                                             plot.margin = margin(5.5, 5.5, 5.5, 5.5), 
                                             plot.title = element_text(size = 5), 
                                             axis.title.x = element_text(size = 7),
                                             axis.title.y = element_text(size = 7),
                                             axis.text.x = element_text(size = 5, angle = 25),
                                             axis.text.y = element_text(size = 5))

# RT histogram
simdat_ddm %>%
  complete(choice, nu_tar) %>%
  ggplot() +
  geom_histogram(aes(x = rt, color = choice, fill = choice), alpha = 0.5) +
  scale_color_discrete(labels = c("target", "distractor")) +
  xlab("RT (s)") + ylab("Count") + coord_cartesian(ylim = c(0, 120000)) +
  scale_x_continuous(breaks = c(0, 1000, 2000, 3000), labels = c("0", "1", "2", "3")) +
  facet_wrap(. ~ nu_tar + eta + sz, ncol = 4,
             labeller = labeller(nu_tar = as_labeller(c("0.0015" = "ν[target] == 0.0015", "0.0025" = "ν[target] == 0.0025", 
                                                           "0.0035" = "ν[target] == 0.0035"), label_parsed),
                                 eta = as_labeller(c("0.0005" = "η == 0.0005", "0.005" = "η == 0.005"), label_parsed),
                                 sz = as_labeller(c("0.1" = "s[z] == 0.1", "0.6" = "s[z] == 0.6"), label_parsed))) +
  scale_fill_npg() + scale_color_npg() + 
  theme(axis.title = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(-0.5, -0.5, -0.5, -0.5, "pt"))) -> p2_ddm
p2_ddm 
p2_ddm <- p2_ddm + 
  theme(strip.text = element_text(size = 6.5, face = "bold"),
        plot.margin = margin(5, 5, 5, 5), 
        plot.title = element_text(size = 8))

# mean RT
simdat_ddm %>%
  group_by(choice, nu_tar, eta, sz) %>%
  summarise(mean_rt = mean(rt)) %>%
  mutate(condition = ifelse(eta == 5e-04 & sz == 0.1, "η == 0.0005 ~ phantom(',') ~ s[z] == 0.1",
                            ifelse(eta == 5e-04 & sz == 0.6, "η == 0.0005 ~ phantom(',') ~ s[z] == 0.6",
                                   ifelse(eta == 0.005 & sz == 0.1, "η == 0.005 ~ phantom(',') ~ s[z] == 0.1",
                                          "η == 0.005 ~ phantom(',') ~ s[z] == 0.6")))) %>%
  ggplot() + geom_line(aes(x = nu_tar, y = mean_rt / 1000, color = factor(choice)), size = 0.32) +
  geom_point(aes(x = nu_tar, y = mean_rt / 1000, color = factor(choice)), size = 0.32) +
  scale_color_discrete(labels = c("target", "distractor")) + labs(color = "choice") +
  scale_x_continuous(breaks = c(0.0015, 0.0025, 0.0035)) + xlab(expression(bold(paste({ν[target]})))) +
  facet_wrap(~ condition, labeller = label_parsed, nrow = 2)  +
  ylab("Mean RT (s)") + scale_color_npg() + coord_cartesian(ylim = c(0.5, 1.5)) -> p3_ddm
p3_ddm
p3_ddm <- p3_ddm + guides(color = F) + theme(legend.text = element_text(size = 3),
                                             legend.position = "top", 
                                             legend.direction = "horizontal",
                                             legend.background = element_rect(fill = NA, colour = NA),
                                             legend.key = element_rect(fill = NA),
                                             strip.text = element_text(size = 5.8, face = "bold", color = "black"),
                                             plot.margin = margin(5.5, 5.5, 5.5, 5.5), 
                                             plot.title = element_text(size = 5), 
                                             axis.title.x = element_text(size = 7),
                                             axis.title.y = element_text(size = 7),
                                             axis.text.x = element_text(size = 5, angle = 25),
                                             axis.text.y = element_text(size = 5))

# mean post confidence
simdat_ddm %>%
  group_by(choice, nu_tar, eta, sz) %>%
  mutate(mean_conf = mean(post_conf)) %>%
  mutate(condition = ifelse(eta == 5e-04 & sz == 0.1, "η == 0.0005 ~ phantom(',') ~ s[z] == 0.1",
                            ifelse(eta == 5e-04 & sz == 0.6, "η == 0.0005 ~ phantom(',') ~ s[z] == 0.6",
                                   ifelse(eta == 0.005 & sz == 0.1, "η == 0.005 ~ phantom(',') ~ s[z] == 0.1",
                                          "η == 0.005 ~ phantom(',') ~ s[z] == 0.6")))) %>%
  ggplot() + geom_line(aes(x = nu_tar, y = mean_conf, color = factor(choice)), size = 0.3) +
  geom_point(aes(x = nu_tar, y = mean_conf, color = factor(choice)), size = 0.3) +
  scale_color_discrete(labels = c("target", "distractor")) + labs(color = "choice") +
  scale_x_continuous(breaks = c(0.0015, 0.0025, 0.0035)) + xlab(expression(bold(paste({ν[target]})))) +
  facet_wrap(~ condition, labeller = label_parsed, nrow = 2)  +
  ylab("Mean confidence") + scale_color_npg() + coord_cartesian(ylim = c(-0.5, 1.2)) -> p4_ddm
p4_ddm
p4_ddm <- p4_ddm + guides(color = F) + theme(legend.text = element_text(size = 3),
                                             legend.position = "top", 
                                             legend.direction = "horizontal",
                                             legend.background = element_rect(fill = NA, colour = NA),
                                             legend.key = element_rect(fill = NA),
                                             strip.text = element_text(size = 5.8, face = "bold", color = "black"),
                                             plot.margin = margin(5.5, 5.5, 5.5, 5.5), 
                                             plot.title = element_text(size = 5), 
                                             axis.title.x = element_text(size = 7),
                                             axis.title.y = element_text(size = 7),
                                             axis.text.x = element_text(size = 5, angle = 25),
                                             axis.text.y = element_text(size = 5))

# sdt measures
sdt_ddm %>%
  pivot_longer(cols = c("dp", "mdp_post", "mdp_logit_post", "mdp_rt"), 
               names_to = "index", values_to = "value") %>%
  mutate(index = fct_recode(index, "d'" = "dp", "meta-d' confidence" = "mdp_post", 
                            "meta-d' RT" = "mdp_rt", "meta-d' logit" = "mdp_logit_post"),
         index = fct_relevel(index, "d'", "meta-d' confidence", "meta-d' RT", "meta-d' logit")) %>%
  mutate(condition = ifelse(eta == 5e-04 & sz == 0.1, "η == 0.0005 ~ phantom(',') ~ s[z] == 0.1",
                            ifelse(eta == 5e-04 & sz == 0.6, "η == 0.0005 ~ phantom(',') ~ s[z] == 0.6",
                                   ifelse(eta == 0.005 & sz == 0.1, "η == 0.005 ~ phantom(',') ~ s[z] == 0.1",
                                          "η == 0.005 ~ phantom(',') ~ s[z] == 0.6")))) %>%
  ggplot() +  geom_point(aes(x = nu_tar, y = value, color = index, size = index), show.legend = FALSE) + 
  scale_size_manual(values = c("d'" = 0.3, "meta-d' confidence" = 0.8, "meta-d' RT" = 0.3, "meta-d' logit" = 0.3)) + labs(color = "Index") +
  geom_line(aes(x = nu_tar, y = value, color = index, linetype = ifelse(index == "d'", "solid", "dashed")), size = 0.3, show.legend = FALSE) + 
  theme(legend.position = "top") + ylab("meta-d'") +
  scale_color_manual(values = c("black", hue_pal()(3))) +
  scale_x_continuous(breaks = c(0.0015, 0.0025, 0.0035)) + xlab(expression(bold(paste({ν[target]})))) +
  facet_wrap(~ condition, labeller = label_parsed, nrow = 2, scales = "free_y") -> p5_ddm
p5_ddm

position_scales <- list(
  scale_y_continuous(limits = c(-1, 4.2), breaks = -1:4),
  scale_y_continuous(limits = c(-1, 4.2), breaks = -1:4),
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)),
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)))

p5_ddm <- p5_ddm + facetted_pos_scales(y = position_scales) + 
  guides(color = F) +
  theme(strip.text = element_text(size = 5.8, face = "bold", color = "black"),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5), 
        plot.title = element_text(size = 5), 
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 5, angle = 25),
        axis.text.y = element_text(size = 5))

# m-ratio
sdt_ddm %>%
  pivot_longer(cols = c("mratio_post", "mratio_rt", "mratio_logit_post"), 
               names_to = "index", values_to = "value") %>%
  mutate(index = fct_recode(index, "m-ratio confidence" = "mratio_post", "m-ratio RT" = "mratio_rt", "m-ratio logit" = "mratio_logit_post"),
         index = fct_relevel(index, "m-ratio confidence", "m-ratio RT", "m-ratio logit")) %>%
  mutate(condition = ifelse(eta == 5e-04 & sz == 0.1, "η == 0.0005 ~ phantom(',') ~ s[z] == 0.1",
                            ifelse(eta == 5e-04 & sz == 0.6, "η == 0.0005 ~ phantom(',') ~ s[z] == 0.6",
                                   ifelse(eta == 0.005 & sz == 0.1, "η == 0.005 ~ phantom(',') ~ s[z] == 0.1",
                                          "η == 0.005 ~ phantom(',') ~ s[z] == 0.6")))) %>%
  ggplot() + geom_hline(yintercept = 1, linetype = "dashed", size = 0.3) +
  geom_point(aes(x = nu_tar, y = value, color = index, size = index), show.legend = FALSE) + 
  scale_size_manual(values = c("m-ratio confidence" = 0.8, "m-ratio RT" = 0.3, "m-ratio logit" = 0.3)) + labs(color = "Index") +
  geom_line(aes(x = nu_tar, y = value, color = index), size = 0.3) +
  scale_color_manual(values = c(hue_pal()(3))) + ylab("m-ratio") +
  theme(legend.position = "top") +
  scale_x_continuous(breaks = c(0.0015, 0.0025, 0.0035)) + xlab(expression(bold(paste({ν[target]})))) +
  facet_wrap(~ condition, labeller = label_parsed, nrow = 2) +
  coord_cartesian(ylim = c(-0.3, 1.1)) + scale_y_continuous(breaks = seq(-0.2, 1, 0.4)) -> p6_ddm
p6_ddm
p6_ddm <- p6_ddm + guides(color = F) +
  theme(strip.text = element_text(size = 5.8, face = "bold", color = "black"),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5), 
        plot.title = element_text(size = 5), 
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 5, angle = 25),
        axis.text.y = element_text(size = 5))

# m-ratio boost by rt
sdt_ddm %>%
  pivot_longer(cols = c("po_logit_to_po"), 
               names_to = "index", values_to = "value") %>%
  mutate(condition = ifelse(eta == 5e-04 & sz == 0.1, "η == 0.0005 ~ phantom(',') ~ s[z] == 0.1",
                            ifelse(eta == 5e-04 & sz == 0.6, "η == 0.0005 ~ phantom(',') ~ s[z] == 0.6",
                                   ifelse(eta == 0.005 & sz == 0.1, "η == 0.005 ~ phantom(',') ~ s[z] == 0.1",
                                          "η == 0.005 ~ phantom(',') ~ s[z] == 0.6")))) %>%
  ggplot() + geom_point(aes(x = nu_tar, y = value, color = index), size = 0.3) +
  geom_line(aes(x = nu_tar, y = value, color = index), size = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 0.3) + theme(legend.position = "top") +
  scale_color_manual(values = c("#E36414")) +
  ylab(expression(bold(paste(meta, "-", d, "'"[confidence+RT]/meta, "-", d, "'"[confidence])))) +
  scale_x_continuous(breaks = c(0.0015, 0.0025, 0.0035)) + xlab(expression(bold(paste({ν[target]})))) +
  facet_wrap(~ condition, labeller = label_parsed, nrow = 2) +
  coord_cartesian(ylim = c(0.95, 1.8)) + scale_y_continuous(breaks = seq(1, 1.8, 0.2)) -> p7_ddm
p7_ddm
p7_ddm <- p7_ddm + guides(color = F) +
  theme(strip.text = element_text(size = 5.8, face = "bold", color = "black"),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5), 
        plot.title = element_text(size = 5), 
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 6.3),
        axis.text.x = element_text(size = 5, angle = 25),
        axis.text.y = element_text(size = 5))

# correctness x rt x post confidence
simdat_ddm$rt_bin <- cut(simdat_ddm$rt, breaks = quantile(simdat_ddm$rt, probs = seq(0, 1, by = 0.25)), include.lowest = TRUE)
levels(simdat_ddm$rt_bin) <- c("q1", "q2", "q3", "q4")
simdat_ddm %>%
  group_by(choice, nu_tar, eta, sz, rt_bin) %>%
  summarise(mean_conf = mean(post_conf), n = n()) %>%
  ggplot() + 
  geom_line(aes(x = rt_bin, y = mean_conf, color = choice, group = choice), size = 0.3) +
  geom_point(aes(x = rt_bin, y = mean_conf, color = choice), size = 0.7) +
  facet_wrap(~ nu_tar + eta + sz, ncol = 4,
             labeller = labeller(nu_tar = as_labeller(c("0.0015" = "ν[target] == 0.0015", "0.0025" = "ν[target] == 0.0025", "0.0035" = "ν[target] == 0.0035"), label_parsed),
                                 eta =  as_labeller(c("0.0005" = "η == 0.0005", "0.005" = "η == 0.005"), label_parsed),
                                 sz =  as_labeller(c("0.1" = "s[z] == 0.1", "0.6" = "s[z] == 0.6"), label_parsed))) +
  xlab("RT quantile") + ylab("Mean confidence") + scale_color_npg(guide = guide_legend(override.aes = list(color = "white"))) +
  theme(axis.title = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(-0.5, -0.5, -0.5, -0.5, "pt"))) -> p8_ddm
p8_ddm
p8_ddm <- p8_ddm + theme(strip.text = element_text(size = 6.5, face = "bold"),
                         plot.margin = margin(5, 5, 5, 5), 
                         plot.title = element_text(size = 8),
                         legend.text = element_text(color = "white"),
                         legend.title = element_text(color = "white"),
                         legend.key = element_rect(fill = "white"))

# accuracy x rt x post confidence
simdat_ddm$rt_bin <- cut(simdat_ddm$rt, breaks = quantile(simdat_ddm$rt, probs = seq(0, 1, 0.25)), include.lowest = TRUE)
levels(simdat_ddm$rt_bin) <- c("q1", "q2", "q3", "q4")
simdat_ddm$conf_bin <- cut(simdat_ddm$post_conf, breaks = quantile(simdat_ddm$post_conf, probs = seq(0, 1, 0.25)), include.lowest = TRUE)
levels(simdat_ddm$conf_bin) <- c("q1", "q2", "q3", "q4")
simdat_ddm %>%
  group_by(nu_tar, eta, sz, rt_bin, conf_bin) %>%
  summarise(accuracy = mean(correct)) %>%
  ggplot() + 
  geom_line(aes(x = conf_bin, y = accuracy, color = rt_bin, group = rt_bin), size = 0.3) +
  geom_point(aes(x = conf_bin, y = accuracy, color = rt_bin), size = 0.7) +
  scale_color_manual(values = c("#3300CC", "#7900FF", "#FF0077", "#CC0033")) +
  facet_wrap(~ nu_tar + eta + sz, ncol = 4,
             labeller = labeller(nu_tar = as_labeller(c("0.0015" = "ν[target] == 0.0015", "0.0025" = "ν[target] == 0.0025", "0.0035" = "ν[target] == 0.0035"), label_parsed),
                                 eta =  as_labeller(c("0.0005" = "η == 0.0005", "0.005" = "η == 0.005"), label_parsed),
                                 sz =  as_labeller(c("0.1" = "s[z] == 0.1", "0.6" = "s[z] == 0.6"), label_parsed))) +
  theme(legend.position = "top") + labs(color = "RT bin") +
  xlab("Confidence bin") + ylab("Accuracy") + coord_cartesian(ylim = c(0.5, 1)) +
  theme(axis.title = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(-0.5, -0.5, -0.5, -0.5, "pt"))) -> p9_ddm
p9_ddm

# accuracy x rt x post confidence (binned for each condition)
simdat_ddm_bin <- c()
for (i in nu_tar) {
  for (j in eta) {
    for (k in sz) {
      subdat_ddm <- filter(simdat_ddm, nu_tar == i, eta == j, sz == k)
      subdat_ddm$rt_bin <- cut(subdat_ddm$rt, breaks = quantile(subdat_ddm$rt, probs = seq(0, 1, 0.25)), include.lowest = TRUE)
      levels(subdat_ddm$rt_bin) <- c("q1", "q2", "q3", "q4")
      subdat_ddm$conf_bin <- cut(subdat_ddm$post_conf, breaks = quantile(subdat_ddm$post_conf, probs = seq(0, 1, 0.25)), include.lowest = TRUE)
      levels(subdat_ddm$conf_bin) <- c("q1", "q2", "q3", "q4")
      simdat_ddm_bin <- rbind(simdat_ddm_bin, subdat_ddm)
    }
  }
}

simdat_ddm_bin %>%
  group_by(nu_tar, eta, sz, rt_bin, conf_bin) %>%
  summarise(accuracy = mean(correct)) %>%
  ggplot() + 
  geom_line(aes(x = conf_bin, y = accuracy, color = rt_bin, group = rt_bin), size = 0.3) +
  geom_point(aes(x = conf_bin, y = accuracy, color = rt_bin), size = 0.7) +
  scale_color_manual(values = c("#3300CC", "#7900FF", "#FF0077", "#CC0033")) +
  facet_wrap(~ nu_tar + eta + sz, ncol = 4,
             labeller = labeller(nu_tar = as_labeller(c("0.0015" = "ν[target] == 0.0015", "0.0025" = "ν[target] == 0.0025", "0.0035" = "ν[target] == 0.0035"), label_parsed),
                                 eta =  as_labeller(c("0.0005" = "η == 0.0005", "0.005" = "η == 0.005"), label_parsed),
                                 sz =  as_labeller(c("0.1" = "s[z] == 0.1", "0.6" = "s[z] == 0.6"), label_parsed))) +
  theme(legend.position = "top") + labs(color = "RT bin") +
  xlab("Confidence bin") + ylab("Accuracy") + coord_cartesian(ylim = c(0.5, 1)) +
  theme(axis.title = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(-0.5, -0.5, -0.5, -0.5, "pt")),
        plot.margin = unit(c(1, 1, 1, 1), "pt")) -> p10_ddm
p10_ddm

# RT only
simdat_ddm_bin %>%
  group_by(nu_tar, eta, sz, rt_bin) %>%
  summarise(accuracy = mean(correct)) %>%
  ggplot() + 
  geom_line(aes(x = rt_bin, y = accuracy, color = rt_bin, group = rt_bin), size = 0.3) +
  geom_point(aes(x = rt_bin, y = accuracy, color = rt_bin), size = 0.7) +
  scale_color_manual(values = c("#3300CC", "#7900FF", "#FF0077", "#CC0033")) +
  facet_wrap(~ nu_tar + eta + sz, ncol = 4,
             labeller = labeller(nu_tar = as_labeller(c("0.0015" = "ν[target] == 0.0015", "0.0025" = "ν[target] == 0.0025", "0.0035" = "ν[target] == 0.0035"), label_parsed),
                                 eta =  as_labeller(c("0.0005" = "η == 0.0005", "0.005" = "η == 0.005"), label_parsed),
                                 sz =  as_labeller(c("0.1" = "s[z] == 0.1", "0.6" = "s[z] == 0.6"), label_parsed))) +
  theme(legend.position = "top") + labs(color = "RT bin") +
  xlab("RT bin") + ylab("Accuracy") + coord_cartesian(ylim = c(0.5, 1)) +
  theme(axis.title = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(-0.5, -0.5, -0.5, -0.5, "pt")),
        plot.margin = unit(c(1, 1, 1, 1), "pt"))

# Confidence only
simdat_ddm_bin %>%
  group_by(nu_tar, eta, sz, conf_bin) %>%
  summarise(accuracy = mean(correct)) %>%
  ggplot() + 
  geom_line(aes(x = conf_bin, y = accuracy, group = conf_bin), size = 0.3) +
  geom_point(aes(x = conf_bin, y = accuracy), size = 0.7) +
  facet_wrap(~ nu_tar + eta + sz, ncol = 4,
             labeller = labeller(nu_tar = as_labeller(c("0.0015" = "ν[target] == 0.0015", "0.0025" = "ν[target] == 0.0025", "0.0035" = "ν[target] == 0.0035"), label_parsed),
                                 eta =  as_labeller(c("0.0005" = "η == 0.0005", "0.005" = "η == 0.005"), label_parsed),
                                 sz =  as_labeller(c("0.1" = "s[z] == 0.1", "0.6" = "s[z] == 0.6"), label_parsed))) +
  theme(legend.position = "top") + labs(color = "RT bin") +
  xlab("Confidence bin") + ylab("Accuracy") + coord_cartesian(ylim = c(0.5, 1)) +
  theme(axis.title = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(-0.5, -0.5, -0.5, -0.5, "pt")),
        plot.margin = unit(c(1, 1, 1, 1), "pt"))

# type-1 roc
roc_1_ddm %>% filter(variable != "be" & variable != "logit_be") %>%
  mutate(variable = fct_recode(variable, "confidence" = "post", "RT" = "RT", "logit" = "logit_post"),
         variable = fct_relevel(variable, "confidence", "RT", "logit")) %>%
  ggplot() + geom_point(aes(x = far_1, y = hr_1, color = factor(variable), size = factor(variable)), alpha = 0.7) +
  guides(size = "none", color = guide_legend(override.aes = list(size = 2))) +
  scale_size_manual(values = c("confidence" = 1.8, "RT" = 0.8, "logit" = 0.8)) +  
  geom_line(roc_1_sdt_ddm, mapping = aes(x = far_1, y = hr_1), linetype = "dashed") +
  scale_color_manual(values = c(hue_pal()(3)), guide = guide_legend(override.aes = list(color = "white"))) +
  facet_wrap(. ~ nu_tar + eta + sz, ncol = 4,
             labeller = labeller(nu_tar = as_labeller(c("0.0015" = "ν[target] == 0.0015", "0.0025" = "ν[target] == 0.0025", 
                                                           "0.0035" = "ν[target] == 0.0035"), label_parsed),
                                 eta = as_labeller(c("0.0005" = "η == 0.0005", "0.005" = "η == 0.005"), label_parsed),
                                 sz = as_labeller(c("0.1" = "s[z] == 0.1", "0.6" = "s[z] == 0.6"), label_parsed))) +
  xlab("p(Resp = S2|Stim = S1)") + ylab("p(Resp = S2|Stim = S2)") + labs(color = "Variable") + 
  scale_x_continuous(breaks = c(0, 0.5, 1)) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme(axis.title = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(-0.5, -0.5, -0.5, -0.5, "pt"))) -> p11_ddm
p11_ddm
p11_ddm <- p11_ddm + theme(axis.text.x = element_text(angle = 25),
                           strip.text = element_text(size = 6.5, face = "bold"),
                           plot.margin = margin(5, 5, 5, 5), 
                           plot.title = element_text(size = 8))

# type-2 roc
roc_2_ddm %>% filter(variable != "be" & variable != "logit_be") %>%
  mutate(variable = fct_recode(variable, "confidence" = "post", "RT" = "RT", "logit" = "logit_post"),
         variable = fct_relevel(variable, "confidence", "RT", "logit")) %>%
  ggplot() + geom_point(aes(x = far_2, y = hr_2, color = factor(variable), size = factor(variable)), alpha = 0.7) +
  guides(size = "none", color = guide_legend(override.aes = list(size = 2))) +
  scale_size_manual(values = c("confidence" = 1.8, "RT" = 0.8, "logit" = 0.8)) +  
  geom_line(roc_2_sdt_ddm, mapping = aes(x = far_2, y = hr_2), linetype = "dashed") +
  scale_color_manual(values = c(hue_pal()(3)), guide = guide_legend(override.aes = list(color = "white"))) +
  facet_wrap(. ~ nu_tar + eta + sz, ncol = 4,
             labeller = labeller(nu_tar = as_labeller(c("0.0015" = "ν[target] == 0.0015", "0.0025" = "ν[target] == 0.0025", 
                                                           "0.0035" = "ν[target] == 0.0035"), label_parsed),
                                 eta = as_labeller(c("0.0005" = "η == 0.0005", "0.005" = "η == 0.005"), label_parsed),
                                 sz = as_labeller(c("0.1" = "s[z] == 0.1", "0.6" = "s[z] == 0.6"), label_parsed))) +
  xlab("p(SV = high|Resp = incorrect)") + ylab("p(SV = high|Resp = correct)") + labs(color = "Variable") + 
  scale_x_continuous(breaks = c(0, 0.5, 1)) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme(axis.title = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(-0.5, -0.5, -0.5, -0.5, "pt"))) -> p12_ddm
p12_ddm
p12_ddm <- p12_ddm + theme(axis.text.x = element_text(angle = 25),
                           strip.text = element_text(size = 6.5, face = "bold"),
                           plot.margin = margin(5, 5, 5, 5), 
                           plot.title = element_text(size = 8))

# mean early confidence
simdat_ddm %>%
  group_by(choice, nu_tar, eta, sz) %>%
  mutate(mean_conf = mean(early_conf)) %>%
  ggplot() + geom_line(aes(x = nu_tar, y = mean_conf, color = factor(choice))) +
  geom_point(aes(x = nu_tar, y = mean_conf, color = factor(choice))) +
  scale_color_discrete(labels = c("target", "distractor")) + labs(color = "choice") +
  facet_wrap(. ~ eta + sz, labeller = labeller(eta = ~ paste("σ_dr:", .),
                                                          sz = ~ paste("α_sp:", .),
                                                          .multi_line = FALSE)) +
  ylab("Mean early confidence") + scale_color_npg() -> p13_ddm
p13_ddm <- p13_ddm + guides(color = F) + theme(axis.text.x  = element_text(angle = 20))
p13_ddm

# interaction
simdat_ddm$rt_bin <- cut(simdat_ddm$rt, breaks = quantile(simdat_ddm$rt, probs = seq(0, 1, 0.05)), labels = 1:20, include.lowest = TRUE)
simdat_ddm$rt_bin <- as.numeric(simdat_ddm$rt_bin)
simdat_ddm$conf_bin <- cut(simdat_ddm$post_conf, breaks = quantile(simdat_ddm$post_conf, probs = seq(0, 1, 1/3)), include.lowest = TRUE)
levels(simdat_ddm$conf_bin) <- c("q1", "q2", "q3")
simdat_ddm %>%
 # group_by(nu_tar, eta, sz, rt_bin, conf_bin) %>%
 #  summarise(accuracy = mean(correct)) %>%
  ggplot(aes(x = rt_bin, y = correct, color = conf_bin)) + 
#  geom_line(aes(x = rt_bin, y = accuracy, color = conf_bin, group = conf_bin), size = 0.3) +
#  geom_point(aes(x = rt_bin, y = accuracy, color = conf_bin), size = 0.7) +
  geom_smooth(method = "loess", se = TRUE, size = 0.5) +
  scale_color_manual(values = c("#3300CC", "#7900FF", "#FF0077")) +
  facet_wrap(~ nu_tar + eta + sz, ncol = 4,
             labeller = labeller(nu_tar = as_labeller(c("0.0015" = "ν[target] == 0.0015", "0.0025" = "ν[target] == 0.0025", "0.0035" = "ν[target] == 0.0035"), label_parsed),
                                 eta =  as_labeller(c("0.0005" = "η == 0.0005", "0.005" = "η == 0.005"), label_parsed),
                                 sz =  as_labeller(c("0.1" = "s[z] == 0.1", "0.6" = "s[z] == 0.6"), label_parsed))) +
  labs(color = "Confidence bin") +
  xlab("RT bin") + ylab("Accuracy") + coord_cartesian(ylim = c(0.5, 1)) +
  theme(axis.title = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(-0.5, -0.5, -0.5, -0.5, "pt"))) -> p14_ddm


# write.csv(simdat_ddm, "simdat_ddm.csv", row.names =  F)

save_plot("figure_6a.jpg", p1_ddm, width = 5, height = 4.1)
save_plot("figure_6b.jpg", p3_ddm, width = 5, height = 4.1)
save_plot("figure_6c.jpg", p4_ddm, width = 5, height = 4.1)

save_plot("figure_7a.jpg", p5_ddm, width = 5.5, height = 4.1)
save_plot("figure_7b.jpg", p6_ddm, width = 5.5, height = 4.1)
save_plot("figure_7c.jpg", p7_ddm, width = 5.5, height = 4.1)

save_plot("figure_8.jpg", p14_ddm, width = 12, height = 12.5)
save_plot("figure_a4.png", p11_ddm, width = 12, height = 12.5)
save_plot("figure_a5.png", p12_ddm, width = 12, height = 12.5)
save_plot("figure_a6.png", p2_ddm)