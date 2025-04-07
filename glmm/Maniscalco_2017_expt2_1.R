library(lme4)
library(lmerTest)
library(tidyverse)
library(data.table)
library(ggeffects)
library(car)
library(metaSDT)

dat <- fread("data_Maniscalco_2017_expt2.csv", header = T)
dat <- subset(dat, dat$ContrastLevel == "1")
dat <- mutate(dat, Correct = ifelse(Stimulus == Response, 1, 0))
dat <- na.omit(dat)

m1 <- glmer(Correct ~ Confidence + (1 + Confidence | Subj_idx), family = binomial, data = dat) 
summary(m1)
Anova(m1)

m2 <- glmer(Correct ~ Confidence + RT_dec + (1 + Confidence + RT_dec | Subj_idx), family = binomial, data = dat) 
summary(m2)
Anova(m2)

m3 <- glmer(Correct ~ Confidence * RT_dec + (1 + Confidence + RT_dec | Subj_idx), family = binomial, data = dat) 
summary(m3)
Anova(m3)

dat <- mutate(dat, logit1 = predict(m1),
              logit2 = predict(m2),
              logit3 = predict(m3))

q1 <- quantile(dat$logit1, probs = seq(0, 1, by = 0.1), na.rm = TRUE)
q2 <- quantile(dat$logit2, probs = seq(0, 1, by = 0.1), na.rm = TRUE)
q3 <- quantile(dat$logit3, probs = seq(0, 1, by = 0.1), na.rm = TRUE)
all_levels <- seq(length(q1) - 1)

dat$logit_bin1 <- factor(cut(dat$logit1, breaks = q1, labels = FALSE, include.lowest = TRUE), levels = all_levels)
dat$logit_bin2 <- factor(cut(dat$logit2, breaks = q2, labels = FALSE, include.lowest = TRUE), levels = all_levels)
dat$logit_bin3 <- factor(cut(dat$logit3, breaks = q3, labels = FALSE, include.lowest = TRUE), levels = all_levels)

roc_1 <- c()
roc_2 <- c()
sdt <- c()
for (i in 0:2) {
  
  column_index <- 14 + i
  dat$logit_bin = dat[, ..column_index]
  
  h1 <- c(rev(table(dat[which(dat$Response == 0 & dat$Stimulus == 0), "logit_bin"])),
          table(dat[which(dat$Response == 1 & dat$Stimulus == 0), "logit_bin"]))
  f1 <- c(rev(table(dat[which(dat$Response == 0 & dat$Stimulus == 1), "logit_bin"])),
          table(dat[which(dat$Response == 1 & dat$Stimulus == 1), "logit_bin"]))
  df1 <- data.frame(hr_1 = cumsum(h1)/sum(h1), far_1 = cumsum(f1)/sum(f1))
  df1$variable <- i + 1
  df1$order <- seq_len(nrow(df1))
  
  h2 <- rev(table(dat[which(dat$Correct == 1), "logit_bin"]))
  f2 <- rev(table(dat[which(dat$Correct == 0), "logit_bin"]))
  df2 <- data.frame(hr_2 = cumsum(h2)/sum(h2), far_2 = cumsum(f2)/sum(f2))
  df2$variable <- i + 1
  df2$order <- seq_len(nrow(df2))
  
  fit <-  fit_meta_d_MLE(h1, f1)
  fit$variable <- i + 1
  
  roc_1 <- rbind(roc_1, df1)
  roc_2 <- rbind(roc_2, df2)
  sdt <- rbind(sdt, fit)
  
}

# sdt prediction
sdt <- distinct(sdt[, 1:5])[, -2]
sdt <- mutate(sdt, variable = factor(c("conf_logit", "additive_logit", "interaction_logit"), 
                                     levels = c("conf_logit", "additive_logit", "interaction_logit")))
dp <- sdt[1, 1]
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

h2 <- resp2[1 : (length(resp2) / 2)] + rev(resp1[(length(resp1) / 2 + 1) : length(resp1)])
f2 <- rev(resp2[(length(resp2) / 2 + 1) : length(resp2)]) + resp1[1: (length(resp1) / 2)]

roc_1_sdt <- data.frame(hr_1 = cumsum(resp2)/sum(resp2), far_1 = cumsum(resp1)/sum(resp1))
roc_2_sdt <- data.frame(hr_2 = cumsum(h2)/sum(h2), far_2 = cumsum(f2)/sum(f2))

# visualization
g1 <- ggplot(roc_1) + geom_point(aes(x = far_1, y = hr_1, color = factor(variable))) +
  geom_line(roc_1_sdt, mapping = aes(x = far_1, y = hr_1), linetype = "dashed")
g1

g2 <- ggplot(roc_2) + geom_point(aes(x = far_2, y = hr_2, color = factor(variable))) +
  geom_line(roc_2_sdt, mapping = aes(x = far_2, y = hr_2), linetype = "dashed")
g2

sdt %>%
  pivot_longer(cols = "meta_da") %>%
  ggplot() + geom_col(aes(x = variable, y = value, fill = variable)) +
  guides(fill = F) -> g3
g3


#'#
m4 <- lmer(Confidence ~ RT_dec * Correct + (1 + RT_dec + Correct | Subj_idx), data = dat) 
summary(m4)
Anova(m4)
plot(ggpredict(m4, terms = c("RT_dec", "Correct")))