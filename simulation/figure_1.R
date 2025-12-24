library(tidyverse)
library(data.table)
library(grid)
library(cowplot)
library(sjPlot)
library(metaSDT)
library(GGally)
library(ppcor)

#'#
dat <- fread("data_Hainguerlot_2018.csv", header = T)
i <- 4
d <- subset(dat, dat$Subj_idx == i)
d <- na.omit(d)
d <- mutate(d, Correct = ifelse(Stimulus == Response, 1, 0))

q <- quantile(d$RT_dec, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
d$RT_bin <- cut(d$RT_dec, breaks = q, labels = FALSE, include.lowest = TRUE)

nR_S1 <- c(sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 1, na.rm = TRUE),
           sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 2, na.rm = TRUE),
           sum(d$Stimulus == "1" & d$Response == "1" & d$RT_bin == 3, na.rm = TRUE),
           sum(d$Stimulus == "1" & d$Response == "2" & d$RT_bin == 3, na.rm = TRUE),
           sum(d$Stimulus == "1" & d$Response == "2" & d$RT_bin == 2, na.rm = TRUE),
           sum(d$Stimulus == "1" & d$Response == "2" & d$RT_bin == 1, na.rm = TRUE))

nR_S2 <- c(sum(d$Stimulus == "2" & d$Response == "1" & d$RT_bin == 1, na.rm = TRUE),
           sum(d$Stimulus == "2" & d$Response == "1" & d$RT_bin == 2, na.rm = TRUE),
           sum(d$Stimulus == "2" & d$Response == "1" & d$RT_bin == 3, na.rm = TRUE),
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

fit_r[1, 1:5]
auc1_r
auc2_r
  
#'#
d <- mutate(d, 
            Stimulus = ifelse(Stimulus == 1, "S1", "S2"),
            Response = ifelse(Response == 1, "S1", "S2"),
            Correct  = ifelse(Correct  == 1, "correct", "incorrect"))

d <- d %>%
  mutate(grp1 = paste0(Response, ", B", RT_bin))
desired_order <- c("S1, B1", "S1, B2","S1, B3",
                   "S2, B3", "S2, B2","S2, B1")
d$grp1 <- factor(d$grp1, levels = desired_order)

ggplot(d, aes(x = grp1, fill = Stimulus)) +
  geom_bar(position = "dodge", alpha = 0.7) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_x_discrete(
    labels = c(
      "Resp = S1\nRT = B1",
      "Resp = S1\nRT = B2",
      "Resp = S1\nRT = B3",
      "Resp = S2\nRT = B3",
      "Resp = S2\nRT = B2",
      "Resp = S2\nRT = B1"
    )
  ) +
  ylab("Frequency") +
  xlab(NULL) +
  labs(fill = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 6,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 6,
      color = "black"
    ),
    legend.position = "none"
  ) -> g1

g1


#'#
cuts <- quantile(d$RT_dec, probs = seq(0, 1, length.out = 4))
cuts <- cuts[2:3]
cuts

labels_1 <- lapply(1:2, function(i) bquote(t[.(i)]))
labels_2 <- c("B1", "B2", "B3") 
x_pos <- c(0.66, 0.88, 1.1)

ggplot(d) +
  geom_histogram(
    aes(x = RT_dec, fill = Stimulus),
    binwidth = diff(range(d$RT_dec)) / 30,
    alpha = 0.7
  ) +
  geom_vline(xintercept = cuts, linetype = "dashed", color = "black", size = 0.55) +
  annotate(
    "text",
    x = cuts,
    y = Inf,
    label = labels_1,
    hjust = 1.3,
    vjust = 1.1,
    parse = TRUE) +
  annotate(
    "text",
    x = x_pos,
    y = 55,
    label = labels_2) +
  scale_fill_manual(values = c("red", "blue"),
                    name = "Stim") +
  theme_minimal(base_size = 8) +
  xlab("RT") +
  ylab("Frequency") +
  theme(
    axis.text.x = element_text(size = 6, color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = c(0.805, 0.875)) + 
  ylim(0, 60) -> g2

g2

#'#
d <- mutate(d, bin = ifelse(RT_bin == 1, "RT = B1",
                            ifelse(RT_bin == 2, "RT = B2", "RT = B3")))
ggplot(d, aes(x = bin, fill = factor(Correct))) +
  geom_bar(position = "dodge", alpha = 1.0) +
  scale_fill_manual(values = c("#1BCAD3", "#3A0F6B"),
                    name = "Resp" ) +
  ylab("Frequency") +
  xlab(NULL) +
  labs(fill = NULL) +
  theme_minimal(base_size = 8) +
  ylim(0, 150) +
  theme(
    legend.position = c(0.84, 0.88),
    axis.text.x = element_text(angle = 0, size = 6, color = "black"),
    axis.text.y = element_text(size = 6, color = "black")) -> g3
g3

c("#1F3A5F", "#E36414")
c("#264653", "#E9C46A")
c("#0B1D2A", "#B08968")
c("#3D405B", "#4EA8DE")
c("#005F73", "#EE9B00")
c("#1B4332", "#F1C27D")
c("#1BCAD3", "#3A0F6B")


#'#
roc2 <- data.frame(far2 = c(0, p_fa2), hr2 = c(0, p_hit2))
  
ggplot(roc2[2:3, ], aes(x = far2, y = hr2)) + 
  geom_point(size = 2) +
  geom_line(roc2, mapping = aes(x = far2, y = hr2)) +
  xlab("p(RT = fast|Resp = incorrect)") + 
  ylab("p(RT = fast|Resp = correct)") +
  coord_equal() +
  xlim(0, 1) + ylim(0, 1) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 0, size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black")) -> g4
g4

g <- cowplot::plot_grid(g2, g1, g3, g4, labels = c("(a)", "(b)", "(c)", "(d)"), 
                        nrow = 2, label_x = -0.02, label_y = 1.01, label_size = 10)
ggsave("figure_1.jpg", g, width = 6, height = 6,    units = "in", dpi = 500)