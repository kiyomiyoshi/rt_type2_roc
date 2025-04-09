library(magrittr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggsci)

theme_set(theme_publication()) 


#'# simulation
mu_dr_tar <- 0.0009
mu_dr_dis <- 0.0006
sigma_dr <-  0.00002
sigma_ev <-  0.015
alpha_sp <-  0.1
tau_nd <-    0.3
theta_vrm <- 2.0
theta_ddm <- 1.2
sample <- 3201

dr <- cbind(rnorm(1, mu_dr_tar, sigma_dr),
            rnorm(1, mu_dr_dis, sigma_dr))
a_1 <- runif(1, -alpha_sp, alpha_sp)
a_2 <- runif(1, -alpha_sp, alpha_sp)

accum <- matrix(0, nrow = sample + 1, ncol = 2)
accum[, 1] <- c(a_1, rep(dr[1, 1], sample))
accum[, 2] <- c(a_2, rep(dr[1, 2], sample))
accum[-1, 1] <- accum[-1, 1] + rnorm(sample, 0, sigma_ev)
accum[-1, 2] <- accum[-1, 2] + rnorm(sample, 0, sigma_ev)
accum <- matrixStats::colCumsums(accum)
be <- matrixStats::rowOrderStats(-accum, which = 2) - matrixStats::rowOrderStats(-accum, which = 1)
accum <- as.data.frame(accum)
accum$be <- be
accum <- mutate(accum, t = row_number())
colnames(accum) <- c("E1", "E2", "BE", "t")

rt_vrm <- which(pmax(accum$E1, accum$E2) > theta_vrm)[1]
rt_vrm
rt_ddm <- which(accum$BE > theta_ddm)[1]
rt_ddm
accum[rt_vrm, ]
accum[rt_ddm, ]


#'# visualization
accum[(rt_vrm + 1):nrow(accum), c("E1", "E2")] <- NA
accum[(rt_ddm + 150):nrow(accum), "BE"] <- NA

g1 <- ggplot(accum[1:2000, ]) + 
  geom_line(aes(x = t, y = E1), color = "red", size = 0.3) + ylab("Evidence for S1") +
  geom_hline(yintercept = 2, linetype = "dashed", size = 0.3) +
  scale_x_continuous(breaks = c(0, 1000, 2000)) + ylim(-0.5, 2.2) + 
  theme(axis.title.x = element_text(size = 5, margin = margin(t = -2, unit = "pt")),
        axis.title.y = element_text(size = 5),
        axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(3.1, 3.1, 3.1, 3.1, unit = "pt"))

g2 <- ggplot(accum[1:2000, ]) + 
  geom_line(aes(x = t, y = E2), color = "blue", size = 0.3) + ylab("Evidence for S2") +
  geom_hline(yintercept = 2, linetype = "dashed", size = 0.3) +
  scale_x_continuous(breaks = c(0, 1000, 2000)) + ylim(-0.5, 2.2) + labs(color = NULL) +
  theme(axis.title.x = element_text(size = 5, margin = margin(t = -2, unit = "pt")),
        axis.title.y = element_text(size = 5),
        axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(3.1, 3.1, 3.1, 3.1, unit = "pt"))

g3 <- ggplot(accum[1:2000, ]) + 
  geom_line(aes(x = t, y = BE), color = "purple", size = 0.3) + ylab("Balance of evidence") +
  geom_hline(yintercept =  1.2, linetype = "dashed", size = 0.3) +
  geom_hline(yintercept = -1.2, linetype = "dashed", size = 0.3) +
  annotate("rect", xmin = rt_ddm, xmax = rt_ddm + 150, ymin = -1.55, ymax = 1.55, alpha = 0.2, fill = "grey") +
  scale_x_continuous(breaks = c(0, 1000, 2000)) + 
  scale_y_continuous(breaks = c(-1.2, 0, 1.2), lim = c(-1.55, 1.55)) + labs(color = NULL) +
  theme(axis.title.x = element_text(size = 5, margin = margin(t = -2, unit = "pt")),
        axis.title.y = element_text(size = 5),
        axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(3.1, 3.1, 3.1, 3.1, unit = "pt"))
g3

figure_1 <- plot_grid(g1, g2, g3, ncol = 3, scale = 1)
figure_1
sjPlot::save_plot("Figure_1.jpg", figure_1, width = 9, height = 2.5, dpi = 300)
# write.csv(accum, "figure_1.csv")