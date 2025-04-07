#+ message = F, warning = F
library(tidyverse)
library(gghalves)
library(GGally)
library(car)
library(doParallel)
library(sjPlot)
library(cowplot)
library(gridExtra)

theme_set(theme_publication()) 

script <- c("Maniscalco_2017_expt1.R",
            "Maniscalco_2017_expt3.R",
            "Maniscalco_2017_expt4.R",
            "Hainguerlot_2018.R", 
            "Hainguerlot_unpub.R",
            "Reyes_2015.R")

exp <- c("Maniscalco_2017_expt1",
         "Maniscalco_2017_expt3",
         "Maniscalco_2017_expt4",
         "Hainguerlot_2018", 
         "Hainguerlot_unpub",
         "Reyes_2015")

df_1 <- c()
df_2 <- c()
figures <- c()

for (j in 1:length(script)) {
  source(script[j]) 
  df1$exp <- exp[j]
  df2$exp <- exp[j]
  
  df_1 <- dplyr::bind_rows(df_1, df1) 
  df_2 <- dplyr::bind_rows(df_2, df2)
}  

# visualization
df_2$RT_bin <- as.numeric(df_2$RT_bin)
df_2$exp <- factor(df_2$exp, levels = exp)

ggplot(df_2, aes(x = RT_bin, y = Correct, color = as.factor(Confidence)), group = RT_bin) +
# geom_jitter(height = 0.01, alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE, size = 0.5) +
  scale_color_manual(values = c("#3300CC", "#7900FF", "#FF0077", "#CC0033")) +
  theme(legend.position = "top") + labs(color = "Confidence bin") +
  labs(x = "RT bin", y = "Accuracy", color = "Confidence") +
  facet_wrap(.~ exp) +
  theme(axis.title = element_text(size = 7), legend.position = "top", plot.title = element_text(size = 7),
        strip.text = element_text(size = 7, margin = margin(0.5, 0.5, 0.5, 0.5, "pt")),
        plot.margin = unit(c(1, 1, 1, 1), "pt")) +  coord_cartesian(ylim = c(0.5, 1)) -> p1
p1
sjPlot::save_plot("figure_3.jpg", p1, width = 12, height = 8, dpi = 300)


  guides(color = F) 
p
save_plot("sliding_window2.jpg", p)
ggsave("sliding_window_2.jpg", p)