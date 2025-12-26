sdt_ddm %>%
  pivot_longer(cols = c("mratio_post", "mratio_rt", "mratio_logit_post"), 
               names_to = "index", values_to = "value") %>%
  mutate(Variable = fct_recode(index, "confidence" = "mratio_post", "RT" = "mratio_rt", "confidence+RT" = "mratio_logit_post"),
         Variable = fct_relevel(Variable, "confidence", "RT", "confidence+RT")) %>%
  ggplot() + geom_point(aes(x = nu_tar, y = value, color = Variable, size = Variable)) + 
  scale_size_manual(values = c("confidence" = 0.8, "RT" = 0.25, "confidence+RT" = 0.25)) + labs(color = "Variable") +
  geom_line(aes(x = nu_tar, y = value, color = Variable), size = 0.25) +
  facet_wrap(~ eta + sz, ncol = 2,
             labeller = labeller( eta =  as_labeller(c("0.0005" = "η == 0.0005", "0.005" = "η == 0.005"), label_parsed),
                                  sz =  as_labeller(c("0.1" = "s[z] == 0.1", "0.6" = "s[z] == 0.6"), label_parsed))) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(legend.text = element_text(size = 7),
        legend.position = "top", 
        legend.direction = "horizontal",
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = NA),
        legend.title = element_text(size = 8),
        strip.text = element_text(size = 7, margin = margin(-0.2, -0.2, -0.2, -0.2, "pt"))) +
  scale_color_manual(values = c(hue_pal()(3))) + 
  scale_x_continuous(breaks = c(0.0015, 0.0025, 0.0035)) + xlab(expression(bold(paste({ν[target]})))) + ylab("m-ratio") +
  coord_cartesian(ylim = c(-0.4, 1.1)) + scale_y_continuous(breaks = seq(0, 1, 0.5)) -> f

save_plot("figure.jpg", f, width = 10*0.76, height = 8.2*0.76)