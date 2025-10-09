#+ message = F, warning = F
library(tidyverse)
library(data.table)
library(gghalves)
library(GGally)
library(car)
library(doParallel)
library(sjPlot)
library(ggdist)
library(metafor)

data_filtered <- fread("data_filtered.csv", header = T)
data_sh_filtered <- fread("data_sh_filtered.csv", header = T)

cor_by_exp <- data_filtered %>%
  group_by(exp) %>%
  summarise(
    n = sum(complete.cases(mdp_conf, mdp_rt)),
    r = cor(mdp_conf, mdp_rt, use = "pairwise.complete.obs") 
  )
cor_by_exp

wide <- data_sh_filtered %>%
  pivot_wider(
    id_cols = c(exp, sub),
    names_from = rn,
    values_from = c(mdp_conf, mdp_rt)
  )
nrow(wide)
summary(wide)

reliability <- wide %>%
  group_by(exp) %>%
  summarise(
    r_conf = cor(mdp_conf_odd, mdp_conf_even, use = "pairwise.complete.obs"),
    r_rt   = cor(mdp_rt_odd, mdp_rt_even, use = "pairwise.complete.obs"),
    r_conf_sb = 2 * r_conf / (1 + r_conf),  # Spearman-Brown correction
    r_rt_sb   = 2 * r_rt / (1 + r_rt)
  )
reliability

# attenuation correction
df_joined <- cor_by_exp %>%
  left_join(reliability, by = "exp") %>%
  select(exp, n, r, r_conf, r_rt, r_conf_sb, r_rt_sb)

df_joined <- df_joined %>%
  mutate(
    r_corrected = r / sqrt(r_conf_sb * r_rt_sb)  )
df_joined

# screening & clipping
df_joined <- df_joined %>%
  mutate(
    # screening
    r_corrected = if_else(
      r_conf_sb > 0.3 & r_rt_sb > 0.3,
      r / sqrt(r_conf_sb * r_rt_sb),
      NA_real_
    ),
    # clipping
    r_corrected = case_when(
      is.na(r_corrected) ~ NA_real_,
      r_corrected >  1 ~ 1,
      r_corrected < -1 ~ -1,
      TRUE ~ r_corrected
    )
  )
df_joined

df_complete <- df_joined %>% 
  filter(!is.na(r_corrected))
df_complete

# meta-analysis
r_to_z <- function(r) 0.5 * log((1 + r)/(1 - r))

# Fisher's z transformation
df_meta <- df_complete %>%
  mutate(
    z_r = r_to_z(r),
    z_r_corrected = r_to_z(r_corrected),
    se_z_r = 1 / sqrt(n - 3),
    se_z_r_corrected = 1 / sqrt(n - 3)
  ) %>%
  select(exp, n, r, r_corrected, z_r, z_r_corrected, se_z_r, se_z_r_corrected)

clip_r <- function(r, eps = 1e-6) {
  pmax(pmin(r, 1 - eps), -1 + eps)
}

eps <- 1e-6
df_meta <- df_meta %>%
  mutate(
    r_corrected_clipped = pmin(pmax(r_corrected, -1 + eps), 1 - eps),
    z_r_corrected_clipped = r_to_z(r_corrected_clipped)
  )
df_meta

df_meta <- df_meta %>%
  mutate(
    r_clipped = pmin(pmax(r, -1 + 1e-6), 1 - 1e-6),
    z_r_clipped = r_to_z(r_clipped)
  )

# zero-order correlation
yi_r <- df_meta$z_r_clipped
sei_r <- df_meta$se_z_r
res_r <- rma(yi = yi_r, sei = sei_r, method = "REML")

# disattenuated correlation
yi_rc <- df_meta$z_r_corrected_clipped
sei_rc <- df_meta$se_z_r_corrected
res_rc <- rma(yi = yi_rc, sei = sei_rc, method = "REML")

z_to_r <- function(z) {
  (exp(2*z) - 1) / (exp(2*z) + 1)
}

list(
  original_r = list(
    r = z_to_r(res_r$b),
    ci_lb = z_to_r(res_r$ci.lb),
    ci_ub = z_to_r(res_r$ci.ub)
  ),
  corrected_r = list(
    r = z_to_r(res_rc$b),
    ci_lb = z_to_r(res_rc$ci.lb),
    ci_ub = z_to_r(res_rc$ci.ub)
  )
)