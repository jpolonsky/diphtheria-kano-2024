# Home-based care vs facility-based care for mild diphtheria, Kano State,
# Nigeria, 2023-24: retrospective matched cohort study.
#
# Reproduces every number, table and figure reported in the manuscript.
# Input is the frozen, de-identified matched cohort sample; see README for
# how to request it.

pacman::p_load(
  survival,
  survey,
  srvyr,
  lme4,
  lmerTest,
  broom,
  broom.mixed,
  gt,
  rms,
  sandwich,
  lmtest,
  EValue,
  scales,
  tidyverse
)

if (!requireNamespace("logistf", quietly = TRUE)) {
  install.packages("logistf")
}
if (!requireNamespace("ordinal", quietly = TRUE)) {
  install.packages("ordinal")
}
if (!requireNamespace("car", quietly = TRUE)) {
  install.packages("car")
}

fig_dir <- "output/figures"
tab_dir <- "output/tables"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 0. Data
# ============================================================
# df_pt: one row per interviewed patient (n=669)
# df_cohort_samp: the originally-sampled cohort before tracing/replacement
#   (n=986), used for the attrition comparison and selection weights below
# df_hh: one row per household contact (n=5043), used for the household
#   transmission analysis

df_pt <- read_rds("data/df_pt.rds")
df_cohort_samp <- read_rds("data/df_cohort_samp.rds")
df_hh <- read_rds("data/df_hh.rds")

# DTC is the reference level throughout this script; HBC is always the
# comparator, matching how effects are reported in the manuscript text
# (e.g. "HBC vs DTC: aOR 0.36").
df_pt <- df_pt |> mutate(tx_mod = fct_relevel(tx_mod, "DTC"))
df_hh <- df_hh |> mutate(tx_mod = fct_relevel(tx_mod, "DTC"))

options(survey.lonely.psu = "adjust")
srv_design <- df_pt |> as_survey_design(strata = subclass, ids = NULL)

# ============================================================
# 1. Participants (Table 1, Figure 1 data, attrition)
# ============================================================

nrow(df_pt)
df_pt |> count(tx_mod)
df_pt |> count(lga)

# matched subclasses vary in size, and those holding one arm only don't
# contribute to the conditional estimates
subclass_sizes <- df_pt |>
  group_by(subclass) |>
  summarise(size = n(), n_arms = n_distinct(tx_mod))
nrow(subclass_sizes)
range(subclass_sizes$size)
subclass_sizes |>
  filter(n_arms == 1) |>
  summarise(n_subclasses = n(), n_patients = sum(size))

df_pt |>
  summarise(
    median_age = median(age, na.rm = TRUE),
    q1 = quantile(age, .25, na.rm = TRUE),
    q3 = quantile(age, .75, na.rm = TRUE)
  )

srv_design |>
  summarise(
    n = n(),
    pct_vaccinated = survey_mean(
      (vaccinated == 1) * 100,
      na.rm = TRUE,
      vartype = "ci"
    )
  )

srv_design |>
  group_by(tx_mod) |>
  summarise(
    pct_vaccinated = survey_mean(
      (vaccinated == 1) * 100,
      na.rm = TRUE,
      vartype = "ci"
    )
  )

df_pt |>
  group_by(tx_mod) |>
  summarise(
    comorbidity_n = sum(comorbidity, na.rm = TRUE),
    complication_n = sum(complication, na.rm = TRUE)
  )

df_pt |>
  group_by(tx_mod) |>
  summarise(
    median_hh = median(hh_size, na.rm = TRUE),
    q1 = quantile(hh_size, .25, na.rm = TRUE),
    q3 = quantile(hh_size, .75, na.rm = TRUE)
  )
wilcox.test(hh_size ~ tx_mod, data = df_pt)

df_pt |>
  group_by(tx_mod) |>
  summarise(
    median_delay = median(delay_tx, na.rm = TRUE),
    q1 = quantile(delay_tx, .25, na.rm = TRUE),
    q3 = quantile(delay_tx, .75, na.rm = TRUE)
  )
wilcox.test(delay_tx ~ tx_mod, data = df_pt)

# Table 1 (gt, Word)
get_distribution <- function(df, var) {
  var_name <- deparse(substitute(var))
  df |>
    count(tx_mod, {{ var }}) |>
    group_by(tx_mod) |>
    mutate(
      Characteristic = var_name,
      value = as.character({{ var }}),
      prop = n / sum(n) * 100
    )
}

table_1 <- bind_rows(
  get_distribution(df_pt, lga),
  get_distribution(df_pt, sex),
  get_distribution(df_pt, age_grp),
  get_distribution(df_pt, vaccinated),
  get_distribution(df_pt, comorbidity),
  get_distribution(df_pt, complication),
  get_distribution(df_pt, outcome_cat)
) |>
  mutate(CountPct = sprintf("%d (%.1f%%)", n, prop)) |>
  select(Characteristic, value, tx_mod, CountPct) |>
  pivot_wider(names_from = tx_mod, values_from = CountPct) |>
  mutate(across(c(HBC, DTC), ~ replace_na(.x, "0 (0.0%)")))

table_1 |>
  gt() |>
  tab_header(
    title = md(
      "**Table 1. Characteristics of enrolled participants by treatment modality**"
    )
  ) |>
  cols_label(
    Characteristic = md("**Characteristic**"),
    value = md("**Category**")
  ) |>
  tab_options(table.font.name = "Times New Roman", table.font.size = 11) |>
  gtsave(file.path(tab_dir, "table1.docx"))

# Attrition: how the enrolled cohort relates to the originally-sampled one
df_attr <- df_cohort_samp |> mutate(traced = case_id %in% df_pt$case_id)

nrow(df_attr)
sum(df_attr$traced)
mean(df_attr$traced)

is_replacement <- !(df_pt$case_id %in% df_cohort_samp$case_id)
sum(is_replacement)
df_pt |>
  mutate(is_replacement = is_replacement) |>
  group_by(tx_mod, is_replacement) |>
  summarise(n = n(), deaths = sum(died), .groups = "drop")

for (v in c("tx_location", "sex", "age_group", "vaccination_status")) {
  tab <- table(df_attr[[v]], df_attr$traced)
  print(tab)
  p <- chisq.test(tab)$p.value
  cat(sprintf("%s: p=%.4f\n", v, p))
}

# S5 Table: standardised mean differences on the matching variables (should
# be ~0 by construction) and on two variables not used for matching, where
# imbalance is expected
smd_binary <- function(x, tx) {
  p1 <- mean(x[tx == "HBC"], na.rm = TRUE)
  p2 <- mean(x[tx == "DTC"], na.rm = TRUE)
  s <- sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
  if (s == 0) {
    return(0)
  }
  (p1 - p2) / s
}
smd_continuous <- function(x, tx) {
  m1 <- mean(x[tx == "HBC"], na.rm = TRUE)
  m2 <- mean(x[tx == "DTC"], na.rm = TRUE)
  v1 <- var(x[tx == "HBC"], na.rm = TRUE)
  v2 <- var(x[tx == "DTC"], na.rm = TRUE)
  (m1 - m2) / sqrt((v1 + v2) / 2)
}

matched_vars <- list(
  sex = "sex_reg",
  age_group = "age_grp_reg",
  vaccination_status = "vacc_status_reg",
  lga = "lga"
)
balance_matched <- map_dfr(names(matched_vars), function(vname) {
  col <- matched_vars[[vname]]
  map_dfr(unique(na.omit(df_pt[[col]])), function(lv) {
    tibble(
      variable = vname,
      level = as.character(lv),
      matched_on = TRUE,
      smd = smd_binary(as.integer(df_pt[[col]] == lv), df_pt$tx_mod)
    )
  })
})
balance_unmatched <- tibble(
  variable = c("household_size", "comorbidity", "complication"),
  level = NA_character_,
  matched_on = FALSE,
  smd = c(
    smd_continuous(df_pt$hh_size, df_pt$tx_mod),
    smd_binary(as.integer(df_pt$comorbidity), df_pt$tx_mod),
    smd_binary(as.integer(df_pt$complication), df_pt$tx_mod)
  )
)
s5_table_balance <- bind_rows(balance_matched, balance_unmatched)
print(s5_table_balance, n = Inf)
write_csv(s5_table_balance, file.path(tab_dir, "S5_table_balance_smd.csv"))

# date of notification by arm (date_onset holds the registry notification
# date); HBC was introduced in July 2023
df_pt |>
  group_by(tx_mod) |>
  summarise(
    median_notified = median(date_onset),
    q1 = quantile(date_onset, .25, type = 1),
    q3 = quantile(date_onset, .75, type = 1),
    before_jul_2023 = sum(date_onset < as.Date("2023-07-01"))
  )

# Figure 2: epidemic curve, one panel per arm on a common date axis
epicurve <- df_pt |>
  mutate(week = floor_date(date_onset, unit = "week")) |>
  count(tx_mod, week) |>
  ggplot(aes(week, n, fill = tx_mod)) +
  geom_col(colour = "black") +
  facet_wrap(~tx_mod, ncol = 1) +
  scale_fill_manual(
    name = "Treatment",
    values = c(HBC = "#1b9e77", DTC = "#d95f02")
  ) +
  scale_x_date(date_labels = "%b\n%Y", date_breaks = "1 month") +
  scale_y_continuous(breaks = breaks_pretty()) +
  labs(x = "Date of notification", y = "Number of patients") +
  theme_classic()
ggsave(
  file.path(fig_dir, "figure2_epicurve.jpg"),
  epicurve,
  width = 8,
  height = 5
)

# ============================================================
# 2. Mortality (Table 2)
# ============================================================

df_mort <- df_pt |> filter(!is.na(subclass), !is.na(delay_tx))
nrow(df_mort)
sum(df_mort$died)

# case fatality risk, design-based with subclass as stratum
srv_design |>
  summarise(cfr = survey_mean((died == TRUE) * 100, vartype = "ci"))
srv_design |>
  group_by(tx_mod) |>
  summarise(cfr = survey_mean((died == TRUE) * 100, vartype = "ci"))
df_pt |>
  group_by(age_grp) |>
  summarise(n = n(), deaths = sum(died))

# design-based risk ratios: population-averaged, so unlike the conditional
# models below they include patients in single-arm subclasses
svyglm(died ~ tx_mod, design = srv_design, family = quasipoisson) |>
  tidy(exponentiate = TRUE, conf.int = TRUE)
svyglm(
  died ~ tx_mod + age + vaccinated + delay_tx,
  design = srv_design,
  family = quasipoisson
) |>
  tidy(exponentiate = TRUE, conf.int = TRUE)

# unadjusted, one term at a time, same sample and stratification as the
# adjusted models below
table2_vars <- c("tx_mod", "age", "vaccinated", "delay_tx", "complication")
table2_unadjusted <- map_dfr(table2_vars, \(v) {
  f <- as.formula(paste0("died ~ ", v, " + strata(subclass)"))
  clogit(f, data = df_mort) |>
    tidy(exponentiate = TRUE, conf.int = TRUE) |>
    mutate(variable = v)
})
print(table2_unadjusted)

# Model A: primary (excludes complications, which are treated as an
# intermediate outcome rather than a baseline confounder)
mod_a <- clogit(
  died ~ tx_mod + age + vaccinated + delay_tx + strata(subclass),
  data = df_mort
)
tidy(mod_a, exponentiate = TRUE, conf.int = TRUE)

# Model B: mediation-decomposition sensitivity, additionally adjusting for
# complications
mod_b <- clogit(
  died ~ tx_mod + age + vaccinated + delay_tx + complication + strata(subclass),
  data = df_mort
)
tidy(mod_b, exponentiate = TRUE, conf.int = TRUE)

# forest plot, Model A
tidy(mod_a, exponentiate = TRUE, conf.int = TRUE) |>
  mutate(
    term = case_match(
      term,
      "tx_modHBC" ~ "Treatment: HBC vs DTC",
      "age" ~ "Age (years)",
      "vaccinated" ~ "Vaccinated",
      "delay_tx" ~ "Treatment delay (per day)",
      .default = term
    ),
    term = fct_rev(fct_inorder(term))
  ) |>
  ggplot(aes(x = estimate, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point(size = 3) +
  scale_x_log10(limits = c(0.3, 150)) +
  labs(x = "Adjusted odds ratio (log scale)", y = NULL) +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "forest_plot_mortality.jpg"), width = 8, height = 5)

# delay: linear vs restricted cubic spline (4 knots), to check whether a
# threshold specification would fit better than a linear one
df_sp <- df_pt |>
  filter(!is.na(died), !is.na(delay_tx), !is.na(age), !is.na(vaccinated))
dd <- datadist(df_sp)
options(datadist = "dd")
mod_delay_linear <- lrm(
  died ~ delay_tx + tx_mod + age + vaccinated,
  data = df_sp,
  x = TRUE,
  y = TRUE
)
mod_delay_spline <- lrm(
  died ~ rcs(delay_tx, 4) + tx_mod + age + vaccinated,
  data = df_sp,
  x = TRUE,
  y = TRUE
)
lrtest(mod_delay_linear, mod_delay_spline)

pred_sp <- Predict(
  mod_delay_spline,
  delay_tx = seq(0, 21, by = 0.5),
  tx_mod = "HBC",
  age = median(df_sp$age, na.rm = TRUE),
  vaccinated = 1
)
ggplot(as.data.frame(pred_sp), aes(delay_tx, yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#2c7fb8") +
  geom_line(colour = "#2c7fb8", linewidth = 1) +
  labs(x = "Treatment delay (days)", y = "Log-odds of death") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "spline_delay.png"), width = 8, height = 5)

# calendar time: check whether adding it improves fit over the primary model
mod_no_cal <- clogit(
  died ~ tx_mod + age + vaccinated + delay_tx + strata(subclass),
  data = df_mort
)
mod_with_cal <- clogit(
  died ~ tx_mod +
    age +
    vaccinated +
    delay_tx +
    rcs(as.numeric(date_onset), 4) +
    strata(subclass),
  data = df_mort
)
anova(mod_no_cal, mod_with_cal)
tidy(mod_with_cal, exponentiate = TRUE, conf.int = TRUE)

# epidemic phase: check whether the treatment effect differs before/after
# the median onset date
df_phase <- df_mort |>
  mutate(
    early_phase = as.integer(date_onset <= median(date_onset, na.rm = TRUE))
  )
median(df_phase$date_onset, na.rm = TRUE)
# no HBC deaths in the early phase, so the interaction term is not estimable
df_phase |>
  group_by(early_phase, tx_mod) |>
  summarise(n = n(), deaths = sum(died), .groups = "drop")
mod_phase <- clogit(
  died ~ tx_mod * early_phase + age + vaccinated + delay_ge4 + strata(subclass),
  data = df_phase
)
tidy(mod_phase, exponentiate = TRUE, conf.int = TRUE)

# ============================================================
# 3. Robustness of the treatment-modality estimate
# ============================================================

# excluding the 8 DTC patients with complications
mod_excl_complications <- clogit(
  died ~ tx_mod + age + vaccinated + delay_tx + strata(subclass),
  data = df_mort |> filter(complication == FALSE | is.na(complication))
)
tidy(mod_excl_complications, exponentiate = TRUE, conf.int = TRUE) |>
  filter(term == "tx_modHBC")

# inverse probability of treatment weighting, complementary to the matched
# design (uses the full unmatched cohort).
df_ps <- df_pt |>
  filter(
    !is.na(died),
    !is.na(age),
    !is.na(vaccinated),
    !is.na(delay_tx),
    !is.na(hh_size),
    !is.na(sex)
  ) |>
  mutate(tx_dtc = as.integer(tx_mod == "DTC"), tx_hbc = 1L - tx_dtc)

ps_mod <- glm(
  tx_dtc ~ age + I(age^2) + sex + vaccinated + lga + hh_size,
  data = df_ps,
  family = binomial
)
df_ps <- df_ps |>
  mutate(
    ps = predict(ps_mod, type = "response"),
    p_tx = mean(tx_dtc),
    iptw_stab = if_else(tx_dtc == 1, p_tx / ps, (1 - p_tx) / (1 - ps)),
    iptw_trimmed = pmin(
      pmax(iptw_stab, quantile(iptw_stab, .01)),
      quantile(iptw_stab, .99)
    )
  )

mod_iptw <- glm(
  died ~ tx_hbc,
  data = df_ps,
  family = quasibinomial,
  weights = iptw_trimmed
)
r_iptw <- coeftest(mod_iptw, vcov = vcovHC(mod_iptw, type = "HC1"))
exp(coef(mod_iptw)["tx_hbc"])
exp(coef(mod_iptw)["tx_hbc"] + c(-1, 1) * 1.96 * r_iptw["tx_hbc", "Std. Error"])
r_iptw["tx_hbc", "Pr(>|z|)"]

mod_dr <- glm(
  died ~ tx_hbc + age + vaccinated + delay_tx,
  data = df_ps,
  family = quasibinomial,
  weights = iptw_trimmed
)
r_dr <- coeftest(mod_dr, vcov = vcovHC(mod_dr, type = "HC1"))
exp(coef(mod_dr)["tx_hbc"])
exp(coef(mod_dr)["tx_hbc"] + c(-1, 1) * 1.96 * r_dr["tx_hbc", "Std. Error"])
r_dr["tx_hbc", "Pr(>|z|)"]

# marginal risk ratio (quasi-Poisson) and risk difference (identity link,
# in percentage points), same weights
mod_iptw_rr <- glm(
  died ~ tx_hbc,
  data = df_ps,
  family = quasipoisson,
  weights = iptw_trimmed
)
r_iptw_rr <- coeftest(mod_iptw_rr, vcov = vcovHC(mod_iptw_rr, type = "HC1"))
exp(r_iptw_rr["tx_hbc", 1] + c(0, -1, 1) * 1.96 * r_iptw_rr["tx_hbc", 2])
r_iptw_rr["tx_hbc", 4]

mod_iptw_rd <- glm(
  died ~ tx_hbc,
  data = df_ps,
  family = gaussian,
  weights = iptw_trimmed
)
r_iptw_rd <- coeftest(mod_iptw_rd, vcov = vcovHC(mod_iptw_rd, type = "HC1"))
100 * (r_iptw_rd["tx_hbc", 1] + c(0, -1, 1) * 1.96 * r_iptw_rd["tx_hbc", 2])
r_iptw_rd["tx_hbc", 4]

df_ps |>
  group_by(tx_mod) |>
  summarise(weighted_mortality = 100 * weighted.mean(died, iptw_trimmed))

# balance and overlap diagnostics for the weighting
range(df_ps$ps)
mean(df_ps$ps < 0.05 | df_ps$ps > 0.95)

ps_range_hbc <- range(df_ps$ps[df_ps$tx_dtc == 0])
ps_range_dtc <- range(df_ps$ps[df_ps$tx_dtc == 1])
ps_range_dtc
ps_range_hbc
common_lo <- max(ps_range_hbc[1], ps_range_dtc[1])
common_hi <- min(ps_range_hbc[2], ps_range_dtc[2])
sum(df_ps$ps < common_lo | df_ps$ps > common_hi)

ggplot(df_ps, aes(ps, fill = tx_mod)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(
    values = c(HBC = "#1b9e77", DTC = "#7570b3"),
    name = "Treatment"
  ) +
  labs(x = "Propensity score (P[DTC])", y = "Density") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "iptw_overlap.jpg"), width = 7, height = 5)

ess <- function(w) sum(w)^2 / sum(w^2)
df_ps |>
  group_by(tx_mod) |>
  summarise(
    n = n(),
    ess_untrimmed = ess(iptw_stab),
    ess_trimmed = ess(iptw_trimmed)
  )

for (v in c("age", "vaccinated", "hh_size")) {
  mn1 <- mean(df_ps[[v]][df_ps$tx_dtc == 1], na.rm = TRUE)
  mn0 <- mean(df_ps[[v]][df_ps$tx_dtc == 0], na.rm = TRUE)
  s <- sqrt(
    (var(df_ps[[v]][df_ps$tx_dtc == 1], na.rm = TRUE) +
      var(df_ps[[v]][df_ps$tx_dtc == 0], na.rm = TRUE)) /
      2
  )
  wmn1 <- weighted.mean(
    df_ps[[v]][df_ps$tx_dtc == 1],
    df_ps$iptw_trimmed[df_ps$tx_dtc == 1],
    na.rm = TRUE
  )
  wmn0 <- weighted.mean(
    df_ps[[v]][df_ps$tx_dtc == 0],
    df_ps$iptw_trimmed[df_ps$tx_dtc == 0],
    na.rm = TRUE
  )
  cat(sprintf(
    "%s: raw SMD=%.3f, weighted SMD=%.3f\n",
    v,
    (mn1 - mn0) / s,
    (wmn1 - wmn0) / s
  ))
}

trim_fit <- function(probs) {
  b <- quantile(df_ps$iptw_stab, probs)
  w <- pmin(pmax(df_ps$iptw_stab, b[1]), b[2])
  m <- glm(died ~ tx_hbc, data = df_ps, family = quasibinomial, weights = w)
  r <- coeftest(m, vcov = vcovHC(m, type = "HC1"))
  c(
    or = exp(coef(m)["tx_hbc"]),
    lo = exp(coef(m)["tx_hbc"] - 1.96 * r["tx_hbc", "Std. Error"]),
    hi = exp(coef(m)["tx_hbc"] + 1.96 * r["tx_hbc", "Std. Error"]),
    p = r["tx_hbc", "Pr(>|z|)"]
  )
}
trim_fit(c(0, 1))
trim_fit(c(.05, .95))
trim_fit(c(.01, .99))

# selection-IPW: reweight the IPTW estimate for differential tracing
# success by sex/vaccination status found in the attrition check above.
# Replacements take the probability predicted from their own registry
# values, assuming they resemble originals with the same values.
sel_mod <- glm(
  traced ~ sex + age_group + vaccination_status + tx_location,
  data = df_attr,
  family = binomial
)
marg_traced <- mean(df_attr$traced)

df_ps_sel <- df_ps |>
  mutate(
    is_replacement = !(case_id %in% df_cohort_samp$case_id),
    p_traced = predict(
      sel_mod,
      newdata = tibble(
        sex = sex_reg,
        age_group = age_grp_reg,
        vaccination_status = vacc_status_reg,
        tx_location = if_else(tx_mod == "HBC", "home", "hosp")
      ),
      type = "response"
    ),
    sel_weight = marg_traced / p_traced,
    combined_weight = iptw_trimmed * sel_weight
  )
range(df_ps_sel$sel_weight)
df_ps_sel |>
  group_by(is_replacement) |>
  summarise(mean_weight = mean(sel_weight))

mod_sel_marginal <- glm(
  died ~ tx_hbc,
  data = df_ps_sel,
  family = quasibinomial,
  weights = combined_weight
)
r_sel_marginal <- coeftest(
  mod_sel_marginal,
  vcov = vcovHC(mod_sel_marginal, type = "HC1")
)
exp(
  r_sel_marginal["tx_hbc", 1] + c(0, -1, 1) * 1.96 * r_sel_marginal["tx_hbc", 2]
)
r_sel_marginal["tx_hbc", 4]

mod_sel_iptw <- glm(
  died ~ tx_hbc + age + vaccinated + delay_tx,
  data = df_ps_sel,
  family = quasibinomial,
  weights = combined_weight
)
r_sel <- coeftest(mod_sel_iptw, vcov = vcovHC(mod_sel_iptw, type = "HC1"))
exp(coef(mod_sel_iptw)["tx_hbc"])
exp(
  coef(mod_sel_iptw)["tx_hbc"] + c(-1, 1) * 1.96 * r_sel["tx_hbc", "Std. Error"]
)
r_sel["tx_hbc", "Pr(>|z|)"]

# originally sampled patients only (no replacements): primary model, Firth
# and IPTW, with the propensity model re-estimated in this subset
df_orig <- df_mort |> filter(case_id %in% df_cohort_samp$case_id)
nrow(df_orig)
sum(df_orig$died)
mod_orig <- clogit(
  died ~ tx_mod + age + vaccinated + delay_tx + strata(subclass),
  data = df_orig
)
tidy(mod_orig, exponentiate = TRUE, conf.int = TRUE) |>
  filter(term == "tx_modHBC")

# deaths in subclasses that inform the conditional estimate (both arms
# present and the outcome varies)
df_orig |>
  group_by(subclass) |>
  filter(n_distinct(tx_mod) == 2, any(died), !all(died)) |>
  ungroup() |>
  summarise(
    n_subclasses = n_distinct(subclass),
    deaths = sum(died),
    deaths_hbc = sum(died & tx_mod == "HBC")
  )

mod_firth_orig <- logistf::logistf(
  died ~ tx_mod + age + vaccinated + delay_tx,
  data = df_orig,
  firth = TRUE
)
c(
  or = exp(coef(mod_firth_orig)[["tx_modHBC"]]),
  lo = exp(mod_firth_orig$ci.lower[["tx_modHBC"]]),
  hi = exp(mod_firth_orig$ci.upper[["tx_modHBC"]]),
  p = mod_firth_orig$prob[["tx_modHBC"]]
)

df_ps_orig <- df_ps |> filter(case_id %in% df_cohort_samp$case_id)
ps_mod_orig <- glm(
  tx_dtc ~ age + I(age^2) + sex + vaccinated + lga + hh_size,
  data = df_ps_orig,
  family = binomial
)
df_ps_orig <- df_ps_orig |>
  mutate(
    ps = predict(ps_mod_orig, type = "response"),
    p_tx = mean(tx_dtc),
    iptw_stab = if_else(tx_dtc == 1, p_tx / ps, (1 - p_tx) / (1 - ps)),
    iptw_trimmed = pmin(
      pmax(iptw_stab, quantile(iptw_stab, .01)),
      quantile(iptw_stab, .99)
    )
  )
mod_iptw_orig <- glm(
  died ~ tx_hbc,
  data = df_ps_orig,
  family = quasibinomial,
  weights = iptw_trimmed
)
r_iptw_orig <- coeftest(mod_iptw_orig, vcov = vcovHC(mod_iptw_orig, type = "HC1"))
exp(r_iptw_orig["tx_hbc", 1] + c(0, -1, 1) * 1.96 * r_iptw_orig["tx_hbc", 2])

# Firth-penalised logistic regression, small-sample bias correction given
# the low events-per-variable ratio (24 deaths / 4 covariates). Firth's
# penalised likelihood doesn't support conditioning on matched strata, so
# this adjusts directly for age and vaccination instead.
mod_firth <- logistf::logistf(
  died ~ tx_mod + age + vaccinated + delay_tx,
  data = df_mort,
  firth = TRUE
)
tibble(
  term = names(coef(mod_firth))[-1],
  aOR = exp(coef(mod_firth))[-1],
  lo = exp(mod_firth$ci.lower)[-1],
  hi = exp(mod_firth$ci.upper)[-1],
  p = mod_firth$prob[-1]
)

# E-values: minimum strength of unmeasured confounding that would explain
# away each estimate. evalues.OR() expects an OR > 1, so the two
# protective factors (treatment, vaccination) are inverted first.
est_a <- tidy(mod_a, exponentiate = TRUE, conf.int = TRUE)
tx_row <- est_a |> filter(term == "tx_modHBC")
vacc_row <- est_a |> filter(term == "vaccinated")
delay_row <- est_a |> filter(term == "delay_tx")

evalues.OR(
  est = 1 / tx_row$estimate,
  lo = 1 / tx_row$conf.high,
  hi = 1 / tx_row$conf.low,
  rare = TRUE
)
evalues.OR(
  est = 1 / vacc_row$estimate,
  lo = 1 / vacc_row$conf.high,
  hi = 1 / vacc_row$conf.low,
  rare = TRUE
)
evalues.OR(
  est = delay_row$estimate,
  lo = delay_row$conf.low,
  hi = delay_row$conf.high,
  rare = TRUE
)

# multicollinearity, discrimination, calibration, influence - clogit's
# conditional likelihood has no absolute intercept, so these are computed
# on an equivalent unconditional model
mod_unconditional <- glm(
  died ~ tx_mod + age + vaccinated + delay_tx,
  data = df_mort,
  family = binomial
)
car::vif(mod_unconditional)

pred <- predict(mod_unconditional, type = "response")
obs <- mod_unconditional$model$died
r <- rank(pred)
n1 <- sum(obs)
n0 <- sum(!obs)
auc <- (sum(r[obs]) - n1 * (n1 + 1) / 2) / (n1 * n0)
auc

# within-subclass concordance of the conditional model itself
conc_a <- concordance(mod_a)
c(concordance = conc_a$concordance, se = sqrt(conc_a$var))

tibble(pred = pred, obs = as.integer(obs)) |>
  mutate(quintile = ntile(pred, 5)) |>
  group_by(quintile) |>
  summarise(n = n(), mean_predicted = mean(pred), observed_rate = mean(obs))

# dfbeta isn't available for method="exact" (the primary model's method) -
# an efron-method twin is fit purely for this diagnostic
mod_diag <- clogit(
  died ~ tx_mod + age + vaccinated + delay_tx + strata(subclass),
  data = df_mort,
  method = "efron"
)
db <- residuals(mod_diag, type = "dfbeta")
colnames(db) <- names(coef(mod_diag))
as.data.frame(db) |>
  mutate(row = row_number()) |>
  arrange(desc(abs(tx_modHBC))) |>
  head(5)

mod_diag_b <- clogit(
  died ~ tx_mod + age + vaccinated + delay_tx + complication + strata(subclass),
  data = df_mort,
  method = "efron"
)
db_b <- residuals(mod_diag_b, type = "dfbeta")
colnames(db_b) <- names(coef(mod_diag_b))
as.data.frame(db_b) |>
  mutate(row = row_number()) |>
  arrange(desc(abs(complicationTRUE))) |>
  head(5)

# LGA: exact-matching variable, so adding it as a covariate should be
# aliased within-stratum; checked directly rather than assumed. A random-
# intercept model (a different, non-matched design) is reported as a
# supportive check.
mod_lga_fixed <- clogit(
  died ~ tx_mod + age + vaccinated + delay_tx + factor(lga) + strata(subclass),
  data = df_mort
)
tidy(mod_lga_fixed, exponentiate = TRUE, conf.int = TRUE) |>
  filter(term == "tx_modHBC")

mod_lga_random <- lme4::glmer(
  died ~ tx_mod + age + vaccinated + delay_tx + (1 | lga),
  data = df_mort,
  family = binomial,
  control = lme4::glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)
broom.mixed::tidy(mod_lga_random, exponentiate = TRUE, conf.int = TRUE) |>
  filter(term == "tx_modHBC")
lme4::isSingular(mod_lga_random)

# card-confirmed vaccination, in case self-report alone is biased by recall
df_card <- df_pt |>
  filter(!is.na(vacc_card_confirmed), !is.na(subclass), !is.na(delay_tx))
nrow(df_card)
sum(df_card$vacc_card_confirmed == 1, na.rm = TRUE)
mod_card <- clogit(
  died ~ tx_mod + age + vacc_card_confirmed + delay_tx + strata(subclass),
  data = df_card
)
tidy(mod_card, exponentiate = TRUE, conf.int = TRUE)

# stricter definition: card seen vs self-reported unvaccinated, everyone
# else excluded
df_card_seen <- df_pt |>
  filter(!is.na(vacc_card_seen), !is.na(subclass), !is.na(delay_tx))
nrow(df_card_seen)
sum(df_card_seen$died)
mod_card_seen <- clogit(
  died ~ tx_mod + age + vacc_card_seen + delay_tx + strata(subclass),
  data = df_card_seen
)
tidy(mod_card_seen, exponentiate = TRUE, conf.int = TRUE)

# effect modification by treatment delay (days), against Model A
mod_tx_delay <- clogit(
  died ~ tx_mod * delay_tx + age + vaccinated + strata(subclass),
  data = df_mort
)
tidy(mod_tx_delay, exponentiate = TRUE, conf.int = TRUE) |>
  filter(term == "tx_modHBC:delay_tx")
anova(mod_a, mod_tx_delay)

# ============================================================
# 4. Long-term sequelae
# ============================================================

# survivors only, Model A covariates and stratification. Complications are
# left out as in Model A; all five complicated survivors were DTC, so they
# can't be estimated in the conditional model anyway.
df_seq <- df_mort |>
  filter(!died) |>
  mutate(has_sequelae = as.integer(!is.na(sequelae)))
df_seq |>
  group_by(tx_mod) |>
  summarise(n = n(), sequelae = sum(has_sequelae))

# subclasses that inform the conditional estimate (both arms present and
# the outcome varies)
df_seq |>
  group_by(subclass) |>
  filter(n_distinct(tx_mod) == 2, any(has_sequelae == 1), !all(has_sequelae == 1)) |>
  ungroup() |>
  summarise(
    n_subclasses = n_distinct(subclass),
    events_dtc = sum(has_sequelae[tx_mod == "DTC"]),
    events_hbc = sum(has_sequelae[tx_mod == "HBC"])
  )
n_distinct(df_seq$subclass)

mod_sequelae <- clogit(
  has_sequelae ~ tx_mod + age + vaccinated + delay_tx + strata(subclass),
  data = df_seq
)
tidy(mod_sequelae, exponentiate = TRUE, conf.int = TRUE)

# sensitivity analyses, given only 14 events
clogit(has_sequelae ~ tx_mod + strata(subclass), data = df_seq) |>
  tidy(exponentiate = TRUE, conf.int = TRUE)

firth_or <- function(m) {
  tibble(
    term = names(coef(m))[-1],
    or = exp(coef(m))[-1],
    lo = exp(m$ci.lower)[-1],
    hi = exp(m$ci.upper)[-1],
    p = m$prob[-1]
  )
}
logistf::logistf(
  has_sequelae ~ tx_mod + age + vaccinated + delay_tx,
  data = df_seq,
  firth = TRUE
) |>
  firth_or()
logistf::logistf(
  has_sequelae ~ tx_mod + age + vaccinated + delay_tx + complication,
  data = df_seq,
  firth = TRUE
) |>
  firth_or()

df_seq_uncompl <- df_seq |> filter(!complication)
nrow(df_seq_uncompl)
clogit(
  has_sequelae ~ tx_mod + age + vaccinated + delay_tx + strata(subclass),
  data = df_seq_uncompl
) |>
  tidy(exponentiate = TRUE, conf.int = TRUE)

# design-based risk ratio and risk difference (percentage points), subclass
# as stratum
srv_seq <- df_seq |> as_survey_design(strata = subclass, ids = NULL)
svyglm(
  has_sequelae ~ tx_mod + age + vaccinated + delay_tx,
  design = srv_seq,
  family = quasipoisson
) |>
  tidy(exponentiate = TRUE, conf.int = TRUE) |>
  filter(term == "tx_modHBC")
svyglm(
  has_sequelae ~ tx_mod + age + vaccinated + delay_tx,
  design = srv_seq
) |>
  tidy(conf.int = TRUE) |>
  filter(term == "tx_modHBC") |>
  mutate(across(c(estimate, conf.low, conf.high), \(x) 100 * x))

# composite of death or long-term sequelae among all interviewed patients
df_composite <- df_mort |>
  mutate(death_or_sequelae = died | !is.na(sequelae))
df_composite |>
  group_by(tx_mod) |>
  summarise(n = n(), events = sum(death_or_sequelae))
clogit(
  death_or_sequelae ~ tx_mod + age + vaccinated + delay_tx + strata(subclass),
  data = df_composite
) |>
  tidy(exponentiate = TRUE, conf.int = TRUE)
logistf::logistf(
  death_or_sequelae ~ tx_mod + age + vaccinated + delay_tx,
  data = df_composite,
  firth = TRUE
) |>
  firth_or()

# ============================================================
# 5. Household transmission (Table 3)
# ============================================================

sar_ci <- function(cases, n) {
  bt <- binom.test(cases, n)
  c(sar = cases / n, lo = bt$conf.int[1], hi = bt$conf.int[2])
}

sar_ci(
  sum(
    df_hh$hh_case == 1 & df_hh$timing_days >= 2 & df_hh$timing_days <= 14,
    na.rm = TRUE
  ),
  nrow(df_hh)
)
sar_ci(
  sum(
    df_hh$hh_case == 1 & df_hh$timing_days >= 2 & df_hh$timing_days <= 30,
    na.rm = TRUE
  ),
  nrow(df_hh)
)
sar_ci(sum(df_hh$hh_case == 1, na.rm = TRUE), nrow(df_hh))

sar_tx <- df_hh |>
  filter(!is.na(tx_mod)) |>
  group_by(tx_mod) |>
  summarise(
    n = n(),
    sec_14 = sum(
      hh_case == 1 & timing_days >= 2 & timing_days <= 14,
      na.rm = TRUE
    ),
    sec_30 = sum(
      hh_case == 1 & timing_days >= 2 & timing_days <= 30,
      na.rm = TRUE
    ),
    sec_any = sum(hh_case == 1, na.rm = TRUE)
  )
sar_tx

# Table 3 by arm, exact binomial CIs
sar_tx |>
  pivot_longer(starts_with("sec_"), names_to = "window", values_to = "cases") |>
  mutate(ci = map2(cases, n, sar_ci)) |>
  unnest_wider(ci)

# index households with a secondary case in the primary window
df_hh |>
  group_by(tx_mod, parent_index) |>
  summarise(
    case_14 = any(hh_case == 1 & timing_days >= 2 & timing_days <= 14, na.rm = TRUE),
    .groups = "drop_last"
  ) |>
  summarise(households = n(), with_case = sum(case_14))

# multivariable model, restricted to the primary 2-14 day window, SEs
# clustered by index-case household
df_hh_mv <- df_hh |>
  filter(
    !is.na(tx_mod),
    !is.na(hh_case),
    !is.na(hh_size),
    !is.na(age_hh_f),
    !is.na(sex_hh_f),
    hh_case == 0 | (!is.na(timing_days) & timing_days >= 2 & timing_days <= 14)
  )
nrow(df_hh_mv)
sum(df_hh_mv$hh_case)

sar_vars <- c(
  "tx_mod",
  "complication",
  "delay_ge4",
  "vaccinated",
  "age_hh_f",
  "sex_hh_f",
  "hh_size",
  "hh_vacc"
)
sar_univariable <- map_dfr(sar_vars, \(v) {
  m <- glm(
    as.formula(paste0("hh_case ~ ", v)),
    data = df_hh_mv,
    family = binomial
  )
  r <- coeftest(m, vcov = vcovCL(m, cluster = df_hh_mv$parent_index))
  tibble(
    term = rownames(r)[-1],
    or = exp(r[-1, 1]),
    lo = exp(r[-1, 1] - 1.96 * r[-1, 2]),
    hi = exp(r[-1, 1] + 1.96 * r[-1, 2]),
    p = r[-1, 4]
  )
})
print(sar_univariable, n = Inf)

mod_sar <- glm(
  hh_case ~ tx_mod +
    complication +
    delay_ge4 +
    vaccinated +
    age_hh_f +
    sex_hh_f +
    hh_size +
    hh_vacc,
  data = df_hh_mv,
  family = binomial
)
r_sar <- coeftest(
  mod_sar,
  vcov = vcovCL(mod_sar, cluster = df_hh_mv$parent_index)
)
sar_multivariable <- tibble(
  term = rownames(r_sar)[-1],
  or = exp(r_sar[-1, 1]),
  lo = exp(r_sar[-1, 1] - 1.96 * r_sar[-1, 2]),
  hi = exp(r_sar[-1, 1] + 1.96 * r_sar[-1, 2]),
  p = r_sar[-1, 4]
)
print(sar_multivariable, n = Inf)

# risk ratios alongside the ORs above - SAR isn't rare enough for the two
# to be assumed equivalent
mod_sar_rr <- glm(
  hh_case ~ tx_mod +
    complication +
    delay_ge4 +
    vaccinated +
    age_hh_f +
    sex_hh_f +
    hh_size +
    hh_vacc,
  data = df_hh_mv,
  family = poisson
)
r_rr <- coeftest(
  mod_sar_rr,
  vcov = vcovCL(mod_sar_rr, cluster = df_hh_mv$parent_index)
)
tibble(
  term = rownames(r_rr)[-1],
  rr = exp(r_rr[-1, 1]),
  lo = exp(r_rr[-1, 1] - 1.96 * r_rr[-1, 2]),
  hi = exp(r_rr[-1, 1] + 1.96 * r_rr[-1, 2])
)

# multilevel model, household random intercept, as a check on within-
# household correlation
mod_sar_multilevel <- lme4::glmer(
  hh_case ~ tx_mod +
    complication +
    delay_ge4 +
    vaccinated +
    hh_size +
    hh_vacc +
    (1 | parent_index),
  data = df_hh_mv,
  family = binomial,
  control = lme4::glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)
broom.mixed::tidy(mod_sar_multilevel, exponentiate = TRUE, conf.int = TRUE)
vc <- as.data.frame(lme4::VarCorr(mod_sar_multilevel))
vc$vcov[1] / (vc$vcov[1] + pi^2 / 3) # household-level ICC

# window sensitivity: same model spec, fit on the 2-30-day and unrestricted
# windows
fit_sar_window <- function(hi_day, restrict_window = TRUE) {
  d <- df_hh |>
    filter(
      !is.na(tx_mod),
      !is.na(hh_case),
      !is.na(hh_size),
      !is.na(age_hh_f),
      !is.na(sex_hh_f)
    )
  if (restrict_window) {
    d <- d |>
      filter(
        hh_case == 0 |
          (!is.na(timing_days) & timing_days >= 2 & timing_days <= hi_day)
      )
  }
  m <- glm(
    hh_case ~ tx_mod +
      complication +
      delay_ge4 +
      vaccinated +
      age_hh_f +
      sex_hh_f +
      hh_size +
      hh_vacc,
    data = d,
    family = binomial
  )
  r <- coeftest(m, vcov = vcovCL(m, cluster = d$parent_index))
  tibble(
    n = nrow(d),
    cases = sum(d$hh_case),
    aor = exp(r["tx_modHBC", 1]),
    lo = exp(r["tx_modHBC", 1] - 1.96 * r["tx_modHBC", 2]),
    hi = exp(r["tx_modHBC", 1] + 1.96 * r["tx_modHBC", 2]),
    p = r["tx_modHBC", 4]
  )
}
fit_sar_window(30)
fit_sar_window(NA, restrict_window = FALSE)

# ============================================================
# 6. Acceptability
# ============================================================

df_pt |>
  group_by(tx_mod) |>
  summarise(
    across(
      c(score_improve, score_acceptable, score_quality, score_mean),
      list(
        mean = \(x) mean(x, na.rm = TRUE),
        sd = \(x) sd(x, na.rm = TRUE),
        median = \(x) median(x, na.rm = TRUE)
      )
    )
  ) |>
  pivot_longer(-tx_mod) |>
  print(n = Inf)

for (v in c("score_improve", "score_acceptable", "score_quality")) {
  wt <- wilcox.test(
    df_pt[[v]][df_pt$tx_mod == "HBC"],
    df_pt[[v]][df_pt$tx_mod == "DTC"]
  )
  cat(sprintf("%s: W=%.0f, p=%.4f\n", v, wt$statistic, wt$p.value))
}

# linear mixed-effects models, matched subclass as a random intercept, one
# per dimension and for the composite. Identity link, so the coefficients
# are mean differences in Likert points (HBC minus DTC), not exponentiated.
acceptability_lmm <- map_dfr(
  c("score_acceptable", "score_quality", "score_improve", "score_mean"),
  \(v) {
    m <- lmer(
      as.formula(paste0(v, " ~ tx_mod + delay_ge4 + complication + (1 | subclass)")),
      data = df_pt
    )
    tidy(m, effects = "fixed", conf.int = TRUE) |>
      filter(term != "(Intercept)") |>
      mutate(
        outcome = v,
        n = nobs(m),
        subclass_var = as.data.frame(VarCorr(m))$vcov[1]
      )
  }
)
print(acceptability_lmm, n = Inf, width = Inf)

# ordinal (cumulative-link) sensitivity models, all three Likert dimensions
df_pt_ord <- df_pt |>
  mutate(
    score_acceptable_ord = factor(score_acceptable, ordered = TRUE),
    score_quality_ord = factor(score_quality, ordered = TRUE),
    score_improve_ord = factor(score_improve, ordered = TRUE)
  )

fit_ordinal <- function(outcome_var, label) {
  m <- ordinal::clm(
    as.formula(paste0(outcome_var, " ~ tx_mod + delay_ge4 + complication")),
    data = df_pt_ord
  )
  s <- summary(m)$coefficients
  s <- s[!str_detect(rownames(s), "\\|"), , drop = FALSE]
  tibble(
    dimension = label,
    term = rownames(s),
    OR = exp(s[, "Estimate"]),
    lo = exp(s[, "Estimate"] - 1.96 * s[, "Std. Error"]),
    hi = exp(s[, "Estimate"] + 1.96 * s[, "Std. Error"]),
    p = s[, "Pr(>|z|)"]
  )
}
ordinal_results <- bind_rows(
  fit_ordinal("score_acceptable_ord", "Acceptability"),
  fit_ordinal("score_quality_ord", "Quality"),
  fit_ordinal("score_improve_ord", "Perceived improvement")
)
print(ordinal_results, n = Inf)

# a subclass random intercept can't be estimated in these models (variance
# ~0), so CIs for treatment and delay come from resampling subclasses
# instead: 2,000 replicates, percentile intervals, non-converged fits
# dropped. Complication is too unstable in resamples to bootstrap.
set.seed(2026)
boot_ordinal <- function(outcome_var, label, n_boot = 2000) {
  d <- df_pt_ord |> filter(!is.na(subclass))
  subclasses <- unique(d$subclass)
  rows <- split(seq_len(nrow(d)), d$subclass)
  f <- as.formula(paste0(outcome_var, " ~ tx_mod + delay_ge4 + complication"))
  map_dfr(seq_len(n_boot), \(b) {
    sampled <- sample(subclasses, length(subclasses), replace = TRUE)
    fit <- tryCatch(
      suppressWarnings(
        ordinal::clm(f, data = d[unlist(rows[as.character(sampled)]), ])
      ),
      error = \(e) NULL
    )
    tibble(
      dimension = label,
      converged = !is.null(fit) && fit$convergence$code == 0,
      tx_modHBC = if (is.null(fit)) NA else coef(fit)["tx_modHBC"],
      delay_ge4TRUE = if (is.null(fit)) NA else coef(fit)["delay_ge4TRUE"]
    )
  })
}
ordinal_boot <- bind_rows(
  boot_ordinal("score_acceptable_ord", "Acceptability"),
  boot_ordinal("score_quality_ord", "Quality"),
  boot_ordinal("score_improve_ord", "Perceived improvement")
)
ordinal_boot |> count(dimension, converged)

ordinal_boot_ci <- ordinal_boot |>
  filter(converged) |>
  pivot_longer(c(tx_modHBC, delay_ge4TRUE), names_to = "term") |>
  group_by(dimension, term) |>
  summarise(
    boot_lo = exp(quantile(value, .025, na.rm = TRUE)),
    boot_hi = exp(quantile(value, .975, na.rm = TRUE)),
    .groups = "drop"
  )
ordinal_results <- ordinal_results |>
  left_join(ordinal_boot_ci, by = c("dimension", "term"))
print(ordinal_results, n = Inf)

# ============================================================
# 7. Supplementary tables
# ============================================================

fmt_ci <- function(est, lo, hi, digits = 2) {
  sprintf("%.*f [%.*f,%.*f]", digits, est, digits, lo, digits, hi)
}
gt_style <- function(g) {
  g |> tab_options(table.font.name = "Times New Roman", table.font.size = 11)
}

sar_label <- function(term) {
  case_match(
    term,
    "tx_modHBC" ~ "Treatment: HBC vs DTC",
    "complicationTRUE" ~ "Complication present",
    "delay_ge4TRUE" ~ "Treatment delay ≥4 days",
    "vaccinated" ~ "Index case vaccinated",
    "age_hh_f5 - 14 years" ~ "Contact age 5-14y (vs 0-4y)",
    "age_hh_f15+ years" ~ "Contact age 15+y (vs 0-4y)",
    "sex_hh_fMale" ~ "Contact sex: male",
    "hh_size" ~ "Household size",
    "hh_vacc" ~ "Contact vaccinated",
    .default = term
  )
}

s_table_sar_rr <- sar_multivariable |>
  select(term, or, or_lo = lo, or_hi = hi) |>
  left_join(
    tibble(
      term = rownames(r_rr)[-1],
      rr = exp(r_rr[-1, 1]),
      rr_lo = exp(r_rr[-1, 1] - 1.96 * r_rr[-1, 2]),
      rr_hi = exp(r_rr[-1, 1] + 1.96 * r_rr[-1, 2])
    ),
    by = "term"
  ) |>
  transmute(
    `Risk factor` = sar_label(term),
    `Adjusted OR (95% CI)` = fmt_ci(or, or_lo, or_hi),
    `Risk ratio (95% CI)` = fmt_ci(rr, rr_lo, rr_hi)
  )
gt(s_table_sar_rr) |>
  tab_header(
    title = md(
      "**S3 Table. Odds ratios and risk ratios for household secondary attack**"
    )
  ) |>
  gt_style() |>
  gtsave(file.path(tab_dir, "S3_table_SAR_risk_ratios.docx"))

s_table_iptw <- bind_rows(
  tibble(
    Section = "Propensity score",
    Metric = "Range (both arms)",
    Value = sprintf("%.3f-%.3f", min(df_ps$ps), max(df_ps$ps))
  ),
  tibble(
    Section = "Propensity score",
    Metric = "Proportion PS <0.05 or >0.95",
    Value = sprintf("%.1f%%", 100 * mean(df_ps$ps < 0.05 | df_ps$ps > 0.95))
  ),
  tibble(
    Section = "Propensity score",
    Metric = "Patients outside common support",
    Value = sprintf(
      "%d/%d",
      sum(df_ps$ps < common_lo | df_ps$ps > common_hi),
      nrow(df_ps)
    )
  )
)
gt(s_table_iptw, groupname_col = "Section") |>
  tab_header(title = md("**S4 Table. Propensity score diagnostics**")) |>
  gt_style() |>
  gtsave(file.path(tab_dir, "S4_table_IPTW_diagnostics.docx"))

# bootstrap CIs for treatment and delay, model-based (Wald) CIs for
# complication
s_table_ordinal <- ordinal_results |>
  transmute(
    Dimension = dimension,
    Term = term,
    `OR (95% CI)` = if_else(
      is.na(boot_lo),
      paste(fmt_ci(OR, lo, hi), "†"),
      fmt_ci(OR, boot_lo, boot_hi)
    ),
    p = sprintf("%.3f", p)
  )
gt(s_table_ordinal, groupname_col = "Dimension") |>
  tab_header(
    title = md("**S6 Table. Ordinal (cumulative-link) acceptability models**")
  ) |>
  tab_source_note(
    "CIs are subclass-bootstrap percentile intervals, except † model-based Wald intervals; p-values are model-based."
  ) |>
  gt_style() |>
  gtsave(file.path(tab_dir, "S6_table_ordinal_acceptability.docx"))
