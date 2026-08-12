# Community burden of diphtheria during the 2023-24 epidemic in Kano
# State, Nigeria: population-based household survey.
#
# Reproduces every number, table and figure reported in the manuscript,
# in the order they appear (Participants -> Community attack rate ->
# Factors associated with diphtheria -> Mortality -> Health-seeking and
# treatment delay -> Vaccination and mortality -> Vaccination sensitivity
# analyses -> Vaccination coverage and campaign performance -> spatial
# figures and correlations). Input is the frozen, de-identified survey
# dataset; see README for how to request it.

pacman::p_load(survey, srvyr, EValue, scales, sf, spatstat, stars, mgcv, metR, patchwork, tidyverse)

fig_dir <- "output/figures"
tab_dir <- "output/tables"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 0. Data
# ============================================================
# df_ind: one row per surveyed individual (n=7,998), unweighted - used for
#   raw counts (Table 1) and as the basis for the survey design below
# df_hh: one row per approached household (n=1,083) - used for the
#   response-rate figure
# map_boundaries: LGA and study-area boundary polygons, used only by the
#   spatial figures at the end of this script

df_ind <- read_rds("data/df_ind.rds")
df_hh <- read_rds("data/df_hh.rds")
maps <- read_rds("data/map_boundaries.rds")
map_nga_ward <- maps$map_nga_ward
map_nga_ward_union <- maps$map_nga_ward_union

options(survey.lonely.psu = "adjust")

df_srvy_ind <-
  df_ind |>
  filter(!is.na(cluster), !is.na(index_hh)) |>
  as_survey(ids = c(cluster, index_hh), strata = lga, weights = weight, nest = TRUE)

# ============================================================
# 1. Participants (Table 1)
# ============================================================

nrow(df_hh)
sum(df_hh$consent == "Yes")
mean(df_hh$consent == "Yes")

nrow(df_ind)
sum(df_ind$died)

df_ind |>
  summarise(pct_female = mean(sex == "Female") * 100, median_age = median(age), q1 = quantile(age, .25), q3 = quantile(age, .75))

table1_vars <- c("lga", "age_grp_diph", "sex", "vacc", "died")
table_1 <-
  table1_vars |>
  set_names() |>
  map(\(v) {
    df_ind |>
      count(diphtheria, across(all_of(v))) |>
      group_by(diphtheria) |>
      mutate(pct = n / sum(n) * 100) |>
      rename(level = 2) |>
      mutate(level = as.character(level))
  }) |>
  list_rbind(names_to = "variable")
print(table_1, n = Inf)

table1_chisq <- map_dfr(table1_vars, \(v) {
  test <- survey::svychisq(as.formula(paste0("~", v, "+diphtheria")), design = df_srvy_ind, statistic = "Chisq")
  tibble(variable = v, p = test$p.value)
})
print(table1_chisq)

write_csv(table_1, file.path(tab_dir, "table1_counts.csv"))
write_csv(table1_chisq, file.path(tab_dir, "table1_chisq.csv"))

# ============================================================
# 2. Community attack rate and comparison with facility surveillance
# ============================================================

df_srvy_ind |>
  summarise(attack_rate = survey_mean(diphtheria * 100, vartype = "ci", deff = TRUE, na.rm = TRUE))

# incidence rate per 1,000 person-years (case count / total person-time)
df_srvy_ind |>
  summarise(
    ir = survey_ratio(numerator = diphtheria * 1000 * 365.25, denominator = pt, vartype = "ci", na.rm = TRUE)
  )

# facility-reported AR for the same 4 LGAs (Abbas et al, ref 13): 262 per
# 100,000 = 0.262%, a fixed external denominator, not survey data
facility_ar <- 262 / 100000
ar_est <- df_srvy_ind |> summarise(ar = survey_mean(diphtheria, vartype = "ci", na.rm = TRUE))
ar_est$ar / facility_ar
ar_est$ar_low / facility_ar
ar_est$ar_upp / facility_ar

df_srvy_ind |>
  group_by(age_grp_diph) |>
  summarise(attack_rate = survey_mean(diphtheria * 100, vartype = "ci", na.rm = TRUE))

df_srvy_ind |>
  group_by(sex) |>
  summarise(attack_rate = survey_mean(diphtheria * 100, vartype = "ci", na.rm = TRUE))

# household secondary attack rate: among households with >=1 case, how
# many additional ("secondary") cases among the remaining household
# members
hh_summary <-
  df_ind |>
  group_by(index_hh, cluster, lga, weight) |>
  summarise(hh_size = n(), num_cases = sum(diphtheria, na.rm = TRUE), .groups = "drop") |>
  filter(num_cases > 0) |>
  mutate(secondary_cases = pmax(num_cases - 1, 0), exposed = pmax(hh_size - 1, 0), sar = secondary_cases / exposed)

nrow(hh_summary)
sum(hh_summary$num_cases >= 2)

hh_design <- hh_summary |> as_survey(ids = c(cluster, index_hh), strata = lga, weights = weight, nest = TRUE)
hh_design |> summarise(SAR = survey_mean(sar * 100, vartype = "ci", deff = TRUE, na.rm = TRUE))

# case-definition specificity sensitivity: restrict to cases additionally
# reported as diagnosed by a health professional
df_ind_corrob <- df_ind |> mutate(diphtheria_corrob = diphtheria == TRUE & diagnosis == "Yes")
df_srvy_corrob <-
  df_ind_corrob |>
  filter(!is.na(cluster), !is.na(index_hh)) |>
  as_survey(ids = c(cluster, index_hh), strata = lga, weights = weight, nest = TRUE)

df_srvy_corrob |> summarise(ar = survey_mean(diphtheria == TRUE, vartype = "ci", na.rm = TRUE))
df_srvy_corrob |> summarise(ar = survey_mean(diphtheria_corrob, vartype = "ci", na.rm = TRUE))

n_corrob <- sum(df_ind_corrob$diphtheria_corrob, na.rm = TRUE)
died_corrob <- sum(df_ind_corrob$diphtheria_corrob & df_ind_corrob$died, na.rm = TRUE)
n_corrob
died_corrob / n_corrob

svyglm(diphtheria_corrob ~ vacc, family = quasibinomial, design = df_srvy_corrob) |>
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) |>
  filter(term == "vaccTRUE")

# ============================================================
# 3. Factors associated with diphtheria (Table 2)
# ============================================================

run_svy_logistic <- function(outcome, vars_uni, vars_multi, design_object) {
  uni_results <- map_dfr(vars_uni, \(v) {
    mod <- svyglm(as.formula(paste(outcome, "~", v)), family = quasibinomial(), design = design_object)
    broom::tidy(mod, exponentiate = TRUE, conf.int = TRUE) |>
      filter(term != "(Intercept)") |>
      mutate(variable = v)
  }) |>
    select(variable, term, estimate, conf.low, conf.high, p.value) |>
    rename(uni_est = estimate, uni_low = conf.low, uni_high = conf.high, uni_p = p.value)

  multi_mod <- svyglm(as.formula(paste(outcome, "~", paste(vars_multi, collapse = " + "))), family = quasibinomial(), design = design_object)
  multi_results <-
    broom::tidy(multi_mod, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term != "(Intercept)") |>
    mutate(
      variable = case_when(
        grepl("age_grp_diph", term) ~ "age_grp_diph",
        grepl("sex", term) ~ "sex",
        grepl("vacc", term) ~ "vacc",
        grepl("ses_quintile", term) ~ "as.numeric(ses_quintile)",
        grepl("^lga", term) ~ "lga",
        TRUE ~ term
      )
    ) |>
    select(variable, term, multi_est = estimate, multi_low = conf.low, multi_high = conf.high, multi_p = p.value)

  left_join(uni_results, multi_results, by = c("variable", "term"))
}

# LGA is the reference-coded model of interest here (ref "16+" for age),
# so releveled before fitting
df_srvy_ind$variables$age_grp_diph <- df_srvy_ind$variables$age_grp_diph |> factor(ordered = FALSE) |> relevel(ref = "16+")

vars_table2 <- c("age_grp_diph", "sex", "vacc", "as.numeric(ses_quintile)", "hh_size", "lga")

table_2 <- run_svy_logistic(outcome = "diphtheria", vars_uni = vars_table2, vars_multi = vars_table2, design_object = df_srvy_ind)
print(table_2, n = Inf)
write_csv(table_2, file.path(tab_dir, "table2_diphtheria_risk_factors.csv"))

# ============================================================
# 4. Mortality
# ============================================================

calculate_mr <- function(data, key, value) {
  data |>
    filter(!is.na(cluster), !is.na(index_hh)) |>
    as_survey(ids = c(cluster, index_hh), strata = lga, weights = weight, nest = TRUE) |>
    summarise(mr = survey_ratio(numerator = {{ key }} %in% {{ value }} * 10000, denominator = pt, vartype = "ci", deff = TRUE, na.rm = TRUE))
}

# crude mortality rate, per 10,000 person-days
calculate_mr(df_ind, died, TRUE)

# under-5 mortality rate
df_ind |> filter(age_grp %in% "0-4") |> calculate_mr(died, TRUE)

# leading causes of death, survey-weighted proportion of all deaths
df_srvy_ind |>
  filter(died == TRUE) |>
  group_by(died_cause) |>
  summarise(pct = survey_mean(vartype = "ci") * 100) |>
  arrange(desc(pct))

# overall CFR among diphtheria cases
df_srvy_ind |> filter(diphtheria == TRUE) |> summarise(cfr = survey_mean(died * 100, vartype = "ci", na.rm = TRUE))

# CFR among children under 5
df_srvy_ind |>
  group_by(age_cat) |>
  filter(diphtheria == TRUE) |>
  summarise(cfr = survey_mean(died * 100, vartype = "ci", na.rm = TRUE)) |>
  filter(age_cat == "0-4")

# timing of diphtheria deaths relative to symptom onset
df_ind |>
  filter(diphtheria == TRUE, died == TRUE, !is.na(date_died), !is.na(date_onset)) |>
  mutate(days_onset_to_death = as.numeric(date_died - date_onset)) |>
  summarise(n = n(), median_days = median(days_onset_to_death), q1 = quantile(days_onset_to_death, .25), q3 = quantile(days_onset_to_death, .75))

# ============================================================
# 5. Health-seeking behaviour and treatment delay
# ============================================================

df_srvy_ind |> filter(diphtheria == TRUE) |> summarise(sought_care = survey_mean(tx_seeking == "Yes", vartype = "ci", na.rm = TRUE))

df_ind |>
  filter(diphtheria == TRUE, !is.na(delay_tx)) |>
  summarise(n = n(), median_days = median(delay_tx), q1 = quantile(delay_tx, .25), q3 = quantile(delay_tx, .75))

# reasons for not seeking care, among cases who did not seek care
df_ind |>
  filter(diphtheria == TRUE, tx_seeking == "No") |>
  summarise(
    low_severity = sum(tx_hf_no_3 == 1, na.rm = TRUE),
    facility_distance = sum(tx_hf_no_4 == 1, na.rm = TRUE),
    traditional_healer = sum(tx_hf_no_5 == 1, na.rm = TRUE),
    dont_know = sum(tx_hf_no_99 == 1, na.rm = TRUE),
    n = n()
  )

# absolute risk difference per additional day of delay (linear-probability
# model, restricted to cases with a recorded delay)
fit_delay_lp <- svyglm(died ~ delay_tx, family = gaussian, design = df_srvy_ind |> filter(diphtheria == TRUE))
broom::tidy(fit_delay_lp, conf.int = TRUE) |> filter(term == "delay_tx") |> mutate(across(c(estimate, conf.low, conf.high), \(x) x * 100))

# relative risk, dichotomised at >=4 days
calculate_rate_ratio <- function(survey_design, numerator, denominator, group_var, grp_ref, grp_compare, scale_factor = 10000) {
  numerator <- rlang::enquo(numerator)
  denominator <- rlang::enquo(denominator)
  group_var <- rlang::enquo(group_var)
  survey_design$variables <- survey_design$variables |> mutate(group = !!group_var)
  rates <-
    survey_design |>
    group_by(group) |>
    summarise(mr = survey_ratio(numerator = (!!numerator) * scale_factor, denominator = !!denominator, vartype = c("se", "ci"), na.rm = TRUE)) |>
    ungroup()
  rate1 <- rates |> filter(group == grp_compare) |> pull(mr)
  se1 <- rates |> filter(group == grp_compare) |> pull(mr_se)
  rate2 <- rates |> filter(group == grp_ref) |> pull(mr)
  se2 <- rates |> filter(group == grp_ref) |> pull(mr_se)
  tibble(
    rate_ratio = rate1 / rate2,
    ci_low = exp(log(rate1 / rate2) - 1.96 * sqrt((se1 / rate1)^2 + (se2 / rate2)^2)),
    ci_high = exp(log(rate1 / rate2) + 1.96 * sqrt((se1 / rate1)^2 + (se2 / rate2)^2))
  )
}

calculate_rate_ratio(df_srvy_ind |> filter(diphtheria == TRUE), died, pt, delay_tx_binary, grp_ref = FALSE, grp_compare = TRUE)

# place of death, any cause (too few diphtheria-specific deaths for a
# reliable estimate on their own)
map(levels(df_ind$died_location), \(lvl) {
  df_srvy_ind |>
    summarise(pct = survey_mean(died_location == lvl, vartype = "ci", na.rm = TRUE) * 100) |>
    mutate(died_location = lvl, .before = 1)
}) |>
  list_rbind()

# ============================================================
# 6. Vaccination and mortality
# ============================================================

# mortality rate ratio, unvaccinated as reference
calculate_rate_ratio(df_srvy_ind, died, pt, vacc, grp_ref = FALSE, grp_compare = TRUE)

# reciprocal (vaccinated as reference), for the E-value, which requires an
# RR >= 1
mrr_vacc <- calculate_rate_ratio(df_srvy_ind, died, pt, vacc, grp_ref = TRUE, grp_compare = FALSE)
mrr_vacc
evalues.RR(est = mrr_vacc$rate_ratio, lo = mrr_vacc$ci_low, hi = mrr_vacc$ci_high)

# crude vaccine effectiveness = 1 - mortality rate ratio (vaccinated vs
# unvaccinated)
mrr_vs_unvacc <- calculate_rate_ratio(df_srvy_ind, died, pt, vacc, grp_ref = FALSE, grp_compare = TRUE)
(1 - mrr_vs_unvacc$rate_ratio) * 100
(1 - mrr_vs_unvacc$ci_high) * 100
(1 - mrr_vs_unvacc$ci_low) * 100

# CFR by vaccination status among cases
df_srvy_ind |> filter(diphtheria == TRUE) |> group_by(vacc) |> summarise(cfr = survey_mean(died * 100, vartype = "ci", na.rm = TRUE))
svyglm(died ~ vacc, family = quasibinomial, design = df_srvy_ind |> filter(diphtheria == TRUE)) |> broom::tidy()

# ============================================================
# 7. Vaccination and diphtheria: sensitivity analyses (Table 3)
# ============================================================
# The comparator for every restricted analysis below is explicit
# self-report-unvaccinated (vacc == FALSE), not a derived binary's own
# FALSE category - e.g. vacc_epi_binary == FALSE means "reports
# vaccinated but no card evidence", not "unvaccinated".

extract_svyglm_or <- function(fit, term) {
  coef_v <- coef(fit)
  se_v <- sqrt(diag(vcov(fit)))
  est <- coef_v[[term]]
  se <- se_v[[term]]
  z <- est / se
  tibble(or = exp(est), lower = exp(est - 1.96 * se), upper = exp(est + 1.96 * se), p = 2 * pnorm(-abs(z)))
}

fit_t3_full <- svyglm(diphtheria ~ vacc, family = quasibinomial, design = df_srvy_ind)
res_t3_full <- extract_svyglm_or(fit_t3_full, "vaccTRUE") |> mutate(analysis = "Full sample (all respondents)")

d_t3_card <- df_ind |> filter(vacc_epi_binary == TRUE | vacc == FALSE, !is.na(diphtheria)) |> mutate(card_vs_unvacc = vacc_epi_binary == TRUE)
svy_t3_card <- svydesign(ids = ~ cluster + index_hh, strata = ~lga, weights = ~weight, data = d_t3_card, nest = TRUE)
fit_t3_card <- svyglm(diphtheria ~ card_vs_unvacc, family = quasibinomial, design = svy_t3_card)
res_t3_card <- extract_svyglm_or(fit_t3_card, "card_vs_unvaccTRUE") |> mutate(analysis = "Card-confirmed vaccination only")

d_t3_routine <- df_ind |>
  mutate(routine_only = vacc_epi_binary == TRUE & (vacc_camp_binary == FALSE | is.na(vacc_camp_binary))) |>
  filter(routine_only == TRUE | vacc == FALSE, !is.na(diphtheria))
svy_t3_routine <- svydesign(ids = ~ cluster + index_hh, strata = ~lga, weights = ~weight, data = d_t3_routine, nest = TRUE)
fit_t3_routine <- svyglm(diphtheria ~ routine_only, family = quasibinomial, design = svy_t3_routine)
res_t3_routine <- extract_svyglm_or(fit_t3_routine, "routine_onlyTRUE") |> mutate(analysis = "Routine vaccination only vs unvaccinated")

d_t3_camp <- df_ind |> filter(vacc_camp_binary == TRUE | vacc == FALSE, !is.na(diphtheria)) |> mutate(camp_vs_unvacc = vacc_camp_binary == TRUE)
svy_t3_camp <- svydesign(ids = ~ cluster + index_hh, strata = ~lga, weights = ~weight, data = d_t3_camp, nest = TRUE)
fit_t3_camp <- svyglm(diphtheria ~ camp_vs_unvacc, family = quasibinomial, design = svy_t3_camp)
res_t3_camp <- extract_svyglm_or(fit_t3_camp, "camp_vs_unvaccTRUE") |> mutate(analysis = "Any campaign vaccination vs unvaccinated")

table_3 <- bind_rows(res_t3_full, res_t3_card, res_t3_routine, res_t3_camp) |> select(analysis, or, lower, upper, p)
print(table_3)
write_csv(table_3, file.path(tab_dir, "table3_vaccination_sensitivity.csv"))

# among vaccinated cases, how many reported vaccination at a care contact
# for this illness (evidence for the reverse-causality mechanism)
df_ind |> filter(diphtheria == TRUE, vacc == TRUE) |> count(vacc_card == "Was vaccinated in the DTC/OPD/Contact clinic")

# ============================================================
# 8. Vaccination coverage and campaign performance (Table 4)
# ============================================================

df_cc <- df_ind |> filter(age_eligible, vacc_source != "unknown", !is.na(vacc_epi_confirmed), !is.na(camp_confirmed))
nrow(df_cc)

table_4 <-
  df_cc |>
  mutate(
    epi_label = if_else(vacc_epi_confirmed, "EPI confirmed", "EPI not confirmed"),
    camp_label = if_else(camp_confirmed, "Campaign confirmed", "Campaign not confirmed")
  ) |>
  count(epi_label, camp_label) |>
  pivot_wider(names_from = camp_label, values_from = n, values_fill = 0L)
print(table_4)
write_csv(table_4, file.path(tab_dir, "table4_epi_vs_campaign.csv"))

svy_derived <- df_srvy_ind

# survey-weighted routine EPI / campaign coverage, among all age-eligible children
svy_derived |>
  filter(age_eligible) |>
  summarise(
    epi_coverage = survey_mean(vacc_epi_confirmed == TRUE, vartype = "ci", na.rm = TRUE),
    camp_coverage = survey_mean(camp_confirmed == TRUE, vartype = "ci", na.rm = TRUE)
  )

# evidence of vaccination from at least one source, among children with
# a resolved vaccination source (Metric 4)
svy_derived |>
  filter(age_eligible, vacc_source != "unknown") |>
  summarise(any_source = survey_mean(vacc_source %in% c("both", "epi_only", "camp_only"), vartype = "ci", na.rm = TRUE))

# no vaccination of any kind, among all age-eligible children with
# self-reported vaccination status observed (Metric 3; vacc is already
# logical, so it's used directly rather than re-deriving vacc_log)
svy_derived |>
  filter(age_eligible, !is.na(vacc)) |>
  summarise(none = survey_mean(vacc == FALSE, vartype = "ci", na.rm = TRUE))

# Metric 1: % of campaign recipients with prior card-confirmed EPI
svy_derived |> filter(age_eligible, camp_confirmed == TRUE) |> summarise(est = survey_mean(vacc_epi_confirmed == TRUE, vartype = "ci", na.rm = TRUE))

# Metric 2: % of EPI-unvaccinated children reached by the campaign
svy_derived |>
  filter(age_eligible, vacc_epi_confirmed == FALSE, !is.na(camp_confirmed)) |>
  summarise(est = survey_mean(camp_confirmed == TRUE, vartype = "ci", na.rm = TRUE))

# ============================================================
# 9. Figure 1: epidemic curve
# ============================================================

cols_zone <- c("Dawakin Tofa" = "#5a4636", "Fagge" = "#599656", "Nassarawa" = "steelblue", "Ungogo" = "#b66d6d")

df_epi <-
  df_ind |>
  filter(diphtheria == TRUE) |>
  mutate(date_onset = lubridate::floor_date(date_onset, "month")) |>
  count(lga, date_onset) |>
  tidyr::complete(date_onset = seq(min(date_onset), max(date_onset), by = "month"), lga, fill = list(n = 0)) |>
  arrange(date_onset) |>
  mutate(month_ym = format(date_onset, "%b %Y") |> fct_inorder())

df_epi |>
  ggplot(aes(month_ym, n, fill = lga)) +
  geom_col(width = 1, colour = "black") +
  facet_wrap(~lga, ncol = 1, strip.position = "right") +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x), expand = c(0, 0)) +
  scale_fill_manual(values = cols_zone, guide = "none") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  labs(x = "Month of disease onset", y = "Number of individuals") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey90", colour = NA), strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "figure1_epicurve.jpg"), width = 12, height = 8, dpi = 300)

# ============================================================
# 10. Spatial figures (Figure 2, S1 Figure) and cluster-level
# correlation tests
# ============================================================
# Cluster-level rates are empirical-Bayes shrunk toward the survey-wide
# value before kernel-smoothing, since raw small-cluster rates would
# otherwise produce spurious hotspots; the formal correlation tests
# below use the raw (unshrunk) cluster estimates instead, so shrinkage
# doesn't artificially inflate them.

df_ind_maps <-
  df_ind |>
  mutate(
    vacc_any_binary = case_when(
      vacc_source %in% c("both", "epi_only", "camp_only") ~ TRUE,
      vacc_source == "neither" ~ FALSE,
      .default = NA
    )
  )

shrink_binom <- function(dat, outcome) {
  d <- dat |> filter(!is.na({{ outcome }})) |> mutate(.outcome = as.numeric({{ outcome }}), .cluster_fct = factor(cluster))
  if (n_distinct(d$.cluster_fct) < 2 || n_distinct(d$.outcome) < 2) return(NULL)
  m <- tryCatch(mgcv::gam(.outcome ~ s(.cluster_fct, bs = "re"), data = d, family = binomial()), error = \(e) NULL)
  if (is.null(m)) return(NULL)
  d |>
    distinct(.cluster_fct) |>
    mutate(shrunk = 100 * predict(m, newdata = tibble(.cluster_fct = .cluster_fct), type = "response")) |>
    rename(cluster = .cluster_fct)
}

shrink_rate <- function(dat, outcome) {
  d <- dat |> filter(!is.na({{ outcome }}), !is.na(pt), pt > 0) |> mutate(.outcome = as.numeric({{ outcome }}), .cluster_fct = factor(cluster))
  if (n_distinct(d$.cluster_fct) < 2 || n_distinct(d$.outcome) < 2) return(NULL)
  m <- tryCatch(mgcv::gam(.outcome ~ s(.cluster_fct, bs = "re") + offset(log(pt)), data = d, family = poisson()), error = \(e) NULL)
  if (is.null(m)) return(NULL)
  d |>
    distinct(.cluster_fct) |>
    mutate(shrunk = 10000 * predict(m, newdata = tibble(.cluster_fct = .cluster_fct, pt = 1), type = "response")) |>
    rename(cluster = .cluster_fct)
}

tbl_cluster_core <-
  df_ind_maps |>
  filter(!is.na(cluster)) |>
  group_by(cluster) |>
  summarise(
    num_ppl = n(),
    cases = sum(diphtheria %in% TRUE, na.rm = TRUE),
    attack_rate = cases / num_ppl * 100,
    died_n = sum(died %in% TRUE, na.rm = TRUE),
    pt_tot = sum(pt, na.rm = TRUE),
    mr = died_n / pt_tot * 10000
  )

tbl_cluster_u5 <-
  df_ind_maps |>
  filter(!is.na(cluster), age_grp %in% "0-4") |>
  group_by(cluster) |>
  summarise(num_ppl_u5 = n(), died_u5_n = sum(died %in% TRUE, na.rm = TRUE), pt_u5_tot = sum(pt, na.rm = TRUE), u5mr = died_u5_n / pt_u5_tot * 10000)

tbl_cluster_vacc <-
  df_ind_maps |>
  filter(!is.na(cluster), age_eligible %in% TRUE) |>
  group_by(cluster) |>
  summarise(
    num_ppl_vacc = n(),
    vacc_epi = sum(vacc_epi_confirmed %in% TRUE, na.rm = TRUE) / sum(!is.na(vacc_epi_confirmed)) * 100,
    vacc_camp = sum(camp_confirmed %in% TRUE, na.rm = TRUE) / sum(!is.na(camp_confirmed)) * 100,
    vacc_any_n = sum(vacc_source %in% c("both", "epi_only", "camp_only")),
    vacc_any_d = sum(vacc_source %in% c("both", "epi_only", "camp_only", "neither")),
    vacc_any = vacc_any_n / vacc_any_d * 100
  ) |>
  select(-vacc_any_n, -vacc_any_d)

tbl_cluster <- tbl_cluster_core |> full_join(tbl_cluster_u5, by = "cluster") |> full_join(tbl_cluster_vacc, by = "cluster")

# unshrunk, for the Spearman correlations below
tbl_cluster_raw <- tbl_cluster

apply_shrinkage <- function(tbl, raw_col, shrunk_tbl) {
  if (is.null(shrunk_tbl)) return(tbl)
  tbl |> left_join(shrunk_tbl, by = "cluster") |> mutate("{raw_col}" := coalesce(shrunk, .data[[raw_col]])) |> select(-shrunk)
}

df_ind_vacc <- df_ind_maps |> filter(age_eligible %in% TRUE)
df_ind_u5 <- df_ind_maps |> filter(age_grp %in% "0-4")

tbl_cluster <-
  tbl_cluster |>
  apply_shrinkage("attack_rate", shrink_binom(df_ind_maps, diphtheria)) |>
  apply_shrinkage("mr", shrink_rate(df_ind_maps, died)) |>
  apply_shrinkage("u5mr", shrink_rate(df_ind_u5, died)) |>
  apply_shrinkage("vacc_epi", shrink_binom(df_ind_vacc, vacc_epi_confirmed)) |>
  apply_shrinkage("vacc_camp", shrink_binom(df_ind_vacc, camp_confirmed)) |>
  apply_shrinkage("vacc_any", shrink_binom(df_ind_vacc, vacc_any_binary))

tbl_cluster_gps <-
  df_ind_maps |>
  filter(!is.na(cluster), !is.na(lat), !is.na(lon)) |>
  group_by(cluster) |>
  slice(1) |>
  ungroup() |>
  select(cluster, lat, lon)

tbl_cluster <-
  tbl_cluster |>
  inner_join(tbl_cluster_gps, by = "cluster") |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(3857)

shp_boundary <- map_nga_ward_union |> st_transform(3857)
shp_lga <- map_nga_ward |> st_transform(3857)

# indicator: column in tbl_cluster; weight_col: cluster-N column used to
# weight the smoothing; high_is_bad: colour direction (TRUE = high value
# red/bad, e.g. attack rate; FALSE = high value green/good, e.g.
# vaccination coverage); scale_max lets two maps share a colour ceiling
# for a single collected legend.
plot_ppp_srvy <- function(dat, shp, indicator, weight_col, legend_name, high_is_bad = TRUE, scale_max = NULL) {
  poly_win <- as.owin(shp)
  filtered_dat <- dat |> filter(!is.na(.data[[indicator]]))
  ppp_xy <- st_coordinates(filtered_dat) |> as_tibble()

  inside <- spatstat.geom::inside.owin(ppp_xy$X, ppp_xy$Y, poly_win)
  ppp_xy <- ppp_xy[inside, ]
  filtered_dat <- filtered_dat[inside, ]

  sm_dat <-
    ppp(x = ppp_xy$X, y = ppp_xy$Y, marks = filtered_dat[[indicator]], window = poly_win) |>
    Smooth(weights = filtered_dat[[weight_col]]) |>
    stars::st_as_stars() |>
    st_as_sf() |>
    st_set_crs(3857)

  grid_points <-
    sm_dat |>
    mutate(lon = st_coordinates(st_centroid(geometry))[, 1], lat = st_coordinates(st_centroid(geometry))[, 2]) |>
    st_drop_geometry()

  own_max <- max(grid_points$v, na.rm = TRUE)
  max_value <- if (!is.null(scale_max)) scale_max else own_max

  plot <-
    ggplot() +
    geom_sf(data = sm_dat, aes(fill = v), colour = NA) +
    geom_sf(data = st_boundary(shp_lga), colour = "grey40", linewidth = 0.4) +
    geom_sf(data = st_boundary(shp), colour = "black", linewidth = 0.8) +
    geom_contour(data = grid_points, aes(x = lon, y = lat, z = v), colour = "grey20", linewidth = 0.3, linetype = "dashed") +
    metR::geom_text_contour(data = grid_points, aes(x = lon, y = lat, z = v), colour = "black", size = 3, skip = 0, stroke = 0.1) +
    scale_fill_viridis_c(option = "plasma", name = legend_name, limits = c(0, max_value)) +
    coord_sf(crs = st_crs(shp)) +
    theme_void()

  lst(plot, max_value = own_max)
}

maps_srvy <- list(
  attack_rate = plot_ppp_srvy(tbl_cluster, shp_boundary, "attack_rate", "num_ppl", "Diphtheria attack rate (%)")$plot,
  u5mr = plot_ppp_srvy(tbl_cluster, shp_boundary, "u5mr", "num_ppl_u5", "Under-5 mortality rate\n(per 10,000 person-days)")$plot,
  vacc_any = plot_ppp_srvy(tbl_cluster, shp_boundary, "vacc_any", "num_ppl_vacc", "Vaccination coverage - any\n(EPI and/or campaign) (%)", high_is_bad = FALSE)$plot
)

# cluster-level correlation tests, raw (unshrunk) cluster estimates
spearman_test <- function(x_col, y_col, label) {
  test <- cor.test(tbl_cluster_raw[[x_col]], tbl_cluster_raw[[y_col]], use = "complete.obs", method = "spearman")
  tibble(pair = label, rho = unname(test$estimate), p = test$p.value)
}

whole_pop_pairs <- tribble(
  ~x_col,      ~y_col,        ~label,
  "vacc_epi",  "vacc_camp",   "EPI vs campaign coverage",
  "vacc_any",  "mr",          "Vaccinated (any) vs all-cause mortality",
  "vacc_any",  "u5mr",        "Vaccinated (any) vs under-5 mortality",
  "vacc_any",  "attack_rate", "Vaccinated (any) vs diphtheria attack rate",
  "vacc_epi",  "attack_rate", "EPI coverage vs attack rate",
  "vacc_camp", "attack_rate", "Campaign coverage vs attack rate",
  "vacc_epi",  "u5mr",        "EPI coverage vs under-5 mortality",
  "vacc_camp", "u5mr",        "Campaign coverage vs under-5 mortality"
)

spatial_correlations <- pmap_dfr(whole_pop_pairs, spearman_test)
print(spatial_correlations, n = Inf)
write_csv(spatial_correlations, file.path(tab_dir, "spatial_correlations.csv"))

# age-restricted checks (campaign-eligible and under-5 populations), cited
# in Limitations as holding under any population restriction tried
tbl_cluster_eligible <-
  df_ind_maps |>
  filter(!is.na(cluster), age_eligible %in% TRUE) |>
  group_by(cluster) |>
  summarise(
    num_ppl_eligible = n(),
    attack_rate_eligible = sum(diphtheria %in% TRUE, na.rm = TRUE) / n() * 100,
    mr_eligible = sum(died %in% TRUE, na.rm = TRUE) / sum(pt, na.rm = TRUE) * 10000
  )

tbl_cluster_u5_ar <-
  df_ind_maps |>
  filter(!is.na(cluster), age_grp %in% "0-4") |>
  group_by(cluster) |>
  summarise(attack_rate_u5 = sum(diphtheria %in% TRUE, na.rm = TRUE) / n() * 100)

tbl_cluster_raw_restricted <-
  tbl_cluster_raw |>
  select(cluster, vacc_epi, vacc_camp, vacc_any) |>
  left_join(tbl_cluster_eligible, by = "cluster") |>
  left_join(tbl_cluster_u5, by = "cluster") |>
  left_join(tbl_cluster_u5_ar, by = "cluster")

spearman_restricted <- function(x_col, y_col, label) {
  test <- cor.test(tbl_cluster_raw_restricted[[x_col]], tbl_cluster_raw_restricted[[y_col]], use = "complete.obs", method = "spearman")
  tibble(pair = label, rho = unname(test$estimate), p = test$p.value)
}

restricted_pairs <- tribble(
  ~x_col,      ~y_col,                 ~label,
  "vacc_epi",  "attack_rate_eligible", "EPI coverage vs attack rate (campaign-eligible)",
  "vacc_camp", "attack_rate_eligible", "Campaign coverage vs attack rate (campaign-eligible)",
  "vacc_any",  "attack_rate_eligible", "Vaccinated (any) vs attack rate (campaign-eligible)",
  "vacc_epi",  "mr_eligible",          "EPI coverage vs mortality (campaign-eligible)",
  "vacc_camp", "mr_eligible",          "Campaign coverage vs mortality (campaign-eligible)",
  "vacc_any",  "mr_eligible",          "Vaccinated (any) vs mortality (campaign-eligible)",
  "vacc_epi",  "attack_rate_u5",       "EPI coverage vs attack rate (under-5)",
  "vacc_camp", "attack_rate_u5",       "Campaign coverage vs attack rate (under-5)",
  "vacc_any",  "attack_rate_u5",       "Vaccinated (any) vs attack rate (under-5)",
  "vacc_epi",  "u5mr",                 "EPI coverage vs mortality (under-5)",
  "vacc_camp", "u5mr",                 "Campaign coverage vs mortality (under-5)",
  "vacc_any",  "u5mr",                 "Vaccinated (any) vs mortality (under-5)"
)

spatial_correlations_restricted <- pmap_dfr(restricted_pairs, spearman_restricted)
print(spatial_correlations_restricted, n = Inf)
write_csv(spatial_correlations_restricted, file.path(tab_dir, "spatial_correlations_age_restricted.csv"))

# Figure 2: EPI vs campaign coverage, shared colour scale, one collected
# legend below the panels
epi_raw <- plot_ppp_srvy(tbl_cluster, shp_boundary, "vacc_epi", "num_ppl_vacc", "Vaccination coverage (%)", high_is_bad = FALSE)
camp_raw <- plot_ppp_srvy(tbl_cluster, shp_boundary, "vacc_camp", "num_ppl_vacc", "Vaccination coverage (%)", high_is_bad = FALSE)
epi_camp_shared_max <- max(epi_raw$max_value, camp_raw$max_value)

epi_shared <- plot_ppp_srvy(tbl_cluster, shp_boundary, "vacc_epi", "num_ppl_vacc", "Vaccination coverage (%)", high_is_bad = FALSE, scale_max = epi_camp_shared_max)$plot
camp_shared <- plot_ppp_srvy(tbl_cluster, shp_boundary, "vacc_camp", "num_ppl_vacc", "Vaccination coverage (%)", high_is_bad = FALSE, scale_max = epi_camp_shared_max)$plot

figure_2 <-
  (epi_shared + labs(title = "a. Routine EPI coverage")) +
  (camp_shared + labs(title = "b. Campaign coverage")) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "figure2_vaccination_coverage.jpg"), figure_2, height = 6.5, width = 12)

# S1 Figure: vaccinated (any), attack rate, under-5 mortality - different
# units, so each panel keeps its own colour scale
s1_figure <-
  (maps_srvy$vacc_any + labs(title = "a. Vaccinated (any)")) +
  (maps_srvy$attack_rate + labs(title = "b. Diphtheria attack rate")) +
  (maps_srvy$u5mr + labs(title = "c. Under-5 mortality")) +
  plot_layout(ncol = 3)

ggsave(file.path(fig_dir, "S1_figure_vaccination_ar_u5mr.jpg"), s1_figure, height = 6, width = 20)
