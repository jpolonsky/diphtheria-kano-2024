# Community burden of diphtheria during the 2023-24 epidemic in Kano
# State, Nigeria: population-based household survey.
#
# Reproduces every number, table and figure reported in the manuscript,
# in the order they appear (Participants -> Community attack rate and
# comparison with facility notifications -> Factors associated with
# diphtheria -> Mortality -> Health-seeking and treatment delay ->
# Vaccination and mortality -> Vaccination sensitivity analyses ->
# Vaccination coverage and campaign performance -> Figure 1 -> spatial
# figures and correlations). Input is the frozen, de-identified survey
# dataset; see README for how to request it. The facility comparator is a
# fixed external count (MSF line list, aggregated in S2 Table), and a few
# descriptive items rest on variables not in the released file (symptom
# checklist, vaccination dates, household roster changes); those are not
# recomputed.

pacman::p_load(survey, srvyr, scales, sf, spatstat, stars, mgcv, metR, patchwork, tidyverse)

fig_dir <- "output/figures"
tab_dir <- "output/tables"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 0. Data
# ============================================================
# df_ind: one row per surveyed individual (n=7,998), unweighted - used for
#   raw counts (Table 1) and as the basis for the survey design below
# df_hh (data/df_hh_srvy.rds): one row per approached household (n=1,083) - used for the
#   response-rate figure
# map_boundaries: LGA and study-area boundary polygons, used only by the
#   spatial figures at the end of this script

df_ind <- read_rds("data/df_ind.rds") |>
  mutate(
    # restricted case definition: cases also reported as diagnosed by a health professional
    diphtheria_corrob = diphtheria == TRUE & diagnosis == "Yes",
    card_chr = as.character(vacc_card)
  )
df_hh <- read_rds("data/df_hh_srvy.rds")
maps <- read_rds("data/map_boundaries.rds")
map_nga_ward <- maps$map_nga_ward
map_nga_ward_union <- maps$map_nga_ward_union

options(survey.lonely.psu = "adjust")

# Design strata: eight clusters have households in two LGAs, but each cluster
# was sampled once, from one LGA, so each cluster is assigned to the LGA
# holding most of its households (lga_samp). The 78 clusters are the primary
# sampling units; weights stay on each household's own LGA.
lga_samp_lookup <-
  df_ind |>
  filter(!is.na(cluster)) |>
  distinct(cluster, index_hh, lga) |>
  count(cluster, lga) |>
  slice_max(n, n = 1, by = cluster, with_ties = FALSE) |>
  select(cluster, lga_samp = lga)

df_ind <- df_ind |> left_join(lga_samp_lookup, by = "cluster")
stopifnot(
  n_distinct(df_ind$cluster[!is.na(df_ind$cluster)]) == 78,
  n_distinct(paste(df_ind$lga_samp, df_ind$cluster)[!is.na(df_ind$cluster)]) == 78
)

# two-stage design (cluster, then household within cluster), stratified by
# lga_samp; no finite population correction
as_kano_design <- function(data) {
  data |>
    filter(!is.na(cluster), !is.na(index_hh)) |>
    as_survey(ids = c(cluster, index_hh), strata = lga_samp, weights = weight, nest = TRUE)
}

df_srvy_ind <- as_kano_design(df_ind)
design_df <- degf(df_srvy_ind)

# ============================================================
# 1. Participants (Table 1)
# ============================================================

nrow(df_hh)
sum(df_hh$consent == "Yes")
mean(df_hh$consent == "Yes")

# Diphtheria status is missing for some respondents. They are left out of the
# diphtheria analyses but kept in the all-cause mortality analyses.
flow_counts <- tibble(
  item = c(
    "Surveyed individuals", "Diphtheria status recorded", "Diphtheria cases", "No diphtheria",
    "Diphtheria status not recorded", "Deaths, all cause", "Deaths with diphtheria status recorded"
  ),
  n = c(
    nrow(df_ind), sum(!is.na(df_ind$diphtheria)), sum(df_ind$diphtheria, na.rm = TRUE), sum(!df_ind$diphtheria, na.rm = TRUE),
    sum(is.na(df_ind$diphtheria)), sum(df_ind$died, na.rm = TRUE), sum(df_ind$died & !is.na(df_ind$diphtheria), na.rm = TRUE)
  )
)
print(flow_counts)
write_csv(flow_counts, file.path(tab_dir, "flow_counts.csv"))

# enrolled households and clusters (the primary sampling units), clusters per
# design stratum, and the design degrees of freedom
df_ind |> summarise(n_households = n_distinct(index_hh), n_clusters = n_distinct(cluster))
lga_samp_lookup |> count(lga_samp, name = "n_clusters")
design_df

df_ind |>
  summarise(
    pct_female = mean(sex == "Female") * 100, median_age = median(age), q1 = quantile(age, .25), q3 = quantile(age, .75),
    pct_age_le15 = mean(age <= 15) * 100
  )

# Table 1: counts are unweighted, column percentages survey-weighted, and
# missing values are excluded from each variable
table1_vars <- c("lga", "age_grp_diph", "sex", "vacc", "died")
table_1 <-
  table1_vars |>
  map_dfr(\(v) {
    df_srvy_ind |>
      filter(!is.na(diphtheria), !is.na(.data[[v]])) |>
      mutate(level = as.character(.data[[v]])) |>
      group_by(diphtheria, level) |>
      summarise(n = unweighted(n()), pct = survey_prop() * 100, .groups = "drop") |>
      mutate(variable = v)
  })
print(table_1, n = Inf)

# Rao-Scott F tests of each characteristic against diphtheria status
table_1_p <- map_dfr(table1_vars, \(v) {
  d <- df_srvy_ind |> filter(!is.na(diphtheria), !is.na(.data[[v]]))
  tibble(variable = v, p = svychisq(as.formula(paste0("~", v, "+diphtheria")), design = d)$p.value)
})
print(table_1_p)

# vaccination status missing, by diphtheria status
df_ind |> filter(!is.na(diphtheria), is.na(vacc)) |> count(diphtheria)

write_csv(table_1, file.path(tab_dir, "table1_weighted.csv"))
write_csv(table_1_p, file.path(tab_dir, "table1_pvalues.csv"))

# ============================================================
# 2. Community attack rate and comparison with facility surveillance
# ============================================================

df_srvy_ind |>
  summarise(attack_rate = survey_mean(diphtheria * 100, vartype = "ci", deff = TRUE, na.rm = TRUE))

# incidence rate, reported cases per 1,000 person-years (ratio estimator over
# person-time pt in days; missing diphtheria status dropped from both terms),
# t CI on the design df
inc_ratio <- svyratio(~ as.numeric(diphtheria), ~ I(pt / 365.25 / 1000), design = df_srvy_ind, na.rm = TRUE)
c(ir = coef(inc_ratio)[[1]], confint(inc_ratio, df = design_df)[1, ])

df_srvy_ind |>
  group_by(age_grp_diph) |>
  summarise(attack_rate = survey_mean(diphtheria * 100, vartype = "ci", na.rm = TRUE))

df_srvy_ind |>
  group_by(sex) |>
  summarise(attack_rate = survey_mean(diphtheria * 100, vartype = "ci", na.rm = TRUE))

svyttest(diphtheria ~ sex, design = df_srvy_ind)

# household secondary attack rate: among households with >=1 case, how
# many additional ("secondary") cases among the remaining household
# members
hh_summary <-
  df_ind |>
  group_by(index_hh, cluster, lga, lga_samp, weight) |>
  summarise(hh_size = n(), num_cases = sum(diphtheria, na.rm = TRUE), .groups = "drop") |>
  filter(num_cases > 0) |>
  mutate(secondary_cases = pmax(num_cases - 1, 0), exposed = pmax(hh_size - 1, 0), sar = secondary_cases / exposed)

nrow(hh_summary)
sum(hh_summary$num_cases >= 2)

hh_design <- as_kano_design(hh_summary)
hh_design |> summarise(SAR = survey_mean(sar * 100, vartype = "ci", deff = TRUE, na.rm = TRUE))

# Facility comparator for the same four LGAs: 9,864 residents notified at the
# three MSF-supported sites (IDH, MMSH, WSGGH), dated by symptom onset, 1 January
# 2023 to 25 October 2024 (the survey recall period), over a population of
# 2,837,111 (WorldPop and National Population Commission estimates), or 347.7 per
# 100,000. The counts come from the MSF line list, which is not public
# (aggregated counts are in S2 Table), so this is a fixed external rate, not
# survey data.
facility_ar <- 9864 / 2837111

# Survey attack rate and its ratio to the facility rate, under the primary case
# definition and the restricted one. The CI of the ratio reflects survey
# sampling uncertainty only, since the facility rate is treated as fixed.
ar_vs_facility <-
  bind_rows(
    df_srvy_ind |>
      summarise(ar = survey_mean(diphtheria, vartype = "ci", na.rm = TRUE)) |>
      mutate(definition = "Primary (suspected diphtheria)"),
    df_srvy_ind |>
      summarise(ar = survey_mean(diphtheria_corrob, vartype = "ci", na.rm = TRUE)) |>
      mutate(definition = "Restricted (diagnosed by a health professional)")
  ) |>
  mutate(ratio = ar / facility_ar, ratio_low = ar_low / facility_ar, ratio_upp = ar_upp / facility_ar) |>
  select(definition, ar, ar_low, ar_upp, ratio, ratio_low, ratio_upp)
print(ar_vs_facility)
write_csv(ar_vs_facility, file.path(tab_dir, "ar_vs_facility.csv"))

# cases meeting the restricted definition
df_ind |>
  summarise(n_cases = sum(diphtheria, na.rm = TRUE), n_restricted = sum(diphtheria_corrob, na.rm = TRUE), pct = n_restricted / n_cases * 100)

# vaccination and diphtheria under the restricted definition
svyglm(diphtheria_corrob ~ vacc, family = quasibinomial, design = df_srvy_ind) |>
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
        str_detect(term, "age_grp_diph") ~ "age_grp_diph",
        str_detect(term, "sex") ~ "sex",
        str_detect(term, "vacc") ~ "vacc",
        str_detect(term, "ses_quintile") ~ "as.numeric(ses_quintile)",
        str_detect(term, "^lga") ~ "lga",
        .default = term
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

# complete-case samples: multivariable model, and the crude vaccination model
length(svyglm(reformulate(vars_table2, "diphtheria"), family = quasibinomial, design = df_srvy_ind)$y)
length(svyglm(diphtheria ~ vacc, family = quasibinomial, design = df_srvy_ind)$y)

# ============================================================
# 4. Mortality
# ============================================================

calculate_mr <- function(data, key, value) {
  data |>
    as_kano_design() |>
    summarise(mr = survey_ratio(numerator = {{ key }} %in% {{ value }} * 10000, denominator = pt, vartype = "ci", deff = TRUE, na.rm = TRUE))
}

# crude mortality rate, per 10,000 person-days
calculate_mr(df_ind, died, TRUE)

# under-5 mortality rate
df_ind |> filter(age_grp %in% "0-4") |> calculate_mr(died, TRUE)

# diphtheria-specific mortality rate (deaths attributed to diphtheria)
calculate_mr(df_ind, died_diph_binary, TRUE)

# leading causes of death, survey-weighted proportion of all deaths
df_srvy_ind |>
  filter(died == TRUE) |>
  group_by(died_cause) |>
  summarise(pct = survey_mean(vartype = "ci") * 100) |>
  arrange(desc(pct))

# Case-fatality ratio (CFR) with survey CIs on the logit scale, which stay
# within 0-100% when deaths are few
cfr_ci <- function(des) {
  ci <- svyciprop(~died, des, method = "logit")
  tibble(
    n_cases = nrow(des$variables), n_deaths = sum(des$variables$died, na.rm = TRUE),
    cfr = 100 * as.numeric(coef(ci)), lower = 100 * confint(ci)[1, 1], upper = 100 * confint(ci)[1, 2]
  )
}
des_cases <- subset(df_srvy_ind, diphtheria == TRUE)
cfr_table <- bind_rows(
  cfr_ci(des_cases) |> mutate(group = "Overall"),
  cfr_ci(subset(des_cases, age_cat == "0-4")) |> mutate(group = "Under 5 years"),
  cfr_ci(subset(des_cases, age_cat == "5+")) |> mutate(group = "5 years and over"),
  cfr_ci(subset(df_srvy_ind, diphtheria_corrob %in% TRUE)) |> mutate(group = "Restricted case definition")
) |>
  select(group, n_cases, n_deaths, cfr, lower, upper)
print(cfr_table)
write_csv(cfr_table, file.path(tab_dir, "cfr_logit_ci.csv"))

# Timing of diphtheria deaths relative to symptom onset. Some records are
# internally inconsistent: onset after the date of death, or a reported
# care-seeking delay longer than the interval from onset to death.
deaths_diph <-
  df_ind |>
  filter(diphtheria == TRUE, died == TRUE, !is.na(date_died), !is.na(date_onset)) |>
  mutate(
    days_onset_to_death = as.numeric(date_died - date_onset),
    onset_after_death = days_onset_to_death < 0,
    delay_exceeds_interval = delay_tx > days_onset_to_death,
    inconsistent = onset_after_death %in% TRUE | delay_exceeds_interval %in% TRUE
  )

deaths_diph |>
  summarise(
    n = n(), median_days = median(days_onset_to_death), q1 = quantile(days_onset_to_death, .25), q3 = quantile(days_onset_to_death, .75),
    n_onset_after_death = sum(onset_after_death), n_delay_exceeds_interval = sum(delay_exceeds_interval & !onset_after_death),
    n_inconsistent = sum(inconsistent), median_days_consistent = median(days_onset_to_death[!inconsistent])
  )

# individuals with an inconsistent record, used for the sensitivity analysis in Section 5
inconsistent_deaths <- deaths_diph |> filter(inconsistent) |> transmute(id = paste(index_hh, index_ind)) |> pull(id)

# place of death, any cause (too few diphtheria-specific deaths for a
# reliable estimate on their own), among deaths with a recorded place
df_srvy_ind |>
  summarise(
    died_home = survey_mean(died_location == "At home", vartype = "ci", na.rm = TRUE),
    died_facility = survey_mean(died_location == "At health facility", vartype = "ci", na.rm = TRUE),
    died_transit = survey_mean(died_location == "In transit to health facility", vartype = "ci", na.rm = TRUE)
  ) |>
  print(width = Inf)

# care sought in the 2 weeks before death, all-cause deaths
deaths_all <- df_ind |> filter(died %in% TRUE)
deaths_all |> count(tx_seeking_hsb)
deaths_all |> count(died_location, tx_seeking_hsb)

# diphtheria deaths: place of death, care sought and where treated
deaths_all |> filter(diphtheria %in% TRUE) |> count(died_location, tx_seeking, tx_location)

# ============================================================
# 5. Health-seeking behaviour and treatment delay
# ============================================================

df_ind |> filter(diphtheria == TRUE) |> count(tx_seeking)
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

# risk of death among cases by treatment delay (>=4 vs <4 days)
des_delay <- subset(df_srvy_ind, diphtheria %in% TRUE & !is.na(delay_tx_binary))
delay_cfr <- bind_rows(
  cfr_ci(subset(des_delay, delay_tx_binary == FALSE)) |> mutate(delay = "<4 days"),
  cfr_ci(subset(des_delay, delay_tx_binary == TRUE)) |> mutate(delay = ">=4 days")
) |>
  select(delay, n_cases, n_deaths, cfr, lower, upper)
print(delay_cfr)

# risk ratio for >=4 vs <4 days (survey-weighted quasi-Poisson, log link)
delay_risk_ratio <- function(des) {
  fit <- svyglm(died ~ delay_tx_binary, family = quasipoisson(link = "log"), design = des)
  tibble(
    risk_ratio = exp(coef(fit)[["delay_tx_binaryTRUE"]]),
    ci_low = exp(confint(fit)["delay_tx_binaryTRUE", 1]),
    ci_high = exp(confint(fit)["delay_tx_binaryTRUE", 2]),
    p = summary(fit)$coefficients["delay_tx_binaryTRUE", "Pr(>|t|)"]
  )
}

# also without the deaths with inconsistent dates or delays
bind_rows(
  delay_risk_ratio(des_delay) |> mutate(records = "All cases with a recorded delay"),
  delay_risk_ratio(subset(des_delay, !(paste(index_hh, index_ind) %in% inconsistent_deaths))) |>
    mutate(records = "Excluding deaths with inconsistent dates or delays")
)

# ============================================================
# 6. Vaccination and mortality
# ============================================================
# Exploratory: the timing of vaccination relative to illness or death is
# unknown. Individuals with missing vaccination status are excluded.

# Ratio of survey-weighted rates (events per person-day) between two groups,
# from a survey-weighted quasi-Poisson model with a log person-time offset,
# fitted on the two groups together. Both rates come from the same clusters,
# so their covariance is carried by the design-based variance instead of
# being assumed zero. The point estimate equals the ratio of the two weighted
# rates; CI and p-value are t-based on the model's residual df. Rates are
# returned per scale_factor person-days.
calculate_rate_ratio <- function(survey_design, numerator, denominator, group_var, grp_ref, grp_compare, scale_factor = 10000) {
  survey_design$variables <-
    survey_design$variables |>
    mutate(group = {{ group_var }}, .rr_num = as.numeric({{ numerator }}), .rr_den = {{ denominator }})
  des <- subset(survey_design, group %in% c(grp_ref, grp_compare) & !is.na(.rr_num) & !is.na(.rr_den) & .rr_den > 0)
  des$variables$.rr_cmp <- as.numeric(des$variables$group == grp_compare)
  fit <- svyglm(.rr_num ~ .rr_cmp + offset(log(.rr_den)), family = quasipoisson(link = "log"), design = des)
  est <- coef(fit)[[".rr_cmp"]]
  se <- sqrt(vcov(fit)[".rr_cmp", ".rr_cmp"])
  df_t <- fit$df.residual
  rates <- svyby(~ I(.rr_num * scale_factor), ~group, denominator = ~.rr_den, design = des, FUN = svyratio)
  tibble(
    rate_ratio = exp(est),
    ci_low = exp(est - qt(0.975, df_t) * se),
    ci_high = exp(est + qt(0.975, df_t) * se),
    p = 2 * pt(abs(est / se), df_t, lower.tail = FALSE),
    rate_compare = coef(rates)[[as.character(grp_compare)]],
    rate_ref = coef(rates)[[as.character(grp_ref)]]
  )
}

# mortality rate ratios, unvaccinated vs vaccinated (reference), whole surveyed
# population
rate_ratios <-
  bind_rows(
    calculate_rate_ratio(df_srvy_ind, died, pt, vacc, grp_ref = TRUE, grp_compare = FALSE) |> mutate(outcome = "All-cause"),
    calculate_rate_ratio(df_srvy_ind, died_diph_binary, pt, vacc, grp_ref = TRUE, grp_compare = FALSE) |> mutate(outcome = "Diphtheria-attributed")
  ) |>
  select(outcome, everything())
print(rate_ratios, width = Inf)
write_csv(rate_ratios, file.path(tab_dir, "mortality_rate_ratios.csv"))

# events and person-days by vaccination status (missing status shown for
# reference; those individuals are excluded from the ratios above)
df_ind |>
  group_by(vacc) |>
  summarise(
    n = n(), person_days = sum(pt, na.rm = TRUE),
    deaths_all_cause = sum(died, na.rm = TRUE),
    deaths_diphtheria = sum(died_diph_binary, na.rm = TRUE)
  )

# all-cause rate ratio adjusted for age group (quasi-Poisson, person-time offset)
fit_mrr_age <- svyglm(
  died ~ unvacc + age_cat + offset(log(pt)), family = quasipoisson,
  design = df_srvy_ind |> filter(!is.na(vacc), pt > 0) |> mutate(unvacc = !vacc)
)
exp(cbind(rr = coef(fit_mrr_age), confint(fit_mrr_age)))["unvaccTRUE", , drop = FALSE]

# CFR by vaccination status among cases, with a design-based test
des_cases_v <- subset(des_cases, !is.na(vacc))
cfr_by_vacc <- bind_rows(
  cfr_ci(subset(des_cases_v, vacc == FALSE)) |> mutate(vacc = "Unvaccinated"),
  cfr_ci(subset(des_cases_v, vacc == TRUE)) |> mutate(vacc = "Vaccinated")
) |>
  select(vacc, n_cases, n_deaths, cfr, lower, upper)
print(cfr_by_vacc)
svychisq(~ died + vacc, design = des_cases_v)
write_csv(cfr_by_vacc, file.path(tab_dir, "cfr_by_vaccination.csv"))

# ============================================================
# 7. Vaccination and diphtheria: sensitivity analyses (Table 3)
# ============================================================
# Crude survey-weighted ORs, each vaccinated group against individuals
# reported as never vaccinated. Row 1 covers all ages. Rows 2-4 use
# categories the questionnaire records, among individuals under 16 at
# interview (the campaign question was asked only below 16), and do not
# overlap: individuals vaccinated at a DTC, OPD or contact clinic (doses given
# on presentation under the outbreak-response protocol) form row 4 whatever
# their campaign answer and are left out of rows 2 and 3. Each row is a
# domain of the full design; CIs and p-values use each model's residual df,
# as in Table 2.

extract_svyglm_or <- function(fit, term) {
  coef_v <- coef(fit)
  se_v <- sqrt(diag(vcov(fit)))
  est <- coef_v[[term]]
  se <- se_v[[term]]
  # t distribution on the design degrees of freedom (as broom::tidy() in Table 2)
  df_t <- fit$df.residual
  z <- est / se
  tibble(or = exp(est), lower = exp(est - qt(0.975, df_t) * se), upper = exp(est + qt(0.975, df_t) * se), p = 2 * pt(-abs(z), df = df_t))
}

svy_t3 <-
  df_ind |>
  mutate(
    t3_group = case_when(
      vacc %in% FALSE ~ "never",
      vacc %in% TRUE & card_chr %in% "Was vaccinated in the DTC/OPD/Contact clinic" ~ "clinic",
      vacc %in% TRUE & vacc_camp_binary %in% TRUE ~ "campaign",
      vacc %in% TRUE & vacc_camp_binary %in% FALSE ~ "not_campaign",
      .default = NA
    ),
    t3_clinic = vacc %in% TRUE & card_chr %in% "Was vaccinated in the DTC/OPD/Contact clinic",
    t3_u16 = age < 16
  ) |>
  as_kano_design()

t3_row <- function(des, exposed, label) {
  d <- des |> mutate(.x = {{ exposed }}) |> filter(!is.na(.x), !is.na(diphtheria))
  fit <- svyglm(diphtheria ~ .x, family = quasibinomial, design = d)
  dat <- d$variables
  extract_svyglm_or(fit, ".xTRUE") |>
    mutate(
      analysis = label,
      n_exposed = sum(dat$.x), cases_exposed = sum(dat$diphtheria[dat$.x]),
      n_never = sum(!dat$.x), cases_never = sum(dat$diphtheria[!dat$.x])
    )
}

svy_t3_u16 <- svy_t3 |> filter(t3_u16 %in% TRUE)

table_3 <-
  bind_rows(
    t3_row(svy_t3, vacc, "Ever vaccinated, all ages"),
    t3_row(svy_t3_u16 |> filter(t3_group %in% c("campaign", "never")), t3_group == "campaign", "Campaign dose reported (no clinic dose), under 16"),
    t3_row(svy_t3_u16 |> filter(t3_group %in% c("not_campaign", "never")), t3_group == "not_campaign", "Vaccinated, no campaign or clinic dose reported, under 16"),
    t3_row(svy_t3_u16 |> filter(t3_clinic | vacc %in% FALSE), t3_clinic, "Vaccinated at a DTC, OPD or contact clinic, under 16")
  ) |>
  select(analysis, n_exposed, cases_exposed, n_never, cases_never, or, lower, upper, p)
print(table_3, width = Inf)
write_csv(table_3, file.path(tab_dir, "table3_vaccination_sensitivity.csv"))

# among vaccinated cases, how many reported vaccination at a care contact
# for this illness (evidence for the reverse-causality mechanism)
df_ind |> filter(diphtheria == TRUE, vacc == TRUE) |> count(vacc_card == "Was vaccinated in the DTC/OPD/Contact clinic")

# ============================================================
# 8. Vaccination coverage and campaign performance (Table 4, S3 Table)
# ============================================================
# Only what the questionnaire records is used. The card question records
# documentation (card seen, card lost, no card, or a dose given at a DTC, OPD
# or contact clinic), not where a dose came from, and the card date is the
# most recent vaccination of any kind, so no routine EPI status is derived
# and card dates classify nothing. Campaign status is the caregiver's answer
# alone, among children reported as vaccinated. On the first survey day the
# campaign question was also put to never-vaccinated children; those answers
# are not used, and never-vaccinated children stay never vaccinated.
#
# Campaign-eligible children: within the campaigns' target ages (6 weeks to
# 14 years) during at least one campaign phase, i.e. aged 6 weeks or more by
# the end of the last phase and 14 completed years or younger at the start of
# the first, and 15 completed years or younger at interview (age_eligible),
# because the campaign question was asked only below 16.

phases <- tribble(
  ~round, ~start,                ~end,
  "r1",   as.Date("2023-03-01"), as.Date("2023-04-30"),
  "r2",   as.Date("2023-08-01"), as.Date("2023-10-31"),
  "r3",   as.Date("2023-11-01"), as.Date("2024-01-31")
)
camp_age_min <- 6 / 52 # years
camp_age_max <- 14     # completed years, so exact age below camp_age_max + 1
q_age_max <- 15        # campaign question asked only below 16 at interview
camp_start <- min(phases$start)
camp_end <- max(phases$end)

hist_levels <- c(
  "Reported never vaccinated", "Card seen", "Card reported lost",
  "Reported vaccinated, no card", "Vaccinated at a DTC, OPD or contact clinic",
  "Vaccination status missing"
)

df_vacc <-
  df_ind |>
  mutate(
    # exact age at interview, taking the midpoint of the completed unit (`age`
    # is in years for records entered in months too)
    age_exact = if_else(age_unit == "Month", age + 0.5 / 12, age + 0.5),
    age_at_camp_end = age_exact - as.numeric(date_int - camp_end) / 365.25,
    age_at_camp_start = age_exact - as.numeric(date_int - camp_start) / 365.25,
    elig_camp = age_at_camp_end >= camp_age_min & age_at_camp_start < camp_age_max + 1 & age_eligible,
    # within the target ages during the campaigns but aged past 15 by the
    # interview: not classifiable, as the campaign question was not asked
    aged_out = age_at_camp_start < camp_age_max + 1 & age > q_age_max,
    history = factor(case_when(
      vacc %in% FALSE ~ "Reported never vaccinated",
      card_chr %in% "Yes" ~ "Card seen",
      card_chr %in% "Yes, but it is missing/lost" ~ "Card reported lost",
      card_chr %in% "No" ~ "Reported vaccinated, no card",
      card_chr %in% "Was vaccinated in the DTC/OPD/Contact clinic" ~ "Vaccinated at a DTC, OPD or contact clinic",
      .default = "Vaccination status missing"
    ), levels = hist_levels),
    camp_rep = if_else(vacc %in% TRUE, vacc_camp_binary, NA),
    camp_lbl = case_when(
      vacc %in% FALSE ~ "Not applicable (never vaccinated)",
      camp_rep %in% TRUE ~ "Campaign dose reported",
      camp_rep %in% FALSE ~ "No campaign dose reported",
      .default = "Campaign answer missing"
    ),
    # three mutually exclusive groups for coverage; NA when unknown
    status3 = case_when(
      vacc %in% FALSE ~ "Reported never vaccinated",
      camp_rep %in% TRUE ~ "Campaign dose reported",
      camp_rep %in% FALSE ~ "Vaccinated, no campaign dose reported",
      .default = NA
    )
  )
stopifnot(!any(is.na(df_vacc$elig_camp) & df_vacc$age_eligible %in% TRUE))

# campaign-eligible children; children aged 15 or younger at interview who
# were too young or too old during the campaigns; children within the target
# ages during the campaigns but aged 16 or more at interview
df_vacc |>
  summarise(
    n_eligible = sum(elig_camp, na.rm = TRUE),
    n_age_eligible_not_at_campaigns = sum(age_eligible & !elig_camp, na.rm = TRUE),
    n_aged_out = sum(aged_out, na.rm = TRUE)
  )

des_vacc <- as_kano_design(df_vacc)
des_vacc_el <- des_vacc |> filter(elig_camp %in% TRUE)

# Survey-weighted proportion with a logit CI on the full design's df (srvyr's
# survey_mean ignores df = on a filtered design, so svyciprop is called)
ciprop <- function(des, var) {
  p <- svyciprop(var, design = des, method = "logit", df = design_df)
  tibble(est = as.numeric(p), est_low = attr(p, "ci")[[1]], est_upp = attr(p, "ci")[[2]])
}

# Coverage: three mutually exclusive groups on one denominator, campaign-
# eligible children whose vaccination status, and campaign answer if
# vaccinated, are recorded
des_vacc_s3 <- des_vacc_el |> filter(!is.na(status3))
coverage <-
  c("Reported never vaccinated", "Campaign dose reported", "Vaccinated, no campaign dose reported") |>
  map_dfr(\(g) {
    ciprop(des_vacc_s3 |> mutate(.g = status3 == g), ~.g) |>
      mutate(indicator = g, n = sum(des_vacc_s3$variables$status3 == g), .before = 1)
  }) |>
  mutate(N = sum(n), n_unknown = sum(df_vacc$elig_camp %in% TRUE & is.na(df_vacc$status3)))
stopifnot(abs(sum(coverage$est) - 1) < 1e-8)
print(coverage, width = Inf)
write_csv(coverage, file.path(tab_dir, "coverage_campaign_eligible.csv"))

# Table 4: recorded vaccination history by the caregiver's campaign answer
# (counts unweighted), with the survey-weighted % reporting a campaign dose
# among vaccinated children with an answer
table_4_counts <-
  df_vacc |>
  filter(elig_camp %in% TRUE) |>
  count(history, camp_lbl) |>
  pivot_wider(names_from = camp_lbl, values_from = n, values_fill = 0) |>
  mutate(Total = rowSums(across(where(is.numeric))))

des_vacc_camp <- des_vacc_el |> filter(vacc %in% TRUE, !is.na(camp_rep))
table_4_pct <-
  bind_rows(
    map_dfr(hist_levels[2:5], \(h) ciprop(des_vacc_camp |> filter(history == h), ~camp_rep) |> mutate(history = h)),
    ciprop(des_vacc_camp, ~camp_rep) |> mutate(history = "All reported vaccinated")
  ) |>
  rename(pct_camp = est, pct_camp_low = est_low, pct_camp_upp = est_upp)

table_4 <-
  table_4_counts |>
  mutate(history = as.character(history)) |>
  bind_rows(
    table_4_counts |>
      filter(!history %in% c("Reported never vaccinated", "Vaccination status missing")) |>
      summarise(history = "All reported vaccinated", across(where(is.numeric), sum))
  ) |>
  left_join(table_4_pct, by = "history")
stopifnot(sum(table_4_counts$Total) == sum(df_vacc$elig_camp, na.rm = TRUE))
print(as.data.frame(table_4))
write_csv(table_4, file.path(tab_dir, "table4_vaccination_history_campaign.csv"))

# Table 4 footnote: campaign answers given outside the skip logic (never
# vaccinated, or vaccination status missing), by interview date; not used
df_vacc |>
  filter(elig_camp %in% TRUE, !(vacc %in% TRUE), !is.na(vacc_camp)) |>
  count(history, date_int, vacc_camp = as.character(vacc_camp))

# S3 Table: main reason given for no campaign dose (vaccinated eligible
# children whose caregiver answered No); unweighted counts and survey-weighted
# %, logit CI on the design df
des_vacc_no <-
  des_vacc_el |>
  filter(camp_rep %in% FALSE) |>
  mutate(reason = coalesce(as.character(vacc_no_reason_001), "Reason missing"))
s3_table <-
  des_vacc_no$variables |>
  count(reason, sort = TRUE) |>
  pull(reason) |>
  map_dfr(\(r) {
    ciprop(des_vacc_no |> mutate(.r = reason == r), ~.r) |>
      mutate(reason = r, n = sum(des_vacc_no$variables$reason == r), N = nrow(des_vacc_no$variables), .before = 1)
  })
print(as.data.frame(s3_table), digits = 3)
write_csv(s3_table, file.path(tab_dir, "S3_table_no_campaign_reasons.csv"))

# Reported vaccination by whether the household respondent believes that
# vaccination protects against disease: counts, then the odds ratio adjusted for
# sex, age group and SES quintile
df_ind |>
  filter(!is.na(vacc), !is.na(vacc_belief_protect)) |>
  group_by(vacc_belief_protect) |>
  summarise(n = n(), n_vaccinated = sum(vacc), pct = n_vaccinated / n * 100)

svyglm(vacc ~ sex + age_grp_diph + vacc_belief_protect + ses_quintile, family = quasibinomial, design = df_srvy_ind) |>
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) |>
  filter(term == "vacc_belief_protectYes")

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

figure_1 <-
  df_epi |>
  ggplot(aes(month_ym, n)) +
  # One colour-blind-safe fill for all LGAs; the facet strips name the LGA
  geom_col(width = 1, colour = "black", fill = "#4575b4") +
  facet_wrap(~lga, ncol = 1, strip.position = "right") +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x), expand = c(0, 0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  labs(x = "Month of disease onset", y = "Number of individuals") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey90", colour = NA), strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "figure1_epicurve.jpg"), figure_1, width = 12, height = 8, dpi = 300)

# ============================================================
# 10. Spatial figures (Figure 2, S1 Figure) and cluster-level
# correlation tests
# ============================================================
# Cluster-level rates are empirical-Bayes shrunk toward the survey-wide
# value before kernel-smoothing, since raw small-cluster rates would
# otherwise produce spurious hotspots; the formal correlation tests
# below use the raw (unshrunk) cluster estimates instead, so shrinkage
# doesn't artificially inflate them. Vaccination indicators are on the
# campaign-eligible population (Section 8), as shares of children whose
# status is recorded: campaign dose reported, vaccinated with no campaign
# dose reported, and never vaccinated sum to 100% in every cluster; ever
# vaccinated is among children with a recorded vaccination status.

df_vacc_map <-
  df_vacc |>
  filter(!is.na(cluster), elig_camp %in% TRUE) |>
  mutate(
    cov_camp_i = if_else(is.na(status3), NA, status3 == "Campaign dose reported"),
    cov_notcamp_i = if_else(is.na(status3), NA, status3 == "Vaccinated, no campaign dose reported"),
    cov_never_i = if_else(is.na(status3), NA, status3 == "Reported never vaccinated"),
    cov_ever_i = vacc
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
  df_ind |>
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
  df_ind |>
  filter(!is.na(cluster), age_grp %in% "0-4") |>
  group_by(cluster) |>
  summarise(num_ppl_u5 = n(), died_u5_n = sum(died %in% TRUE, na.rm = TRUE), pt_u5_tot = sum(pt, na.rm = TRUE), u5mr = died_u5_n / pt_u5_tot * 10000)

tbl_cluster_vacc <-
  df_vacc_map |>
  group_by(cluster) |>
  summarise(
    num_ppl_vacc = n(),
    cov_camp = mean(cov_camp_i, na.rm = TRUE) * 100,
    cov_notcamp = mean(cov_notcamp_i, na.rm = TRUE) * 100,
    cov_never = mean(cov_never_i, na.rm = TRUE) * 100,
    cov_ever = mean(cov_ever_i, na.rm = TRUE) * 100
  )

tbl_cluster <- tbl_cluster_core |> full_join(tbl_cluster_u5, by = "cluster") |> left_join(tbl_cluster_vacc, by = "cluster")
stopifnot(nrow(tbl_cluster) == 78)

# unshrunk, for the Spearman correlations below
tbl_cluster_raw <- tbl_cluster

apply_shrinkage <- function(tbl, raw_col, shrunk_tbl) {
  if (is.null(shrunk_tbl)) return(tbl)
  tbl |> left_join(shrunk_tbl, by = "cluster") |> mutate("{raw_col}" := coalesce(shrunk, .data[[raw_col]])) |> select(-shrunk)
}

df_ind_u5 <- df_ind |> filter(age_grp %in% "0-4")

tbl_cluster <-
  tbl_cluster |>
  apply_shrinkage("attack_rate", shrink_binom(df_ind, diphtheria)) |>
  apply_shrinkage("mr", shrink_rate(df_ind, died)) |>
  apply_shrinkage("u5mr", shrink_rate(df_ind_u5, died)) |>
  apply_shrinkage("cov_camp", shrink_binom(df_vacc_map, cov_camp_i)) |>
  apply_shrinkage("cov_notcamp", shrink_binom(df_vacc_map, cov_notcamp_i)) |>
  apply_shrinkage("cov_never", shrink_binom(df_vacc_map, cov_never_i))

tbl_cluster_gps <-
  df_ind |>
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

# cluster-level correlation tests, raw (unshrunk) cluster estimates: each
# vaccination indicator against the attack rate, all-cause and under-5
# mortality (never vaccinated gives the ever-vaccinated values with the sign
# reversed, so it is not tested separately)
spearman_test <- function(x_col, y_col, label) {
  test <- suppressWarnings(cor.test(tbl_cluster_raw[[x_col]], tbl_cluster_raw[[y_col]], method = "spearman"))
  tibble(pair = label, rho = unname(test$estimate), p = test$p.value)
}

spearman_pairs <- tribble(
  ~x_col,        ~y_col,        ~label,
  "cov_ever",    "attack_rate", "Ever vaccinated vs diphtheria attack rate",
  "cov_ever",    "mr",          "Ever vaccinated vs all-cause mortality",
  "cov_ever",    "u5mr",        "Ever vaccinated vs under-5 mortality",
  "cov_notcamp", "attack_rate", "Vaccinated, no campaign dose vs diphtheria attack rate",
  "cov_notcamp", "u5mr",        "Vaccinated, no campaign dose vs under-5 mortality",
  "cov_camp",    "attack_rate", "Campaign dose reported vs diphtheria attack rate",
  "cov_camp",    "u5mr",        "Campaign dose reported vs under-5 mortality"
)

spatial_correlations <- pmap_dfr(spearman_pairs, spearman_test)
print(spatial_correlations, n = Inf)
write_csv(spatial_correlations, file.path(tab_dir, "spatial_correlations.csv"))

# indicator: column in tbl_cluster; weight_col: cluster-N column used to
# weight the smoothing. Viridis throughout (begin/end trim the extremes so
# contour labels stay legible); direction = -1 makes high values dark.
# contour_by: contour interval; overlapping contour labels are dropped.
plot_ppp_srvy <- function(dat, shp, indicator, weight_col, legend_name, contour_by, direction = 1, begin = 0.15, end = 0.97) {
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
  max_v <- max(grid_points$v, na.rm = TRUE)
  contour_breaks <- seq(0, ceiling(max_v / contour_by) * contour_by, by = contour_by)

  ggplot() +
    geom_sf(data = sm_dat, aes(fill = v), colour = NA) +
    geom_sf(data = st_boundary(shp_lga), colour = "grey40", linewidth = 0.4) +
    geom_sf(data = st_boundary(shp), colour = "black", linewidth = 0.8) +
    geom_contour(data = grid_points, aes(x = lon, y = lat, z = v), breaks = contour_breaks, colour = "grey20", linewidth = 0.3, linetype = "dashed") +
    metR::geom_text_contour(data = grid_points, aes(x = lon, y = lat, z = v), breaks = contour_breaks, colour = "black", size = 2, skip = 0, stroke = 0.1, check_overlap = TRUE) +
    scale_fill_viridis_c(
      option = "viridis", direction = direction, begin = begin, end = end,
      name = legend_name, limits = c(0, max_v)
    ) +
    coord_sf(crs = st_crs(shp)) +
    theme_void()
}

# PLOS-ready figures: TIFF, LZW, 300 dpi, 7.5 in wide, Arial 8-10 pt, no title
# inside the figure file
plos_fig_dir <- file.path(fig_dir, "plos")
dir.create(plos_fig_dir, showWarnings = FALSE, recursive = TRUE)

map_theme <-
  theme(
    legend.position = "bottom", text = element_text(family = "Arial", size = 9),
    plot.title = element_text(family = "Arial", size = 10, face = "bold"),
    legend.text = element_text(family = "Arial", size = 8),
    legend.title = element_text(family = "Arial", size = 8)
  )

save_plos_tif <- function(plot, file, height) {
  ggsave(file.path(plos_fig_dir, file), plot, device = ragg::agg_tiff, compression = "lzw", width = 7.5, height = height, units = "in", dpi = 300, bg = "white")
}

# Figure 2: the three mutually exclusive vaccination groups among
# campaign-eligible children, each panel on its own scale; darker always
# means less vaccinated (the two vaccinated groups run dark = low, never
# vaccinated dark = high)
fig2_panel <- \(ind, lab, direction, contour_by) {
  plot_ppp_srvy(tbl_cluster, shp_boundary, ind, "num_ppl_vacc", lab, contour_by, direction = direction) +
    guides(fill = guide_colourbar(title.position = "top", barwidth = unit(4.2, "cm"), barheight = unit(0.3, "cm")))
}
figure_2 <-
  (fig2_panel("cov_camp", "Campaign dose reported (%)", 1, 10) + labs(title = "a")) +
  (fig2_panel("cov_notcamp", "Vaccinated, no campaign dose (%)", 1, 5) + labs(title = "b")) +
  (fig2_panel("cov_never", "Never vaccinated (%)", -1, 10) + labs(title = "c")) +
  plot_layout(ncol = 3) &
  map_theme

ggsave(file.path(fig_dir, "figure2_vaccination_coverage.jpg"), figure_2, width = 7.5, height = 3.3, units = "in", dpi = 300, bg = "white")
save_plos_tif(figure_2, "Fig2.tif", height = 3.3)

# S1 Figure: diphtheria attack rate and under-5 mortality, each panel on its
# own scale, lighter = higher. No vaccination panel: ever vaccinated is 100%
# minus Figure 2c.
s1_panel <- \(ind, w, lab, contour_by) {
  plot_ppp_srvy(tbl_cluster, shp_boundary, ind, w, lab, contour_by, direction = 1) +
    guides(fill = guide_colourbar(title.position = "top", barwidth = unit(5.5, "cm"), barheight = unit(0.3, "cm")))
}
s1_figure <-
  (s1_panel("attack_rate", "num_ppl", "Diphtheria attack rate (%)", 0.5) + labs(title = "a")) +
  (s1_panel("u5mr", "num_ppl_u5", "Under-5 mortality rate (per 10,000 person-days)", 0.2) + labs(title = "b")) +
  plot_layout(ncol = 2) &
  map_theme

ggsave(file.path(fig_dir, "S1_figure_attack_rate_u5mr.jpg"), s1_figure, width = 7.5, height = 3.9, units = "in", dpi = 300, bg = "white")
save_plos_tif(s1_figure, "S1_Fig.tif", height = 3.9)

# PLOS-ready Figure 1, drawn at final size from the epicurve above
plos_theme <-
  theme(
    text = element_text(family = "Arial", size = 9),
    plot.title = element_text(family = "Arial", size = 10, face = "plain"),
    axis.text = element_text(family = "Arial", size = 8),
    axis.title = element_text(family = "Arial", size = 9),
    legend.text = element_text(family = "Arial", size = 8),
    legend.title = element_text(family = "Arial", size = 9),
    strip.text = element_text(family = "Arial", size = 9)
  )

save_plos_tif(figure_1 + plos_theme, "Fig1.tif", height = 5.5)
