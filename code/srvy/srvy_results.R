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
# df_hh: one row per approached household (n=1,083) - used for the
#   response-rate figure
# map_boundaries: LGA and study-area boundary polygons, used only by the
#   spatial figures at the end of this script

df_ind <- read_rds("data/df_ind.rds") |>
  mutate(
    # Table 4 EPI status: a card dated inside a campaign window (vacc_epi_confirmed
    # NA, vacc_source camp_only) counts as a campaign dose, not a routine one, so
    # those children sit in the "No reported routine EPI vaccination" row and
    # Table 4 holds all 3,287 children with a classified vaccination source.
    epi_t4 = if_else(is.na(vacc_epi_confirmed) & vacc_source == "camp_only", FALSE, vacc_epi_confirmed),
    # restricted case definition: cases also reported as diagnosed by a health professional
    diphtheria_corrob = diphtheria == TRUE & diagnosis == "Yes"
  )
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

# enrolled households, clusters, and LGA-by-cluster units (the primary sampling units)
df_ind |> summarise(n_households = n_distinct(index_hh), n_clusters = n_distinct(cluster), n_lga_clusters = n_distinct(paste(lga, cluster)))

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

# incidence rate per 1,000 person-years (case count / total person-time)
df_srvy_ind |>
  summarise(
    ir = survey_ratio(numerator = diphtheria * 1000 * 365.25, denominator = pt, vartype = "ci", na.rm = TRUE)
  )

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
  group_by(index_hh, cluster, lga, weight) |>
  summarise(hh_size = n(), num_cases = sum(diphtheria, na.rm = TRUE), .groups = "drop") |>
  filter(num_cases > 0) |>
  mutate(secondary_cases = pmax(num_cases - 1, 0), exposed = pmax(hh_size - 1, 0), sar = secondary_cases / exposed)

nrow(hh_summary)
sum(hh_summary$num_cases >= 2)

hh_design <- hh_summary |> as_survey(ids = c(cluster, index_hh), strata = lga, weights = weight, nest = TRUE)
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
    filter(!is.na(cluster), !is.na(index_hh)) |>
    as_survey(ids = c(cluster, index_hh), strata = lga, weights = weight, nest = TRUE) |>
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
# reliable estimate on their own)
df_srvy_ind |>
  filter(died == TRUE) |>
  group_by(died_location) |>
  summarise(pct = survey_mean(vartype = "ci") * 100)

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
# with a delta-method CI on the log scale that treats the groups as independent.
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

# mortality rate ratios, unvaccinated vs vaccinated (reference), whole surveyed
# population; p-values are Wald z on the log scale, from the same SE as the CI
rate_ratios <-
  bind_rows(
    calculate_rate_ratio(df_srvy_ind, died, pt, vacc, grp_ref = TRUE, grp_compare = FALSE) |> mutate(outcome = "All-cause"),
    calculate_rate_ratio(df_srvy_ind, died_diph_binary, pt, vacc, grp_ref = TRUE, grp_compare = FALSE) |> mutate(outcome = "Diphtheria-attributed")
  ) |>
  mutate(p = 2 * pnorm(-abs(log(rate_ratio) / ((log(ci_high) - log(ci_low)) / (2 * 1.96))))) |>
  select(outcome, everything())
print(rate_ratios)

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
# The comparator for every restricted analysis below is explicit
# self-report-unvaccinated (vacc == FALSE), not a derived binary's own
# FALSE category - e.g. vacc_epi_binary == FALSE means "reports
# vaccinated but no card evidence", not "unvaccinated".

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
# Campaign-eligible children are those aged 6 months to 15 years (age_eligible).

df_ind |>
  filter(age_eligible) |>
  summarise(
    n_eligible = n(),
    n_status_missing = sum(is.na(vacc)),
    n_reported_unvaccinated = sum(vacc == FALSE, na.rm = TRUE),
    n_reported_vaccinated = sum(vacc == TRUE, na.rm = TRUE),
    n_epi_status_known = sum(!is.na(vacc_epi_confirmed)),
    n_camp_status_known = sum(!is.na(camp_confirmed)),
    n_source_unknown = sum(vacc_source == "unknown", na.rm = TRUE),
    pct_source_unknown = n_source_unknown / n_eligible * 100
  )

# Vaccination evidence, and the Table 4 row each group falls in. Children with
# missing status, or reported as vaccinated with no card and no campaign dose,
# have no classifiable source and are excluded from Table 4.
vacc_evidence <-
  df_ind |>
  filter(age_eligible) |>
  mutate(
    evidence = case_when(
      vacc_card == "Yes" & is.na(vacc_epi_confirmed) ~ "Card seen, dated inside a campaign window",
      vacc_card == "Yes" ~ "Card seen, dated outside campaign windows",
      vacc_card == "Yes, but it is missing/lost" ~ "Reported vaccinated, card lost",
      vacc_card == "Was vaccinated in the DTC/OPD/Contact clinic" ~ "Reported vaccinated at DTC/OPD/contact clinic",
      vacc_card == "No" ~ "Vaccinated, no card",
      vacc == FALSE ~ "Self-reported unvaccinated (card not asked)",
      .default = "Vaccination status missing"
    ),
    table4_group = case_when(
      is.na(vacc) ~ "Excluded: vaccination status missing",
      vacc_source == "unknown" ~ "Excluded: source unclassifiable",
      epi_t4 %in% TRUE ~ "Reported routine EPI vaccination",
      .default = "No reported routine EPI vaccination"
    )
  ) |>
  count(evidence, table4_group) |>
  mutate(pct = round(100 * n / sum(n), 1))
print(vacc_evidence)
write_csv(vacc_evidence, file.path(tab_dir, "vaccination_evidence.csv"))

# survey-weighted routine EPI / campaign coverage, among age-eligible children
# with the status observed
df_srvy_ind |>
  filter(age_eligible) |>
  summarise(
    epi_coverage = survey_mean(vacc_epi_confirmed == TRUE, vartype = "ci", na.rm = TRUE),
    camp_coverage = survey_mean(camp_confirmed == TRUE, vartype = "ci", na.rm = TRUE)
  )

# evidence of vaccination from at least one source, among children with
# a classified vaccination source
df_srvy_ind |>
  filter(age_eligible, vacc_source != "unknown") |>
  summarise(any_source = survey_mean(vacc_source %in% c("both", "epi_only", "camp_only"), vartype = "ci", na.rm = TRUE))

# no vaccination of any kind, among age-eligible children with
# self-reported vaccination status observed
df_srvy_ind |>
  filter(age_eligible, !is.na(vacc)) |>
  summarise(n = unweighted(n()), none = survey_mean(vacc == FALSE, vartype = "ci", na.rm = TRUE))

# Table 4: reported routine EPI vaccination against campaign vaccination, among
# children with a classified vaccination source (counts unweighted, row
# percentages of campaign-vaccinated)
df_cc <- df_ind |> filter(age_eligible, vacc_source != "unknown", !is.na(epi_t4), !is.na(camp_confirmed))
nrow(df_cc)

table_4 <-
  df_cc |>
  mutate(
    epi_label = factor(epi_t4, levels = c(TRUE, FALSE), labels = c("Reported routine EPI vaccination", "No reported routine EPI vaccination")),
    camp_label = if_else(camp_confirmed, "Campaign confirmed", "Campaign not confirmed")
  ) |>
  count(epi_label, camp_label) |>
  pivot_wider(names_from = camp_label, values_from = n, values_fill = 0L)
table_4 <-
  table_4 |>
  bind_rows(table_4 |> summarise(epi_label = "Total", across(where(is.numeric), sum))) |>
  mutate(total = `Campaign confirmed` + `Campaign not confirmed`, pct_campaign = `Campaign confirmed` / total * 100)
print(table_4)
write_csv(table_4, file.path(tab_dir, "table4_epi_vs_campaign.csv"))

# design-based (Rao-Scott) test of routine EPI status against campaign status
des_cc <- df_srvy_ind |> filter(age_eligible, vacc_source != "unknown", !is.na(epi_t4), !is.na(camp_confirmed))
svychisq(~ epi_t4 + camp_confirmed, design = des_cc)

# survey-weighted: campaign recipients with reported routine EPI vaccination
df_srvy_ind |> filter(age_eligible, camp_confirmed == TRUE) |> summarise(epi_among_campaign = survey_mean(epi_t4 == TRUE, vartype = "ci", na.rm = TRUE))

# survey-weighted: children without reported routine EPI vaccination who
# received a campaign dose
des_cc |>
  filter(epi_t4 == FALSE) |>
  summarise(campaign_among_no_epi = survey_mean(camp_confirmed == TRUE, vartype = "ci", na.rm = TRUE))

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

# PLOS-ready Figure 1 (TIFF, LZW, 300 dpi, 7.5 in wide, Arial 8-10 pt), drawn at
# final size from the epicurve above.
plos_fig_dir <- file.path(fig_dir, "plos")
dir.create(plos_fig_dir, showWarnings = FALSE, recursive = TRUE)

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

ggsave(
  file.path(plos_fig_dir, "Fig1.tif"), figure_1 + plos_theme,
  device = ragg::agg_tiff, compression = "lzw",
  width = 7.5, height = 5.5, units = "in", dpi = 300, bg = "white"
)
