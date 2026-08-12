# Builds the three de-identified analysis-ready files hbc_results.R reads
# (data/df_pt.rds, data/df_cohort_samp.rds, data/df_hh.rds) from this
# project's internal frozen cohort snapshot. Not part of the public
# release itself - run once internally, then share the resulting rds
# files under the MSF data access process described in the README.

pacman::p_load(here, tidyverse, haven, dm)

load(here("output/hbc/data/cohort_sample_final.rds"))
db_raw <- read_rds(here("data/hbc/cohort_study_kobo_snapshot_2026-07-23.rds"))

df_subclass <- df_samp |> select(case_id, subclass, weights)

df_pt <-
  db_raw$main |>
  select(
    index = `_index`,
    date_int,
    case_id = id_pt,
    age,
    sex,
    education,
    vacc,
    vacc_camp,
    vacc_card,
    starts_with(c("symptoms_", "sequelae_")),
    sequelae,
    tx_received,
    tx_location,
    tx_improve,
    tx_acceptable,
    tx_quality,
    outcome,
    delay_tx,
    hh_size,
    room_count
  ) |>
  mutate(
    across(where(is.character), as_factor),
    case_id = str_extract(case_id, "^[^:]+"),
    tx_location = case_match(tx_location, "hosp" ~ "hospital", .default = tx_location),
    tx_mod = factor(if_else(tx_location == "home", "HBC", "DTC"), levels = c("HBC", "DTC")),
    died = case_match(
      as.character(outcome),
      "Died" ~ TRUE, "Recovered" ~ FALSE, "Still unwell/complications" ~ FALSE,
      .default = NA
    ),
    outcome_cat = as.character(outcome),
    vaccinated = case_match(as.character(vacc), "Yes" ~ 1L, "No" ~ 0L, .default = NA_integer_),
    vacc_card_confirmed = case_when(
      as.character(vacc) == "No" ~ 0L,
      as.character(vacc) == "Yes" &
        as.character(vacc_card) %in% c("Yes", "Yes, but it is missing/lost", "Was vaccinated in the DTC/OPD/Contact clinic") ~ 1L,
      as.character(vacc) == "Yes" & as.character(vacc_card) == "No" ~ 0L,
      TRUE ~ NA_integer_
    ),
    age_grp = cut(age, breaks = c(0, 5, 10, 15, Inf), labels = c("0-5", "6-10", "11-15", "16+"), right = TRUE, include.lowest = TRUE),
    delay_ge4 = case_when(delay_tx >= 4 ~ TRUE, delay_tx < 4 ~ FALSE, TRUE ~ NA),
    score_improve = case_match(as.character(tx_improve), "Strongly disagree" ~ 1, "Disagree" ~ 2, "Neutral" ~ 3, "Agree" ~ 4, "Strongly agree" ~ 5, .default = NA_real_),
    score_acceptable = case_match(as.character(tx_acceptable), "Completely unacceptable" ~ 1, "Unacceptable" ~ 2, "Neutral/no opinion" ~ 3, "Acceptable" ~ 4, "Completely acceptable" ~ 5, .default = NA_real_),
    score_quality = case_match(as.character(tx_quality), "Strongly negative" ~ 1, "Quite negative" ~ 2, "Neutral" ~ 3, "Quite positive" ~ 4, "Strongly positive" ~ 5, .default = NA_real_),
    score_mean = (score_improve + score_acceptable + score_quality) / 3
  ) |>
  left_join(df_subclass, by = "case_id")

df_pt <- df_pt |>
  left_join(
    df_cohort |>
      mutate(complication = !is.na(complication_1), comorbidity = !is.na(comorbidity_1)) |>
      select(
        case_id,
        lga = ref_adm2_name,
        date_onset = date_notification,
        complication,
        comorbidity,
        vacc_db = vacci_diphtheria_doses,
        sex_reg = sex,
        age_grp_reg = age_group,
        vacc_status_reg = vacci_diphtheria_status
      ),
    by = "case_id"
  ) |>
  mutate(
    vaccinated = case_when(!is.na(vaccinated) ~ vaccinated, vacc_db == 0 ~ 0L, vacc_db > 0 ~ 1L, TRUE ~ NA_integer_),
    lga = case_when(
      !is.na(lga) ~ lga,
      str_detect(case_id, "^MM") ~ "Ungogo",
      str_detect(case_id, "^ID") ~ "Fagge",
      str_detect(case_id, "^DT") ~ "Dawakin Tofa",
      str_detect(case_id, "^NS") ~ "Nassarawa",
      TRUE ~ NA_character_
    )
  ) |>
  select(-tx_received, -education, -vacc_card, -tx_improve, -tx_acceptable, -tx_quality, -outcome, -vacc, -vacc_camp, -vacc_db, -date_int)

df_hh <-
  db_raw$group_hh |>
  select(parent_index = `_parent_index`, diphtheria_hh, timing_post, vacc_hh, age_hh, sex_hh) |>
  left_join(
    df_pt |> select(index, tx_mod, died, vaccinated, delay_tx, delay_ge4, complication, lga, hh_size, room_count, case_id),
    by = c("parent_index" = "index")
  ) |>
  mutate(
    hh_case = as.integer(diphtheria_hh == "1"),
    hh_vacc = as.integer(vacc_hh == "1"),
    timing_days = as.numeric(timing_post),
    age_hh_f = na_if(haven::as_factor(age_hh), "Do not know"),
    sex_hh_f = haven::as_factor(sex_hh)
  ) |>
  select(-diphtheria_hh, -timing_post, -vacc_hh, -age_hh, -sex_hh)

df_pt_public <- df_pt |> select(-index)

df_cohort_samp_public <- df_cohort_samp |>
  select(case_id, tx_location, sex, age_group, vaccination_status = vacci_diphtheria_status)

dir.create(here("code/public_release/data"), showWarnings = FALSE, recursive = TRUE)
write_rds(df_pt_public, here("code/public_release/data/df_pt.rds"))
write_rds(df_cohort_samp_public, here("code/public_release/data/df_cohort_samp.rds"))
write_rds(df_hh, here("code/public_release/data/df_hh.rds"))

cat(sprintf(
  "Wrote df_pt (%d rows), df_cohort_samp (%d rows), df_hh (%d rows) to code/public_release/data/\n",
  nrow(df_pt_public), nrow(df_cohort_samp_public), nrow(df_hh)
))
