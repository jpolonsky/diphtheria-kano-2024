# Builds the de-identified analysis-ready files srvy_results.R and
# srvy_maps_public.R read (data/df_ind.rds, data/df_hh.rds,
# data/map_boundaries.rds) from this project's internal raw Kobo pull.
# Not part of the public release itself - run once internally, then
# share the resulting rds files under the MSF data access process
# described in the README.

pacman::p_load(here, sf, janitor, mice, dm, haven, tidyverse)

db_raw <- read_rds(here("data/srvy/db_raw_svy.rds"))

df_hh_raw <- db_raw$main |> clean_names() |> select(!where(is.list))

df_hh <-
  df_hh_raw |>
  select(
    index_hh = index,
    date_int,
    team = id_team,
    cluster,
    consent,
    hh = hh_num,
    starts_with("hh_effects_"),
    property,
    land,
    starts_with("animals_"),
    education:vacc_refused,
    vacc_refused_reasons,
    vacc_prevent_distance:vacc_prevent_cost,
    vacc_groups:vacc_comm_leaders,
    hp_visit:diphtheria_count,
    hh_size,
    num_left = hh_size_left,
    num_died = hh_size_died,
    lat = gps_latitude,
    lon = gps_longitude
  ) |>
  mutate(
    across(where(\(x) !is.null(attr(x, "labels"))), \(x) x |> as_factor()),
    team = team |> fct_relevel(str_c("Team ", 1:13)),
    hh = hh |> as.numeric()
  )

map_nga_ward <-
  here("data/srvy/maps/nga_adm3/Nigeria_-_Ward_Boundaries.shp") |>
  st_read()

map_nga_ward <-
  c("Dawakin Tofa", "Fagge", "Nassarawa", "Ungogo") |>
  set_names() |>
  map(\(x) {
    map_nga_ward |>
      filter(statename %in% c("Kano"), lganame %in% x) |>
      st_union() |>
      st_as_sf()
  }) |>
  list_rbind(names_to = "lga") |>
  st_as_sf()

map_nga_ward_union <- map_nga_ward |> st_union() |> st_as_sf()

df_gis_hh <-
  df_hh |>
  select(index_hh, lat, lon) |>
  filter(!is.na(lat), !is.na(lon)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

df_hh <- df_hh |> left_join(df_gis_hh) |> st_as_sf()

df_present_raw <-
  db_raw$Demography_present |>
  clean_names() |>
  select(index_ind = index, index_hh = parent_index, age:sequelae, date_onset:symptoms, diagnosis:tx_hf_no, starts_with("tx_hf_no_")) |>
  rename(join = join_new) |>
  add_column(status = "present")

df_left_raw <-
  db_raw$Demography_departed |>
  clean_names() |>
  rename_with(\(x) x |> str_remove("_left")) |>
  select(index_ind = index, index_hh = parent_index, age:sequelae, date_onset:symptoms, diagnosis:tx_hf_no) |>
  add_column(status = "left")

df_died_raw <-
  db_raw$Demography_deceased |>
  clean_names() |>
  rename(date_death = date_died) |>
  rename_with(\(x) x |> str_remove("_died")) |>
  select(index_ind = index, index_hh = parent_index, age:sequelae, date_onset:symptoms, diagnosis:tx_hf_no) |>
  rename(date_died = date_death) |>
  add_column(status = "died")

df_ind <- bind_rows(df_present_raw, df_left_raw, df_died_raw)
df_ind <- left_join(df_hh, df_ind)
df_ind <- df_ind |> mutate(across(where(\(x) !is.null(attr(x, "labels"))), \(x) x |> as_factor()))

df_ind <-
  df_ind |>
  mutate(
    education_binary = education |>
      case_match(
        c("No education", "Non-formal eductaion (e.g. Islamic)", "Some primary", "Completed primary") ~ 0,
        c("Some secondary", "Completed secondary", "More than secondary") ~ 1,
        .default = NA_real_
      ),
    # vacc/vacc_pre coerced to logical FIRST, before anything downstream
    # compares them to FALSE/TRUE - reordering this fixed a real bug (see
    # CLAUDE.md "Key analytical decisions", SUPERSEDED 2026-08-11) where
    # vacc_camp_binary/vacc_epi_binary silently never matched self-report-
    # unvaccinated individuals, because the comparison ran against the
    # still-uncoerced Kobo factor.
    vacc = vacc |> case_match("No" ~ FALSE, "Yes" ~ TRUE, .default = NA),
    vacc_pre = vacc_pre |> case_match("No" ~ FALSE, "Yes" ~ TRUE, .default = NA),
    vacc_epi_binary = case_when(
      vacc == FALSE ~ FALSE,
      vacc_card == "No" ~ FALSE,
      vacc_card %in% c("Yes", "Yes, but it is missing/lost", "Was vaccinated in the DTC/OPD/Contact clinic") ~ TRUE,
      .default = NA
    ),
    vacc_camp_binary = case_when(
      vacc == FALSE ~ FALSE,
      vacc_camp == "No" ~ FALSE,
      vacc_camp |> str_detect("Yes") ~ TRUE,
      .default = NA
    ),
    vacc_ever = case_when(
      vacc_epi_binary == TRUE | vacc_camp_binary == TRUE ~ TRUE,
      vacc_epi_binary == FALSE & vacc_camp_binary == FALSE ~ FALSE,
      .default = NA
    ),
    diphtheria = case_when(
      diphtheria == "Yes" | died_diph == "Yes" ~ TRUE,
      diphtheria == "No" ~ FALSE,
      .default = NA
    ),
    died = case_when(!is.na(date_died) ~ TRUE, .default = FALSE),
    died_cause = case_when(died_diph == "Yes" ~ "Diphtheria", .default = died_cause) |>
      fct_relevel(c("Diphtheria", levels(df_ind$died_cause))),
    died_diph_binary = case_when(died_cause == "Diphtheria" ~ TRUE, .default = FALSE),
    tx_seeking_binary = case_when(
      died_diph == "Yes" ~ tx_seeking,
      status == "died" ~ tx_seeking_hsb,
      .default = tx_seeking
    ),
    tx_seeking_binary = tx_seeking_binary |> case_match("Do not know" ~ NA, .default = tx_seeking_binary),
    tx_received_binary = tx_received |>
      case_match(
        c("None", "Other/traditional treatment", "Spritual assistance") ~ FALSE,
        c("Orthodox medicine", "Immediate referral") ~ TRUE,
        .default = NA
      ),
    delay_tx_binary = case_when(delay_tx >= 4 ~ TRUE, delay_tx < 4 ~ FALSE, .default = NA)
  )

recall_begin <- as_date("2023-01-01")

df_ind <-
  df_ind |>
  mutate(
    age = case_when(age_unit == "Month" ~ age / 12, TRUE ~ age),
    age_cat = age |> AMR::age_groups(split_at = 5),
    age_grp = age |> AMR::age_groups(split_at = "fives"),
    age_grp_diph = age |> AMR::age_groups(split_at = c(6, 11, 16)),
    cause_start = case_when(
      join == "Joined household during recall period" ~ "Joined",
      born == "Yes" ~ "Born",
      .default = "Present at start"
    ),
    cause_end = case_when(
      status == "left" ~ "Left",
      status == "died" ~ "Died",
      status == "present" ~ "Present at end",
      .default = NA
    ),
    pt_start = case_when(
      cause_start == "Present at start" ~ recall_begin,
      cause_start %in% c("Joined", "Born") & cause_end == "Left" ~ recall_begin + (date_int - recall_begin) / 4,
      cause_start %in% c("Joined", "Born") & cause_end == "Died" ~ recall_begin + (date_died - recall_begin) / 2,
      cause_start %in% c("Joined", "Born") ~ recall_begin + (date_int - recall_begin) / 2,
      .default = NA
    ),
    pt_end = case_when(
      cause_end == "Present at end" ~ date_int,
      cause_start %in% c("Joined", "Born") & cause_end == "Left" ~ date_int - (date_int - recall_begin) / 4,
      cause_end == "Left" ~ date_int - (date_int - recall_begin) / 2,
      cause_end == "Died" ~ date_died,
      .default = NA
    ),
    pt = as.numeric(pt_end - pt_start)
  )

df_hh <-
  df_hh |>
  st_join(map_nga_ward, join = st_within) |>
  mutate(
    lga = index_hh |>
      case_match(
        c(152, 153, 615) ~ "Nassarawa",
        731 ~ "Fagge",
        c(873:892) ~ "Dawakin Tofa",
        .default = lga
      ) |>
      fct_relevel("Dawakin Tofa", "Fagge", "Nassarawa", "Ungogo")
  )

df_ind <- df_ind |> left_join(df_hh |> as_tibble() |> select(index_hh, lga))

cluster_allocation <-
  df_hh |>
  filter(consent == "Yes") |>
  as_tibble() |>
  group_by(lga) |>
  count(cluster) |>
  summarise(num_clusters = n_distinct(cluster))

df_ind <- df_ind |> left_join(cluster_allocation)

tbl_pop <- tribble(
  ~lga,           ~pop_geo_tot,
  "Dawakin Tofa",  461798,
  "Fagge",         360368,
  "Nassarawa",     943502,
  "Ungogo",       1071443
)

df_ind <-
  df_ind |>
  left_join(tbl_pop) |>
  add_count(index_hh, name = "hh_size") |>
  add_count(lga, name = "strata_size") |>
  add_count(cluster, name = "cluster_size") |>
  mutate(lga = lga |> fct_relevel("Dawakin Tofa", "Fagge", "Nassarawa", "Ungogo"))

# stratum-specific weight: LGA's total population / individuals actually
# sampled in that LGA
df_ind <-
  df_ind |>
  filter(!is.na(lga)) |>
  group_by(lga) |>
  mutate(pop_geo_samp = n(), weight = pop_geo_tot / pop_geo_samp) |>
  ungroup()

df_hh <- df_hh |> as_tibble()
df_ind <- df_ind |> as_tibble()

# ---- Socioeconomic status (PCA on household assets) ----

df_assets <-
  df_ind |>
  select(starts_with("hh_effects_"), starts_with("animals_"), property, land, education, water_source, toilet_type) |>
  mutate(
    property_score = case_when(property == "Owns" ~ 1, property == "Rents" ~ 0, .default = NA_real_),
    land_score = case_when(land == "Yes" ~ 1, land == "No" ~ 0, .default = NA_real_),
    education_score = case_when(
      education == "No education" ~ 0,
      education == "Some primary" ~ 1,
      education == "Completed primary" ~ 2,
      education == "Some secondary" ~ 3,
      education == "Completed secondary" ~ 4,
      education == "More than secondary" ~ 5,
      education == "Non-formal eductaion (e.g. Islamic)" ~ 1,
      education == "Do not know" ~ NA_real_,
      .default = NA_real_
    ),
    water_source_score = case_when(
      water_source == "Surface" ~ 0, water_source == "Unimproved" ~ 1, water_source == "Improved" ~ 2,
      water_source == "Do not know" ~ NA_real_, .default = NA_real_
    ),
    toilet_type_score = case_when(
      toilet_type == "Open defecation" ~ 0, toilet_type == "Unimproved" ~ 1, toilet_type == "Improved" ~ 2,
      toilet_type == "Do not know" ~ NA_real_, .default = NA_real_
    )
  ) |>
  select(-animals_0, -animals_99, -property, -land, -education, -water_source, -toilet_type) |>
  drop_na()

df_assets_clean <-
  df_assets |>
  select(where(~ is.numeric(.) && !all(is.na(.)) && var(., na.rm = TRUE) > 0)) |>
  select(-starts_with("hh_effects"), -c(education_score, water_source_score, toilet_type_score))

row_indices <- attr(df_assets_clean, "row.names")
pca_results <- prcomp(df_assets_clean, center = TRUE, scale. = TRUE)

# ses_quintile is joined back onto the FULL df_ind (not a row-restricting
# slice) - the original script keeps a separate df_ind_filtered object for
# this and never overwrites df_ind itself, since df_ind's full N=7998 is
# what Table 1 and other unweighted descriptive tabulations are built
# from. Rows without complete asset data (excluded from the PCA) get
# ses_quintile = NA here; svyglm()'s own listwise deletion handles that
# correctly wherever ses_quintile enters a model.
df_ses <- tibble(
  .row = row_indices,
  ses_composite = pca_results$x[, 1],
  ses_quintile = pca_results$x[, 1] |> ntile(5) |> factor(levels = 1:5, ordered = TRUE)
)

df_ind <-
  df_ind |>
  mutate(.row = row_number()) |>
  left_join(df_ses, by = ".row") |>
  select(-.row)

# ---- Vaccination source triangulation (routine EPI vs campaign), used
# for Tables 3-4 and the coverage maps ----

coerce_logical <- function(x) {
  if (is.logical(x)) return(x)
  x_chr <- as.character(x)
  case_when(
    is.na(x_chr) ~ NA,
    x_chr %in% c("TRUE", "True", "true", "T", "Yes", "YES", "Y", "Yes (card observed)",
                 "Yes (card not observed)", "Yes (card)", "Yes, but it is missing/lost",
                 "Was vaccinated in the DTC/OPD/Contact clinic") ~ TRUE,
    x_chr %in% c("FALSE", "False", "false", "F", "No", "NO", "N") ~ FALSE,
    .default = NA
  )
}

in_any_phase <- function(date_vec, phases) {
  map_lgl(date_vec, \(d) {
    if (is.na(d)) return(NA)
    any(d >= phases$start & d <= phases$end, na.rm = TRUE)
  })
}

phases <- tribble(
  ~round, ~start,                ~end,
  "r1",   as.Date("2023-03-01"), as.Date("2023-04-30"),
  "r2",   as.Date("2023-08-01"), as.Date("2023-10-31"),
  "r3",   as.Date("2023-11-01"), as.Date("2024-01-31")
)

df_ind <-
  df_ind |>
  mutate(
    age_yrs = if_else(age_unit == "Month", as.numeric(age) / 12, as.numeric(age)),
    age_eligible = age_yrs >= 0.5 & age_yrs <= 15,
    vacc_log = coerce_logical(vacc),
    card_chr = as.character(vacc_card),
    card_present = card_chr %in% c("Yes", "Yes, but it is missing/lost", "Was vaccinated in the DTC/OPD/Contact clinic"),
    card_no = card_chr == "No",
    vacc_date = as.Date(date_vacc_card),
    date_in_campaign = in_any_phase(vacc_date, phases),
    vacc_epi_confirmed = case_when(
      card_present & (is.na(date_in_campaign) | !date_in_campaign) ~ TRUE,
      card_present & date_in_campaign ~ NA,
      card_no ~ FALSE,
      vacc_log == FALSE ~ FALSE,
      .default = NA
    ),
    camp_binary_flag = if_else(!is.na(vacc_camp_binary), as.logical(vacc_camp_binary), NA),
    camp_date_flag = date_in_campaign,
    camp_text_flag = if_else(
      !is.na(vacc_camp) & str_detect(as.character(vacc_camp), regex("round|campaign|mass|mop-up", ignore_case = TRUE)),
      TRUE, NA
    ),
    camp_confirmed = case_when(
      (camp_binary_flag == TRUE) | (camp_date_flag == TRUE) | (camp_text_flag == TRUE) ~ TRUE,
      camp_binary_flag == FALSE ~ FALSE,
      .default = NA
    ),
    vacc_source = case_when(
      vacc_epi_confirmed == TRUE & camp_confirmed == TRUE ~ "both",
      vacc_epi_confirmed == TRUE & (camp_confirmed == FALSE | is.na(camp_confirmed)) ~ "epi_only",
      (vacc_epi_confirmed == FALSE | is.na(vacc_epi_confirmed)) & camp_confirmed == TRUE ~ "camp_only",
      vacc_log == FALSE & (camp_confirmed == FALSE | is.na(camp_confirmed)) ~ "neither",
      .default = "unknown"
    ) |>
      factor(levels = c("both", "epi_only", "camp_only", "neither", "unknown"))
  )

# ---- De-identify and select down to what srvy_results.R / srvy_maps_public.R need ----

df_ind_public <-
  df_ind |>
  select(
    index_hh, index_ind, cluster, lga, weight, hh_size,
    age, age_cat, age_grp, age_grp_diph, sex,
    vacc, vacc_pre, vacc_card, vacc_epi_binary, vacc_camp_binary, vacc_ever,
    vacc_epi_confirmed, camp_confirmed, vacc_source, age_eligible,
    diphtheria, diagnosis, died, died_cause, died_diph_binary, died_location,
    date_onset, date_died, pt,
    tx_seeking, tx_seeking_binary, tx_seeking_hsb, delay_tx, delay_tx_binary,
    tx_hf_no_1, tx_hf_no_2, tx_hf_no_3, tx_hf_no_4, tx_hf_no_5, tx_hf_no_99,
    sequelae, ses_quintile, lat, lon
  )

df_hh_public <-
  df_hh |>
  select(index_hh, cluster, lga, consent, hh_size, num_left, num_died)

dir.create(here("code/public_release/data"), showWarnings = FALSE, recursive = TRUE)
write_rds(df_ind_public, here("code/public_release/data/df_ind.rds"))
write_rds(df_hh_public, here("code/public_release/data/df_hh.rds"))
write_rds(
  list(map_nga_ward = map_nga_ward, map_nga_ward_union = map_nga_ward_union),
  here("code/public_release/data/map_boundaries.rds")
)

cat(sprintf(
  "Wrote df_ind (%d rows), df_hh (%d rows), map_boundaries.rds to code/public_release/data/\n",
  nrow(df_ind_public), nrow(df_hh_public)
))
