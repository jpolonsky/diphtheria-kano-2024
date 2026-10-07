# diphtheria-kano-2024

Analysis code for two companion studies of the 2023–24 diphtheria epidemic response in Kano State, Nigeria (MSF Operational Centre West and Central Africa, in collaboration with the Kano State Ministry of Health).

- **Home-based care (HBC) study** (retrospective matched cohort comparing home-based care with facility-based treatment for mild diphtheria): `code/hbc/hbc_results.R`. Manuscript under revision at PLOS Medicine; preprint: [medrxiv.org/content/10.64898/2026.04.10.26350586v1](https://www.medrxiv.org/content/10.64898/2026.04.10.26350586v1)
- **Household survey study** (population-based survey of community diphtheria burden): `code/srvy/srvy_results.R`. Preprint: [medrxiv.org/content/10.64898/2026.04.10.26348327v2](https://www.medrxiv.org/content/10.64898/2026.04.10.26348327v2)

`srvy_results.R` reproduces every number, table and figure reported in the survey manuscript. `hbc_results.R` reproduces the results reported in the HBC manuscript from the released input files, except for quantities that depend on information not in those files: the telephone tracing outcomes and reasons for non-interview, the link between replacement patients and the patients they replaced, the interval from notification to interview, the breakdown of sequelae, complication and comorbidity types, Figure 1 (the flow diagram), and the Total column of Table 1. S1 Table is printed rather than written to a file, and the balance, effective sample size and trimming results behind S4 Table are printed rather than tabulated.

## Structure

```
code/
  hbc/
    hbc_results.R        analysis script
    hbc_prep.R    builds the 3 rds files hbc_results.R reads, from the study team's internal data store
  srvy/
    srvy_results.R       tables, non-spatial results, and spatial figures
    srvy_prep.R   builds the rds files srvy_results.R reads, from the study team's internal data store
data/                     not included in this repo (see Data access below)
output/
  figures/                written by the scripts when run
  tables/
```

## Running

Each results script (`hbc_results.R`, `srvy_results.R`) expects to be run from the repository root, with R (version 4.3 or later) and the packages listed at the top of the script installed, and its required input files in place under `data/` (see Data access below). Running a script writes all figures and tables under `output/`.

- `hbc_results.R` reads `data/df_pt.rds`, `data/df_cohort_samp.rds`, `data/df_hh.rds`.
- `srvy_results.R` reads `data/df_ind.rds`, `data/df_hh.rds`, `data/map_boundaries.rds`.

`hbc_prep.R` and `srvy_prep.R` are included for transparency: they document exactly how those files are derived from the study team's internal data store, but they read from paths internal to that data store and will not run outside it.

## Data access

Data used in this study are subject to restrictions arising from participant informed consent and GDPR-aligned privacy protections for identifying and potentially stigmatising health information, and cannot be deposited in a public repository. Anonymised individual-level data are available to qualified researchers upon request to the [Epicentre Data Sharing Committee](https://epicentre.msf.org/en/our-projects/study-data-access-request), which reviews all requests independently of the study authors in accordance with [Epicentre's Data Sharing Policy](https://epicentre.msf.org/sites/default/files/2020-02/Epicentre-Data%20Sharing%20Policy-EN-2020.pdf).

## Citation

If you use this code, please cite the corresponding manuscript(s) (see the preprint links above; citation details will be updated on publication) and, for the code itself:

Jonathan Polonsky. (2026). jpolonsky/diphtheria-kano-2024 [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.21907933
