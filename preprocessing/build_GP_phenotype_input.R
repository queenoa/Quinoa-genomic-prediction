# Build the phenotype input file for the genomic prediction models.
#
# One row per (location, year, accession) in wide format — one column per
# trait. Raw per-accession means are used uniformly across all trials:
#
#   - AUS 2017/2018/2019: raw per-accession means from data/AUS_phenotypes_raw.csv
#   - PAK 2019-20, PAK 2020-21, PAK 2021-22: raw per-accession means from
#                  data/PAK_phenotypes_raw.csv
#
# Using means uniformly across all trials is justified because:
#   - Unreplicated trials: mean = single observation
#   - PAK 2021-22 (2 replicates): BLUEs equal raw means for fully replicated
#     accessions (241/284), and the adjustment for partially replicated
#     accessions (43/284) is negligible given only 2 replicates.
#
# Only sequenced accessions (non-NA SampleName) are kept.
#
# Run from preprocessing/:  Rscript build_GP_phenotype_input.R

suppressPackageStartupMessages(library(dplyr))

TRAITS <- c("DTF", "DTH", "PtHt", "PcleLng", "SdLen", "TGW", "SdW_z")
AUS_RAW <- "../data/AUS_phenotypes_raw.csv"
PAK_RAW <- "../data/PAK_phenotypes_raw.csv"
OUTPUT_FILE <- "../data/AUSPAK_phenotypes_GP_input.csv"

aus_raw <- read.csv(AUS_RAW, na.strings = c("", "NA"), stringsAsFactors = FALSE)
pak_raw <- read.csv(PAK_RAW, na.strings = c("", "NA"), stringsAsFactors = FALSE)

aus_means <- aus_raw %>%
  filter(!is.na(SampleName)) %>%
  mutate(location = "AUS", year = as.character(year)) %>%
  group_by(location, year, accession, SampleName) %>%
  summarise(across(all_of(TRAITS), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(across(all_of(TRAITS), ~ ifelse(is.nan(.x), NA_real_, .x)))

pak_means <- pak_raw %>%
  filter(!is.na(SampleName)) %>%
  mutate(
    location = "PAK",
    year = case_when(
      year == "2019-20" ~ "2019",
      year == "2020-21" ~ "2020",
      year == "2021-22" ~ "2021",
      TRUE ~ as.character(year)
    )
  ) %>%
  group_by(location, year, accession, SampleName) %>%
  summarise(across(all_of(TRAITS), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(across(all_of(TRAITS), ~ ifelse(is.nan(.x), NA_real_, .x)))

combined <- bind_rows(aus_means, pak_means) %>%
  arrange(location, year, accession) %>%
  mutate(
    across(all_of(TRAITS), ~ round(.x, 2)),
    location_year = paste(location, year, sep = "_")
  ) %>%
  select(location_year, everything()) %>%
  rename(sample.id = SampleName)

message(sprintf("AUS raw means rows: %d", nrow(aus_means)))
message(sprintf("PAK raw means rows: %d", nrow(pak_means)))
message(sprintf("combined rows:      %d", nrow(combined)))

write.csv(combined, file = OUTPUT_FILE, row.names = FALSE)