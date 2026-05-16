# Merge AUS and PAK per-location BLUEs into a single long-format CSV used as
# input for phenotype summaries (correlation plots, descriptive statistics).
# These BLUEs collapse across years within a location and are NOT the input
# to the genomic prediction models.
# Run from preprocessing/BLUEs/:  Rscript merge_BLUEs_for_phenotype_plots.R

library(tidyr)

AUS_FILE <- "BLUEs_results/AUS/trait_BLUEs_results.csv"
PAK_FILE <- "BLUEs_results/PAK/trait_BLUEs_results.csv"
OUTPUT_FILE <- "BLUEs_results/AUSPAK_BLUEs_for_phenotype_plots.csv"

aus <- read.csv(AUS_FILE, stringsAsFactors = FALSE)
pak <- read.csv(PAK_FILE, stringsAsFactors = FALSE)

stopifnot(identical(sort(names(aus)), sort(names(pak))))

combined_long <- rbind(aus, pak[, names(aus)])

combined_long$BLUE <- round(combined_long$BLUE, 2)

combined <- pivot_wider(
  combined_long[, c("location", "accession", "SampleName", "trait", "BLUE")],
  names_from = trait,
  values_from = BLUE
)

message(sprintf("AUS rows (long): %d", nrow(aus)))
message(sprintf("PAK rows (long): %d", nrow(pak)))
message(sprintf("combined rows (wide): %d", nrow(combined)))

write.csv(combined, file = OUTPUT_FILE, row.names = FALSE)
