# Load libraries
library(data.table)  # Fast VCF file reading
library(AGHmatrix)   # Kinship matrix calculation

# Read VCF file, skipping header lines until the column names line (#CHROM) - file can be downloaded
vcf <- fread("quinoa_551accessions_genomic_prediction.vcf", skip="#CHROM", sep="\t")

# Extract sample IDs from columns 10 onwards (first 9 columns are VCF metadata)
samples <- colnames(vcf)[10:ncol(vcf)]  
# Extract genotype calls into a matrix (all columns after the first 9 VCF standard columns)
snp_matrix <- as.matrix(vcf[, ..samples]) 

### recode genotypes to numeric format
# convert diploid genotype calls to 0,1,2
# handle unphased (/) and phased (|) genotypes
snp_matrix[snp_matrix %in% c("0/0", "0|0")] <- 0
snp_matrix[snp_matrix %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
snp_matrix[snp_matrix %in% c("1/1", "1|1")] <- 2

# rename column one in vcf to "CHROM" (it was read as "#CHROM")
colnames(vcf)[1] <- "CHROM"

# Convert genotype matrix from character to numeric
snp_matrix <- matrix(as.numeric(snp_matrix), nrow=nrow(snp_matrix))

# create unique SNP identifiers as "chromosome_position"
rownames(snp_matrix) <- paste0(vcf$CHROM, "_", vcf$POS)
colnames(snp_matrix) <- samples

### Create kinship matrix using VanRaden method
# Transpose matrix as needed for AGHmatrix where samples are rows and SNPs are columns
snp_matrix <- t(snp_matrix)

# Calculate additive relationship matrix with VanRaden method 
kinship_matrix_V2 <- Gmatrix(snp_matrix, 
                            method = "VanRaden", 
                            ploidy = 2,
                            thresh.missing = 0.2)     # our max missing 20% threshold instead of 50%

# save the kinship matrix as RData for future use
save(kinship_matrix_V2, file = "kinship_matrix_VanRaden_auspak_maxmissing20.RData")