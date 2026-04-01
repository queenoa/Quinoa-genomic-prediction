# Load library
library(SNPRelate)  # VCF to GDS conversion and PCA

# if .gds file is not available yet
### convert VCF to GDS for PCA
# SNPRelate requires filepath to a VCF file - file can be downloaded
# vcf.fn <- "quinoa_551accessions_genomic_prediction.vcf"


# # convert VCF to GDS format with error handling
# tryCatch({
#   snpgdsVCF2GDS(vcf.fn, "quinoa_auspak.gds", 
#                 method="copy.num.of.ref",
#                 verbose=TRUE)
# }, error = function(e) {
#   print(paste("Error occurred:", e$message))
# })

# open gds file for reading and PCA
pca_genofile<- snpgdsOpen("quinoa_AUSPAK.gds")
# check summary statistics
# number of samples, SNPs, chromosomes
snpgdsSummary(pca_genofile)


##### perform PCA

pca_quinoa <- snpgdsPCA(pca_genofile, 
                        eigen.cnt = 0,        # add this because the default would be only to return 32 eigenvectors
                        autosome.only=FALSE,  # include all chromosomes, not just autosomes and otherwise SNPRelate defaults to human autosome filtering
                        num.thread=4,        # use 4 CPU threads for parallel processing 
                        verbose=TRUE)         # show progress

### check variance explained by each PC
# check variance proportion (%)
pc.percent <- pca_quinoa$varprop*100
round(pc.percent, 2)

### extract all PCs
# Automatically determine number of PCs (all available)
n_pcs <- ncol(pca_quinoa$eigenvect)
# create data frame with sample IDs and all PCs
quinoa_PCs <- data.frame(sample.id = pca_quinoa$sample.id,
                  setNames(as.data.frame(pca_quinoa$eigenvect[,1:n_pcs]), 
                          paste0("PC", 1:n_pcs)),
                  stringsAsFactors = FALSE)

# save pca_quinoa object for later use in R
save(pca_quinoa, file = "pca_quinoa.RData")
# save PCs to CSV
write.csv(quinoa_PCs, "AUSPAK_PCs_all.csv", row.names = FALSE)
