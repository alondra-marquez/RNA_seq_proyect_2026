library(recount3)
library(SummarizedExperiment)

# Load raw RSE object downloaded in previous step
rse <- readRDS("output/rse/rse_SRP151726.rds")

# Parse sra.sample_attributes into sra_attribute.* columns
rse <- expand_sra_attributes(rse)

# Clean treatment text
rse$sra_attribute.treatment <- factor(tolower(rse$sra_attribute.treatment))

# Create DE group with saline as baseline

rse$group <- ifelse(
  grepl("cocaine", rse$sra_attribute.treatment),
  "cocaine",
  "saline"
)
rse$group <- factor(rse$group, levels = c("saline", "cocaine"))

stopifnot(all(table(rse$group) == c(10, 10)))

# Add QC covariates to colData
meta <- as.data.frame(colData(rse))

# QC metric: proportion assigned to genes
colData(rse)$assigned_gene_prop <- meta$recount_qc.gene_fc_count_all.assigned /
  meta$recount_qc.gene_fc_count_all.total

# Library size from counts assay
colData(rse)$lib_size <- colSums(assay(rse, "counts"))
colData(rse)$log_lib_size <- log10(colData(rse)$lib_size + 1)

# Save prepared object for DE
saveRDS(rse, "output/rse/rse_prepared.rds")
