## Explore recount3 projects (study selection)
## Course: Intro RNA-seq (LCG-UNAM 2026)

library(recount3)

## Use AWS S3 mirror (more stable for downloads)
options(
  recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release"
)

## 1) List available projects for a given organism
mouse_projects <- available_projects(organism = "mouse")

## Save table for transparency
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(
  mouse_projects,
  file = "output/tables/recount3_mouse_projects.csv",
  row.names = FALSE
)

## 2) Filter for a candidate study of interest
project_info <- subset(mouse_projects, project == "SRP151726")

## Check to avoid errors in downstream steps
stopifnot(nrow(project_info) == 1)

## 3) Download and construct RangedSummarizedExperiment
rse_gene <- create_rse(project_info)

## 4) Add "counts" assay for downstream DE workflows
assay(rse_gene, "counts") <- compute_read_counts(rse_gene)

## 5) Save object as .rds for reproducibility and future use
dir.create("output/rse", recursive = TRUE, showWarnings = FALSE)
saveRDS(rse_gene, file = "output/rse/rse_SRP151726.rds")

## Show key sample metadata columns
colnames(colData(rse_gene))
