# Diferential expresion analysis: "Transcriptional and physiological adaptations in nucleus accumbens somatostatin interneurons that regulate behavioral responses to cocaine"

**Course project:** *Intro RNA-seq (LCG-UNAM 2026)* — Leonardo Collado-Torres  

**Objective:** Perform a reproducible differential gene expression analysis using a public study available through **recount3**, including at least 3 figures and a biological interpretation of the results.

## Study Information
- recount3 project ID: SRP151726
- Original study title: "Transcriptional and physiological adaptations in nucleus accumbens somatostatin interneurons that regulate behavioral responses to cocaine"
- Organism: Mus musculus
- Experimental design:
    Saline (control)
    Cocaine (1 hour post-injection)
    10 samples per group

The data were accessed through the **recount3** resource, which provides processed RNA-seq gene-level expression summaries.

## Analysis Overview
The workflow includes:

1. Downloading gene-level counts from recount3
2. Metadata preparation and QC
3. Gene filtering (filterByExpr)
4. TMM normalization
5. Mean–variance modeling (voom)
6. Linear modeling (~ assigned_gene_prop + group)
7. Differential expression testing


## Repository structure
This repository follows the folder structure from `LieberInstitute/template_project`.

```text
RNA_seq_project/
├── code/
│   ├── 01_download_recount3/
│   ├── 02_prepare_metadata/
│   └── 03_differential_expression/
├── output/
│   ├── rse/
│   └── tables/
├── plots/
├── Reports/
│   └── Report_reanalysis.md
├── README.md
└── template_project/
```

## Report
The full analysis and interpretation are available here:
[Report_results.md](Reports/Report_results.md)

Key results:
[Top 50 nominal genes](output/tables/top50_nominal_p.csv)
[Plots directory](plots/)

## Reproducibility
To reproduce the analysis:
1. Open the project in RStudio.
2. Run scripts sequentially:
  - 01_download_recount3
  - 02_prepare_metadata
  - 03_differential_expression

All tables are saved under [output/](output/) and all figures under [plots/](plots/)

## References
This project follows the workflow developed by Dr. Leonardo Collado-Torres and uses the recount3 resource for gene-level RNA-seq summaries.

**Software and Resources**
> Collado-Torres L (2023). Explore and download data from the recount3 project.Bioconductor package recount3, version 1.8.0. https://doi.org/10.18129/B9.bioc.recount3 https://github.com/LieberInstitute/recount3

> Wilks C, Zheng SC, Chen FY, et al. (2021).recount3: summaries and queries for large-scale RNA-seq expression and splicing. Genome Biology, 22, 323. https://doi.org/10.1186/s13059-021-02533-6

**Original Study**

> Ribeiro EA, Salery M, Scarpa JR, et al. (2018). Transcriptional and physiological adaptations in nucleus accumbens somatostatin interneurons that regulate behavioral responses to cocaine. Nature Communications, 9, 3149.https://doi.org/10.1038/s41467-018-05657-9