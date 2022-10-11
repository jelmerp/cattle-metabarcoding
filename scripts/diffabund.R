## Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
packages <- c("BiocManager", "phyloseq", "tidyverse", "microbiome",
              "DESeq2")
pacman::p_load(char = packages)

## Load input files
ps <- readRDS("results/dada/ps_filt_poonam.rds")

## Get count matrix
counts <- t(otu_table(ps)) %>% as(., "matrix")
## Prevalence thresholds - features that do not pass these thresholds are excluded from for statistical testing
## 10% threshold from https://www.biorxiv.org/content/10.1101/2021.05.10.443486v1.full
PREV_FRAC <- 0.05        # Fraction of samples the ASV/pathway should be present in at prevalence `PREV_N`
PREV_N <- 1              # Min count for the ASV/pathway for at least `PREV_FRAC` of the samples
prev_fun <- function(x) sum(x > PREV_N) > (PREV_FRAC * length(x))
counts <- counts[apply(counts, 1, prev_fun), ]

## Make DESeq object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta(ps),
                              design = ~1)

## Analysis by age
## http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
dds_age <- dds
design(dds_age) <- formula("~ Age")
dds_age <- subset(dds_age, select = !is.na(Age))
dds_age <- DESeq(dds_age)
results(dds_age, contrast = x)
