## Settings
ASV_THRESHOLD <- 1000

## Packages
packages <- c("tidyverse", "phyloseq", "decontam")
pacman::p_load(char = packages)

## Define input file
ps_infile <- "results/dada/ps_raw.rds"

## Define output files
ps_outfile <- "results/dada/ps_filt.rds"
contam_df_file <- "results/dada/contanimants.txt"

## Load the phyloseq object
ps_raw <- readRDS(ps_infile)


# CONTAMINATION DETECTION & REMOVAL --------------------------------------------
## Run the contamination detection function
contam_df <- isContaminant(ps_raw, method = "prevalence", neg = "neg_control")
table(contam_df$contaminant)
ps_noncontam <- prune_taxa(!contam_df$contaminant, ps_raw)

## Check contaminants
contam_df$abund <- taxa_sums(ps_raw)
contam_tax <- tax_table(prune_taxa(contam_df$contaminant, ps_raw))
contam_df <- contam_df %>%
  filter(contaminant == TRUE) %>%
  merge(., contam_tax, by = "row.names")
colnames(contam_df)[1] <- "ASV" 
write_tsv(contam_df, contam_df_file)

## This % of the total ASV count was removed:
(sum(contam_df$abund) / sum(sample_sums(ps_raw))) * 100
## This % of unique ASVs was removed
(nrow(contam_df) / ntaxa(ps_raw)) * 100


# OFF-TARGET TAXA REMOVAL ------------------------------------------------------
ps_taxremoval <- ps_noncontam

if (any(ps_taxremoval@tax_table[, "order"] == "Chloroplast", na.rm = TRUE)) {
  #ps_chloroplast <- subset_taxa(ps_taxremoval, order == "Chloroplast")
  ps_taxremoval <- subset_taxa(ps_taxremoval, order != "Chloroplast")
}
if (any(ps_taxremoval@tax_table[, "family"] == "Mitochondria", na.rm = TRUE)) {
  ps_taxremoval <- subset_taxa(ps_taxremoval, family != "Mitochondria")
}
if (any(ps_taxremoval@tax_table[, "domain"] == "Eukaryota", na.rm = TRUE)) {
  ps_taxremoval <- subset_taxa(ps_taxremoval, domain != "Eukaryota")
}

## Number of taxa removed:
ntaxa(ps_noncontam) - ntaxa(ps_taxremoval)


# FILTER SAMPLES ---------------------------------------------------------------
## Remove samples with too few ASVs
ps <- subset_samples(ps_taxremoval, sample_sums(ps_taxremoval) > ASV_THRESHOLD)
## Number of samples removed:
nsamples(ps) - nsamples(ps_taxremoval)

## Remove negative control samples
ps <- subset_samples(ps, neg_control != TRUE | is.na(neg_control))

## Remove any ASVs that now don't occur in the dataset 
ps <- subset_taxa(ps, taxa_sums(ps) > 0)


saveRDS(ps, ps_outfile)
