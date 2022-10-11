## Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
packages <- c("BiocManager", "dada2", "phyloseq", "DECIPHER", "tidyverse")
pacman::p_load(char = packages)

## Input files
seqtab_rds <- "analysis/ASV_inference/seqtab_V4.rds"
taxa_rds <- "taxa_rds"
tree_rds <- "results/dada/tree.rds"
meta_file <- "meta_copy.csv"

## Output file
ps_rds <- "results/dada/ps_raw.rds"


# LOAD INPUT DATA --------------------------------------------------------------
seqtab <- readRDS(seqtab_rds)
taxa <- readRDS(taxa_rds)
tree <- readRDS(tree_rds)

meta <- read.table(meta_file, sep = ",", header = TRUE, row.names = 1)
write.csv(meta, "meta1.csv")


# CHECK SAMPLE IDs -------------------------------------------------------------
message("\n## Are the sample IDs from the metadata and the seqtab the same?")
identical(sort(meta$sampleID), rownames(seqtab))

message("## Are any samples missing from the seqtab?")
setdiff(sort(meta$sampleID), rownames(seqtab))

message("## Are any samples missing from the metadata?")
setdiff(rownames(seqtab), sort(meta$sampleID))


# CREATE PHYLOSEQ OBJECT -------------------------------------------------------
ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE),
               phy_tree(tree$tree),
               sample_data(meta),
               tax_table(taxa))


# RENAME ASVs AND STORE SEQS SEPARATELY ----------------------------------------
## Extract ASV sequences (which are the taxa_names)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)

## Merge sequence object into the phyloseq object:
ps <- merge_phyloseq(ps, dna)

## Rename ASVs from full sequences to ASV1...ASVx
taxa_names(ps) <- paste("ASV", 1:ntaxa(ps), sep = "_")


# SAVE AND WRAP UP -------------------------------------------------------------
## Save phyloseq object in an RDS file
saveRDS(ps, ps_rds)


###### ps_filter

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


#######Phyloseq plot


## Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
packages <- c("BiocManager", "dada2", "phyloseq", "DECIPHER", "tidyverse", 
              "microbiome", "plotly", "ampvis2")
pacman::p_load(char = packages)
#install.packages("ggplot2", type="source")
##remotes::install_github("MadsAlbertsen/ampvis2")
##library("plyr"); packageVersion("plyr")
theme_set(theme_bw())
##install.packages("https://cran.r-project.org/src/contrib/Archive/rlang/rlang_0.4.10.tar.gz", repos = NULL, type="source")
library("ggplot2")
library("plotly")
library("ampvis2")

##install.packages("rlang")
theme_set(theme_bw())


## Prep the phyloseq object
CS <- readRDS("results/dada/ps_filt.rds")
CS <- subset_samples(CS, !is.na(truck_nr))
CS <- prune_taxa(taxa_sums(CS) > 0, CS)
sample_data(CS)$sampling_time <- factor(sample_data(CS)$sampling_time,
                                        levels = c(0, AT, 24, 48, 72))
sample_data(CS)$travel_time <- factor(sample_data(CS)$travel_time,
                                      levels = c(6, 12, 16))

## NMDS plot
CS.ord <- ordinate(CS, "NMDS", "bray")
plot_df <- plot_ordination(CS, CS.ord,
                           type = "samples", color = "travel_time",
                           title = "taxa", justDF = TRUE)


plot_df <- na.omit(plot_df)
p1 <- ggplot(data = plot_df,
             mapping = aes(x = NMDS1, y = NMDS2, color = travel_time)) +
  facet_wrap(vars(sampling_time)) +
  geom_point(aes(text = paste("Sample ID:", SampleID))) 

ggplotly(p1)
meta
