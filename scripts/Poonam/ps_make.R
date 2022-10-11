## Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
packages <- c("BiocManager", "dada2", "phyloseq", "DECIPHER", "tidyverse")
pacman::p_load(char = packages)

## Input files
seqtab_rds <- "analysis/ASV_inference/seqtab_V4.rds"
taxa_rds <- "taxa_rds"
tree_rds <- "results/dada/tree.rds"
meta_file <- "metadata/metadata.csv"

## Output file
ps_rds <- "results/dada/ps_raw.rds"


# LOAD INPUT DATA --------------------------------------------------------------
seqtab <- readRDS(seqtab_rds)
taxa <- readRDS(taxa_rds)
tree <- readRDS(tree_rds)

meta <- read.table(meta_file, sep = ",", header = TRUE, row.names = 1)
write.csv(meta, "meta.csv")


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
