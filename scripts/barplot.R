## Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
packages <- c("BiocManager", "dada2", "phyloseq", "DECIPHER", "tidyverse", 
              "microbiome", "plotly", "ampvis2")
pacman::p_load(char = packages)
##remotes::install_github("MadsAlbertsen/ampvis2")

## Settings
theme_set(theme_bw())
mycols <- c("#a74bb4", "#62b54f", "#7064d3", "#b5b348", "#dd6fc5",
            "#4db18a", "#ce416e", "#45aecf", "#d55035", "#7784cb",
            "#cc8e3e", "#ac6090", "#647934", "#df837e", "#9c5736")

## Run the script with functions for a grouped barplot
source("scripts/barplot_funs.R")

## Prep the phyloseq object
CS <- readRDS("results/dada/ps_filt.rds")
CS <- subset_samples(CS, !is.na(truck_nr))
CS <- prune_taxa(taxa_sums(CS) > 0, CS)
sample_data(CS)$sampling_time <- factor(sample_data(CS)$sampling_time,
                                        levels = c(0, 6, 12, 16, 24, 48, 72))
sample_data(CS)$travel_time <- factor(sample_data(CS)$travel_time,
                                      levels = c(6, 12, 16))
sample_data(CS)$truck_nr <- factor(sample_data(CS)$truck_nr,
                                   levels = seq(1, 7))
CS_prop <- microbiome::transform(CS, "compositional")


# PLOT BY TRAVEL TIME ----------------------------------------------------------
## Remove samples with NA for travel_time
CS_traveltime <- subset_samples(CS_prop, !is.na(travel_time))

## Barplot for phylum level
fCS <- tax_glom(CS_traveltime, taxrank = "phylum")
pbar(ps = fCS, taxrank = "phylum",
     xvar = "travel_time", facetvar = "sampling_time",
     abund_tres = 0.001, cols = mycols, xlab = "Travel time")

pbar(ps = fCS, taxrank = "phylum",
     xvar = "truck_nr", facetvar = "travel_time",
     abund_tres = 0.001, cols = mycols,
     xlab = "Truck number")

## Barplot for family level
fCS <- tax_glom(CS_traveltime, taxrank = "family")
pbar(ps = fCS, taxrank = "family",
     xvar = "travel_time", facetvar = "truck_nr",
     abund_tres = 0.02, cols = mycols,
     xlab = "Travel time")


# PLOT BY AGE ------------------------------------------------------------------
## Remove samples with NA for Age
CS_Age <- subset_samples(CS_prop, !is.na(Age))

## Barplot for phylum level
fCS <- tax_glom(CS_Age, taxrank = "phylum")
pbar(ps = fCS, taxrank = "phylum",
     xvar = "Age", #facetvar = "truck_nr",
     abund_tres = 0.001, cols = mycols, xlab = "Age")

## Barplot for family level
fCS <- tax_glom(CS_Age, taxrank = "family")
pbar(ps = fCS, taxrank = "phylum",
     xvar = "Age", #facetvar = "truck_nr",
     abund_tres = 0.001, cols = mycols, xlab = "Age")

