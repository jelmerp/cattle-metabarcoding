#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --output=slurm-dada-%j.out
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8

## Settings
n_cores <- 8
n_subset <- 4 # Set to 'NA' to disable subsetting

## Install and load packages
if (!require(pacman)) install.packages("pacman")
packages <- c('tidyverse', 'gridExtra', 'dada2',
              'phyloseq', 'DECIPHER', 'phangorn',
              "readxl")
pacman::p_load(char = packages)

## Set directories
indir <- "results/cutadapt"

filter_dir <- "results/dada/fq_filtered"
outdir <- 'results/dada'

dir.create(filter_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

stopifnot(indir != filter_dir)

## Load the metadata
metadata_df <- readxl::read_xlsx("metadata/Metadata.xlsx") %>%
  as.data.frame()
colnames(metadata_df)[1] <- 'SampleID'
rownames(metadata_df) <- metadata_df$SampleID

## FASTQ files
fastqs_raw_F <- sort(list.files(indir, pattern = '_R1_001.fastq.gz', full.names = TRUE))
fastqs_raw_R <- sort(list.files(indir, pattern = '_R2_001.fastq.gz', full.names = TRUE))

sampleIDs <- sub("_S.*", "", basename(fastqs_raw_F))

print('Are the sample IDs from the metadata and the fastq files the same?')
identical(sort(metadata_df$SampleID), sampleIDs)

print('Are any samples missing from the fastq files?')
setdiff(sort(metadata_df$SampleID), sampleIDs)

print('Are any samples missing from the metadata?')
setdiff(sampleIDs, sort(metadata_df$SampleID))

## Remove FASTQ filenames for samples that aren't in the metadata
fastqs_raw_F <- fastqs_raw_F[match(sort(metadata_df$SampleID), sampleIDs)]
fastqs_raw_R <- fastqs_raw_R[match(sort(metadata_df$SampleID), sampleIDs)]

## Subset FASTQ files
if (!is.na(n_subset)) {
  message("\n## NOTE: Taking only the first ", n_subset, " samples!\n")
  fastqs_raw_F <- fastqs_raw_F[1:n_subset]
  fastqs_raw_R <- fastqs_raw_R[1:n_subset]
  sampleIDs <- sub("_S.*", "", basename(fastqs_raw_F))
}

# FILTERING AND TRIMMING -------------------------------------------------------
fastqs_filt_F <- file.path(filter_dir, basename(fastqs_raw_F))
fastqs_filt_R <- file.path(filter_dir, basename(fastqs_raw_R))

filter_results <-
  filterAndTrim(fastqs_raw_F, fastqs_filt_F,
                fastqs_raw_R, fastqs_filt_R,
                truncLen = c(200, 150),
                trimLeft = 0, # default
                maxN = 0, # default
                maxEE = Inf, # default
                truncQ = 2, # default
                rm.phix = TRUE, # default
                multithread = n_cores, 
                compress = FALSE, verbose = TRUE)

# DEREPLICATION AND ERROR TRAINING
fq_derep_f <- derepFastq(fastqs_filt_F, verbose = FALSE)
fq_derep_r <- derepFastq(fastqs_filt_R, verbose = FALSE)

names(fq_derep_f) <- sampleIDs
names(fq_derep_r) <- sampleIDs

err_f <- learnErrors(fq_derep_f, multithread = n_cores)
err_r <- learnErrors(fq_derep_r, multithread = n_cores)

# ASV INFERENCE
dada_Fs <- dada(fq_derep_f, err = err_f, pool = "pseudo", multithread = n_cores)
saveRDS(dada_Fs, "results/dada/dada_F.rds")

# REMOVE CHIMERAS
# seqtab <- removeBimeraDenovo(seqtab_all,
#                              method = 'pooled',
#                              multithread = n_cores,
#                              verbose = TRUE)


message("Done with dada script")
