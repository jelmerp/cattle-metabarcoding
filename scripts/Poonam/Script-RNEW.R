##Settings
n_cores <- 8

##INSTALL AND LOAD PACKAGES-----------------------------------------------------------------------------------------------
if (!require(pacman)) install.packages("pacman")

packages <- c('tidyverse', 'gridExtra', 'dada2',
              'phyloseq', 'DECIPHER', 'phangorn')
pacman::p_load(char = packages)
#devtools::install_version("phangorn", version = "2.7.0")


##SET FILE PATHS-----------------------------------------------------------------------------------------------------------
# Dir with input fastq files:
indir <- 'results/cutadapt'

# Dirs for output:
filter_dir <- 'results/dada/filtered-fastq'
outdir <- 'analysis/ASV_inference'

dir.create(filter_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

stopifnot(indir != filter_dir)

# File with sample metadata:
metadata_file <- 'Metadata.xlsx'


##ASSIGN FASTQ FILES TO FORWARD AND REVERSE READS-----------------------------------------------------------------------
fastqs_raw_F <- sort(list.files(indir, pattern = '_R1_001.fastq.gz', full.names = TRUE))
fastqs_raw_R <- sort(list.files(indir, pattern = '_R2_001.fastq.gz', full.names = TRUE))

print('First fastq files:')
head(fastqs_raw_F)

## CHECK SAMPLE IDS-----------------------------------------------------------------------------------------------------
print() 'Load and prepare sample metadata...'
library(readxl)
metadata_df <- read_xlsx(path = metadata_file) %>% 
  as.data.frame()


colnames(metadata_df)[1] <- 'SampleID'
rownames(metadata_df) <- metadata_df$SampleID

print('First fastq files:')
head(metadata_df)

## Compare Sample ID from the metadata
sampleIDs <- sub("_S.*", "", basename(fastqs_raw_F))
sampleIDs

print('Are the sample IDs from the metadata and the fastq files the same?')
identical(sort(metadata_df$SampleID), sampleIDs)

print('Are any samples missing from the fastq files?')
setdiff(sort(metadata_df$SampleID), sampleIDs)

print('Are any samples missing from the metadata?')
setdiff(sampleIDs, sort(metadata_df$SampleID))

#Remove fastq files that we don't need
fastqs_raw_F <-fastqs_raw_F [match(sort(metadata_df$SampleID), sampleIDs)]
fastqs_raw_R <- fastqs_raw_R[match(sort(metadata_df$SampleID), sampleIDs)]


##fastqs_raw_F <- fastqs_raw_F[1:10]
##fastqs_raw_R <- fastqs_raw_R[1:10]


#FILTERING AND QUALITY TRIMMING--------------------------------------------------------------------------------------------------
fastqs_filt_F <- file.path(filter_dir, basename(fastqs_raw_F))
fastqs_filt_R <- file.path(filter_dir, basename(fastqs_raw_R))


filter_results <-
  filterAndTrim(fastqs_raw_F, fastqs_filt_F,
                fastqs_raw_R, fastqs_filt_R,
                truncLen = c(200, 150), # check amplicon length
                trimLeft = 0,#default
                maxN = 0, #default
                maxEE = Inf, #default
                truncQ = 2, #default
                rm.phix = TRUE, #default
                multithread = n_cores, 
                compress = FALSE, verbose = TRUE) 

head(filter_results)

##DEREPLICATION AND ERROR TRAINING---------------------------------------------------------------------------------------
fastqs_derep_F <- derepFastq(fastqs_filt_F, verbose = FALSE)
fastqs_derep_R <- derepFastq(fastqs_filt_R, verbose = FALSE)

names(fastqs_derep_F) <- sampleIDs
names(fastqs_derep_R) <- sampleIDs

err_F <- learnErrors(fastqs_derep_F, multithread = n_cores, verbose = TRUE)
err_R <- learnErrors(fastqs_derep_R, multithread = n_cores, verbose = TRUE)

print('...Done!')
Sys.time()

plotErrors(err_F, nominalQ = TRUE)
plotErrors(err_R, nominalQ = TRUE)


##Step 5 INFER ASVS-------------------------------------------------------------------------------------------------------
dada_Fs <- dada(fastqs_derep_F, err = err_F, pool = "pseudo", multithread = n_cores)
dada_Rs <- dada(fastqs_derep_R, err = err_R, pool = "pseudo", multithread = n_cores)

dada_Fs[[1]]


##STEP6 MERGE READ PAIRS---------------------------------------------------------------------------------------------------
mergers <- mergePairs(dada_Fs, fastqs_derep_F,
                      dada_Rs, fastqs_derep_R,
                      verbose = TRUE)

saveRDS(fastqs_derep_F, file = file.path(outdir, 'fastqs_derep_F.rds'))
saveRDS(fastqs_derep_R, file = file.path(outdir, 'fastqs_derep_R.rds'))

rm(fastqs_derep_F, fastqs_derep_R) # Remove objects from environment





#STEP7 CONSTRUCT A SEQUENCE TABLE-----------------------------------------------------------------------------------------
seqtab_all <- makeSequenceTable(mergers)

# The dimensions of the object are the nr of samples (rows) and the nr of ASVs (columns):
dim(seqtab_all)


table(nchar(getSequences(seqtab_all)))


# If you need to remove sequences of a particular length (e.g. too long):
# seqtab2 <- seqtab[, nchar(colnames(seqtab_all)) %in% seq(250,256)]


#STEP8 REMOVE CHIMERAS
seqtab <- removeBimeraDenovo(seqtab_all,
                             method = 'pooled',
                             multithread = n_cores,
                             verbose = TRUE)
ncol(seqtab)

saveRDS(seqtab, file = file.path(outdir, 'seqtab_V4.rds'))



##STEP9 GENERATE A SUMMARY TABLE-----------------------------------------------------------------------------------------
getN <- function(x) {
  sum(getUniques(x))
}

denoised_F <- sapply(dada_Fs, getN)
denoised_R <- sapply(dada_Rs, getN)
merged <- sapply(mergers, getN)

nreads_summary <- data.frame(filter_results,
                             denoised_F,
                             denoised_R,
                             merged,
                             nonchim = rowSums(seqtab),
                             row.names = sampleIDs)
colnames(nreads_summary)[1:2] <- c('input', 'filtered')

# Have a look at the first few rows:
head(nreads_summary)

write.table(nreads_summary, file = file.path(outdir, 'nreads_summary.txt'),
            sep = "\t", quote = FALSE, row.names = TRUE)


##STEP10 ASSIGN TAXONOMY TO ASVS
# SET-UP -----------------------------------------------------------------------
## Report
message("\n## Starting script tax_assign_dada.R")
Sys.time()
message()

## Parse command-line arguments
if (!require(argparse)) install.packages("argparse", repos = "https://cran.rstudio.com/")
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-i", "--seqtab",
                    type = "character", required = TRUE,
                    help = "Input file (sequence table RDS) (REQUIRED)")
parser$add_argument("-o", "--taxa",
                    type = "character", required = TRUE,
                    help = "Output file (taxa RDS file) (REQUIRED)")
parser$add_argument("-r", "--ref_url",
                    type = "character",
                    default = "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz",
                    help = "Taxonomic reference URL [default %(default)s]")
parser$add_argument("-s", "--species_url",
                    type = "character",
                    default = "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz",
                    help = "Taxonomic reference URL [default %(default)s]")
parser$add_argument("-c", "--cores",
                    type = "integer", default = 1,
                    help = "Number of cores (threads) to use [default %(default)s]")
args <- parser$parse_args()

seqtab_rds <- args$seqtab
taxa_rds <- args$taxa
n_cores <- args$cores
ref_url <- args$ref_url
species_url <- args$species_url

## Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
packages <- c("BiocManager", "dada2", "DECIPHER", "tidyverse", "argparse")
pacman::p_load(char = packages)

## Constants
TAX_LEVELS <- c("domain", "phylum", "class", "order", "family", "genus", "species")

## Create output dir if needed
outdir <- dirname(taxa_rds)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Define output files
tax_file <- file.path(outdir, basename(ref_url))
species_file <- file.path(outdir, basename(species_url))
qc_file <- file.path(outdir, "tax_prop_assigned_dada.txt")
plot_file <- file.path(outdir, "tax_prop_assigned_dada.png")

## Report
message("## Sequence table RDS file (input):               ", seqtab_rds)
message("## Taxa RDS file (output):                        ", taxa_rds)
message("## Proportion-assigned QC file (output):          ", qc_file)
message("## Number of cores:                               ", n_cores)
message()
message("## Taxonomic assignment file (downloaded input):  ", tax_file)
message("## Species assignment file (downloaded input):    ", species_file)
message()


# FUNCTIONS --------------------------------------------------------------------
## Function to get the proportion of ASVs assigned to taxa
qc_tax <- function(taxa, TAX_LEVELS) {
  prop <- apply(taxa, 2,
                function(x) round(length(which(!is.na(x))) / nrow(taxa), 4))
  
  prop <- data.frame(prop) %>%
    rownames_to_column("tax_level") %>%
    mutate(tax_level = factor(tax_level, levels = TAX_LEVELS))
  
  return(prop)
}


# PREPARE INPUT DATA -----------------------------------------------------------
## Create a DNAStringSet from the ASVs
seqtab <- readRDS(seqtab_rds)
dna <- DNAStringSet(getSequences(seqtab))


# DADA2 TAX. ASSIGNMENT --------------------------------------------------------
## Get and load DADA training set
## (Check for an up-to-date version at <https://benjjneb.github.io/dada2/training.html>)
if (!file.exists(tax_file)) download.file(url = ref_url, destfile = tax_file)
if (!file.exists(species_file)) download.file(url = species_url, destfile = species_file)

## Assign taxonomy
message("\n## Assigning taxonomy...")
taxa <- assignTaxonomy(seqtab, tax_file, multithread = n_cores)

message("\n## Adding species-level assignments...")
taxa <- addSpecies(taxa, species_file)
colnames(taxa) <- TAX_LEVELS

## Save RDS file
saveRDS(taxa, taxa_rds)


# QC TAX. ASSIGNMENTS ----------------------------------------------------------
## Create df with proportion assigned
prop_assigned <- qc_tax(taxa, TAX_LEVELS = TAX_LEVELS)
write_tsv(prop_assigned, qc_file)

message("\n## Prop. ASVs assigned to taxonomy:")
print(prop_assigned)

## Create barplot
p <- ggplot(prop_assigned) +
  geom_col(aes(x = tax_level, y = prop, fill = tax_level),
           color = "grey20") +
  scale_fill_brewer(palette = "Greens") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(y = "Proportion of ASVs assigned", x = NULL) +
  guides(fill = "none") +
  theme_bw(base_size = 14)
ggsave(plot_file, p, width = 7, height = 7)


# WRAP UP ----------------------------------------------------------------------
message("\n## Listing output files:")
system(paste("ls -lh", taxa_rds))
system(paste("ls -lh", qc_file))
system(paste("ls -lh", plot_file))

message("\n## Done with script tax_assign_dada.R")
Sys.time()
message()
