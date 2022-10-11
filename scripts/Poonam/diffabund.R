## Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
packages <- c("BiocManager",  "phyloseq", "tidyverse", 
              "microbiome", "DESeq2")
pacman::p_load(char = packages)
##BiocManager::install("DESeq2")
##library("DESeq2")
##install.packages("locfit")

##Load input files
ps <- readRDS("results/dada/ps_filt.rds")
#ps <- tax_glom(ps, taxrank = "phylum")

tax_df <- tax_table(ps) %>% as.data.frame()


## Get count matrix
counts <- t(otu_table(ps)) %>% as(., "matrix")
## Prevalence thresholds - features that do not pass these thresholds are excluded from for statistical testing
## 10% threshold from https://www.biorxiv.org/content/10.1101/2021.05.10.443486v1.full
PREV_FRAC <- 0.05        # Fraction of samples the ASV/pathway should be present in at prevalence `PREV_N`
PREV_N <- 1              # Min count for the ASV/pathway for at least `PREV_FRAC` of the samples
prev_fun <- function(x) sum(x > PREV_N) > (PREV_FRAC * length(x))

counts <- counts[apply(counts, 1, prev_fun), ]

## Make DESeq object
meta_df <- meta(ps)%>% 
  mutate(truck_nr = factor(truck_nr))
##head(meta_df)
##class(meta_df$calf_id)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta_df,
                              design = ~1) ## can change Age

## Analysis by age
## http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
dds_age <- dds
design(dds_age) <- formula("~ truck_nr + Age")
design(dds_age) <- formula("~ Age")
dds_age <- subset(dds_age, select = !is.na(Age))
dds_age <- DESeq(dds_age)
res_age <- results(dds_age, name = "Age") %>%
  as.data.frame() %>%
  arrange(padj)

sum(res_age$padj <0.05, na.rm =TRUE)
res_sig_age <- res_age %>% filter(padj < 0.05)
resultsNames(dds_age)
#colData(dds_age)
#?results

plot_counts <- function(ASV_ID, dds_age, tax_df) {
  d <- plotCounts(dds, gene=ASV_ID,
                  intgroup = c("Age", "travel_time"), 
                  returnData=TRUE)
  tax_row <- tax_df[rownames(tax_df) == ASV_ID, ]
  tax_info <- paste(tax_row$family, ":", tax_row$genus)
  
  ggplot(d, aes(x=Age, y=count, color = travel_time)) + 
    geom_point(position=position_jitter(w=0.1,h=0)) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = ASV_ID, subtitle = tax_info) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  #ggsave()
}

plot_counts(ASV_ID = rownames(res_age)[1], dds = dds_age, tax_df = tax_df)
map(rownames(res_sig_age), plot_counts, dds = dds_age, tax_df = tax_df)




###travel time 6 vs 16
dds_traveltime <- dds
design(dds_traveltime) <- formula("~ travel_time")
dds_traveltime$travel_time <- factor(dds_traveltime$travel_time, levels = c("6","12", "16")) # suppose 3 groups
dds_traveltime <- subset(dds_traveltime, select = !is.na(travel_time))
dds_traveltime <- DESeq(dds_traveltime)
resultsNames(dds_traveltime)
res_traveltime1 <- results(dds_traveltime, contrast=c("travel_time","6","16")) %>%
  as.data.frame() %>%
  arrange(padj)
sum(res_traveltime$padj <0.05, na.rm =TRUE)
res_sig_traveltime1 <- res_traveltime %>% filter(padj < 0.05)

##Plots
sigtab = cbind(as(res_sig_traveltime1, "data.frame"), as((tax_df)[rownames(res_sig_traveltime1), ], "matrix"))

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$phylum = factor(as.character(sigtab$phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$genus = factor(as.character(sigtab$genus), levels=names(x))
ggplot(sigtab, aes(x=genus, y=log2FoldChange, color=phylum)) + geom_point(size=2) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


plot_counts(ASV_ID = "ASV_160", dds = dds_traveltime, tax_df = tax_df)
sigtab %>% filter(genus == "Fusobacterium")


###travel time 6 vs 12
res_traveltime2 <- results(dds_traveltime, contrast=c("travel_time","6","12")) %>%
  as.data.frame() %>%
  arrange(padj)
sum(res_traveltime2$padj <0.05, na.rm =TRUE)
res_sig_traveltime2 <- res_traveltime2 %>% filter(padj < 0.05)

##Plots
sigtab1 = cbind(as(res_sig_traveltime2, "data.frame"), as((tax_df)[rownames(res_sig_traveltime2), ], "matrix"))

# Phylum order
x = tapply(sigtab1$log2FoldChange, sigtab1$phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab1$phylum = factor(as.character(sigtab1$phylum), levels=names(x))
# Genus order
x = tapply(sigtab1$log2FoldChange, sigtab1$genus, function(x) max(x))
x = sort(x, TRUE)
sigtab1$genus = factor(as.character(sigtab1$genus), levels=names(x))
ggplot(sigtab1, aes(x=genus, y=log2FoldChange, color=phylum)) + geom_point(size=2) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


plot_counts(ASV_ID = "ASV_151", dds = dds_traveltime1, tax_df = tax_df)
sigtab %>% filter(genus == "Fusobacterium")

###travel time 12 vs 16
res_traveltime3 <- results(dds_traveltime, contrast=c("travel_time","12","16")) %>%
  as.data.frame() %>%
  arrange(padj)
sum(res_traveltime3$padj <0.05, na.rm =TRUE)
res_sig_traveltime3 <- res_traveltime3 %>% filter(padj < 0.05)

##Plots
sigtab2 = cbind(as(res_sig_traveltime3, "data.frame"), as((tax_df)[rownames(res_sig_traveltime3), ], "matrix"))

# Phylum order
x = tapply(sigtab2$log2FoldChange, sigtab2$phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab2$phylum = factor(as.character(sigtab2$phylum), levels=names(x))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$genus = factor(as.character(sigtab2$genus), levels=names(x))
ggplot(sigtab2, aes(x=genus, y=log2FoldChange, color=phylum)) + geom_point(size=2) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


plot_counts(ASV_ID = "ASV_151", dds = dds_traveltime1, tax_df = tax_df)
sigtab %>% filter(genus == "Fusobacterium")




plot_counts <- function(ASV_ID, dds, tax_df) {
  d <- plotCounts(dds, gene=ASV_ID,
                  intgroup = c("travel_time"), 
                  returnData=TRUE)
  tax_row <- tax_df[rownames(tax_df) == ASV_ID, ]
  tax_info <- paste(tax_row$family, ":", tax_row$genus)
  
  ggplot(d, aes(x=travel_time, y=count, color = travel_time)) + 
    geom_point(position=position_jitter(w=0.1,h=0)) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = ASV_ID, subtitle = tax_info) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  #ggsave()
}

plot_counts(ASV_ID = rownames(res_sig_traveltime1)[1], dds = dds_traveltime, tax_df = tax_df)
map(rownames(res_sig_traveltime1), plot_counts, dds = dds_traveltime, tax_df = tax_df)

## http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
dds_traveltime <- dds
dds$travel_time <- factor(dds$travel_time, levels = c("6","12", "16"))
dds_traveltime <- DESeq(dds_traveltime)
res_travel <- results(dds_traveltime)
res_travel
plot_counts(ASV_ID = rownames(res_travel)[1], dds = dds_traveltime, tax_df = tax_df)


