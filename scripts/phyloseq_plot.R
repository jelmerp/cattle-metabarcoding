## Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
packages <- c("BiocManager", "dada2", "phyloseq", "DECIPHER",
              "tidyverse", "microbiome", "plotly", "ampvis2",
              "agricolae")
pacman::p_load(char = packages)

## Settings
theme_set(theme_bw())

## Prep the phyloseq object
CS <- readRDS("results/dada/ps_filt_poonam.rds")
CS <- subset_samples(CS, !is.na(truck_nr))
CS <- prune_taxa(taxa_sums(CS) > 0, CS)
sample_data(CS)$sampling_time <- factor(sample_data(CS)$sampling_time,
                                        levels = c(0, 6, 12, 16, 24, 48, 72))
sample_data(CS)$travel_time <- factor(sample_data(CS)$travel_time,
                                      levels = c(6, 12, 16))
CS_rar <- rarefy_even_depth(CS, sample.size = min(sample_sums(CS)))

## NMDS plot
CS.ord <- ordinate(CS, "NMDS", "bray")
plot_df <- plot_ordination(CS, CS.ord,
                           type = "samples", color = "travel_time",
                           title = "taxa", justDF = TRUE)

p1 <- ggplot(data = plot_df,
             mapping = aes(x = NMDS1, y = NMDS2, color = travel_time)) +
  facet_wrap(vars(Age)) +
  geom_point(aes(text = paste("Sample ID:", SampleID))) 

ggplotly(p1)

## Plotting alpha diversity
#plot_richness(CS_rar, measures=c("Chao1", "Shannon"))
div_df <- estimate_richness(CS_rar, measures = "Shannon")%>%
  merge(., meta(CS_rar) , by = "row.names")

ggplot(data = div_df,
       mapping = aes(x = sampling_time, y = Shannon)) +
  facet_wrap(vars(Age)) +
  geom_boxplot() +
  geom_point(aes(color = travel_time))

## Make a rarefaction plot
asv_table <- data.frame(t(otu_table(CS_rar)), tax_table(CS_rar))
amp <- amp_load(asv_table, meta(CS_rar))
amp_rarecurve(amp, color_by = "travel_time") +
  facet_wrap(vars(Age))
amp_rarecurve(amp, color_by = "sampling_time") +
  facet_wrap(vars(Age))  

## Testing alpha diversity
set.seed(100)
kruskal.test(Shannon ~ travel_time, data = div_df)
kruskal.test(Shannon ~ sampling_time, data = div_df)
kruskal.test(Shannon ~ Age, data = div_df)

## Kruskal-wallis test
div_df_samp6 <- div_df %>% filter(sampling_time == 6)
#kruskal.test(Shannon ~ Age, data = div_df_samp6)
kruskal_res <- kruskal(y = div_df_samp6$Shannon, trt = div_df_samp6$Age,
                       p.adj = "hochberg")
kruskal_res


## Plotting the tree
# ntaxa(CS)
# physeq <- prune_taxa(taxa_names(CS)[1:50], CS)
# plot_tree(physeq)
# plot_tree(physeq, "treeonly")
# plot_tree(physeq, nodelabf=nodeplotboot(), ladderize="left", color="travel_time")
