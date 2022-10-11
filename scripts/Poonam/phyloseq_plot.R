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
                                        levels = c(0, 6, 12, 16, 24, 48, 72))
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
  facet_wrap(vars(Age)) +
  geom_point(aes(text = paste("Sample ID:", SampleID))) 

ggplotly(p1)

CS.ord <- ordinate(CS, "NMDS", "bray")
plot_df <- plot_ordination(CS, CS.ord,
                           type = "samples", color = "travel_time",
                           title = "taxa", justDF = TRUE)
plot_df <- na.omit(plot_df)
p2 <- ggplot(data = plot_df,
             mapping = aes(x = NMDS1, y = NMDS2, color = travel_time)) +
  geom_point(aes(text = paste("Sample ID:", SampleID))) 

ggplotly(p2)

### NMDS plot with changes
CS.ord <- ordinate(CS, "NMDS", "bray")
plot_df <- plot_ordination(CS, CS.ord,
                           type = "samples", color = "Age",
                           title = "taxa", justDF = TRUE)
 
p3 <- ggplot(data = plot_df,
             mapping = aes(x = NMDS1, y = NMDS2, color = Age)) +
  facet_wrap(vars(travel_time)) +
  geom_point(aes(text = paste("Sample ID:", SampleID)))
ggplotly(p3)

## Plotting alpha diversity - 
plot_richness(CS, measures=c("Chao1", "Shannon", "Simpson"))

## Plotting alpha diversity - Shannon
##Sampling time - Age
div_df <- estimate_richness(CS, measures = "Shannon")%>%
  merge(., meta(CS) , by = "row.names")
div_df <- na.omit(div_df)
a1 <- ggplot(data=div_df,
       mapping = aes(x = Age, y = Shannon)) +
  facet_wrap(vars(Age)) +
  geom_boxplot() +
  geom_point(aes(color = travel_time))

a1 
a1 <- ggplot(data=div_df,
             mapping = aes(x = Age, y = Shannon)) +
  geom_boxplot() +
  geom_point(aes(color = travel_time))

a1 

##Sampling time 
plot_richness(CS, measures=c("Chao1", "Shannon"))
div_df <- estimate_richness(CS, measures = "Shannon")%>%
  merge(., meta(CS) , by = "row.names")
div_df <- na.omit(div_df)
a2 <- ggplot(data=div_df,
             mapping = aes(x = sampling_time, y = Shannon)) +
  geom_boxplot() +
  geom_point(aes(color = travel_time))

a2 
##Travel time - Age
a3 <- ggplot(data=div_df,
             mapping = aes(x = travel_time, y = Shannon)) +
  facet_wrap(vars(Age)) +
  geom_boxplot() +
  geom_point(aes(color = sampling_time))

a3

##Travel time 
a4 <- ggplot(data=div_df,
             mapping = aes(x = travel_time, y = Shannon)) +
  geom_boxplot() +
  geom_point(aes(color = sampling_time))

a4

## Plotting alpha diversity - Simpson
plot_richness(CS, measures=c("Chao1", "Shannon", "Simpson"))
div_df <- estimate_richness(CS, measures = "Simpson")%>%
  merge(., meta(CS) , by = "row.names")
div_df <- na.omit(div_df)
## Travel time
a5 <- ggplot(data=div_df,
             mapping = aes(x = travel_time, y = Simpson)) +
  geom_boxplot() +
  geom_point(aes(color = sampling_time))

a5 
## Sampling time
a6 <- ggplot(data=div_df,
             mapping = aes(x = sampling_time, y = Simpson)) +
  geom_boxplot() +
  geom_point(aes(color = travel_time))

a6 
## Age
a7 <- ggplot(data=div_df,
             mapping = aes(x = Age, y = Simpson)) +
  geom_boxplot() +
  geom_point(aes(color = travel_time))

a7 

## Plotting alpha diversity - Chao1
plot_richness(CS, measures=c("Chao1", "Shannon", "Simpson"))
div_df <- estimate_richness(CS, measures = "Chao1")%>%
  merge(., meta(CS) , by = "row.names")
div_df <- na.omit(div_df)
## Travel time
a8 <- ggplot(data=div_df,
             mapping = aes(x = travel_time, y = Chao1)) +
  geom_boxplot() +
  geom_point(aes(color = sampling_time))

a8 
## Sampling time
a9 <- ggplot(data=div_df,
             mapping = aes(x = sampling_time, y = Chao1)) +
  geom_boxplot() +
  geom_point(aes(color = travel_time))

a9 
## Age
a10 <- ggplot(data=div_df,
             mapping = aes(x = Age, y = Chao1)) +
  geom_boxplot() +
  geom_point(aes(color = travel_time))

a10 

## Make a rarefaction plot
library(ampvis2)
asv_table <- data.frame(t(otu_table(CS)), tax_table(CS))
amp <- amp_load(asv_table, meta(CS))
amp_rarecurve(amp, color_by = "travel_time") +
  facet_wrap(vars(Age))
amp_rarecurve(amp, color_by = "sampling_time") +
  facet_wrap(vars(travel_time))  
amp_rarecurve(amp, color_by = "Age") +
  facet_wrap(vars(travel_time))  

## Testing alpha diversity
set.seed(100)
CS_rar <- rarefy_even_depth(CS, sample.size = min(sample_sums(CS))) ## rarefying data
kruskal.test(Shannon ~ travel_time, data = div_df)
kruskal.test(Shannon ~ sampling_time, data = div_df)
kruskal.test(Shannon ~ Age, data = div_df)

div_df_samp6 <- div_df %>% filter(sampling_time == 6)
#kruskal.test(Shannon ~ Age, data = div_df_samp6)
kruskal(y = div_df_samp6$Shannon, trt = div_df_samp6$Age, 
        p.adj = "hochberg")

kruskal_res





## Plotting the tree
# ntaxa(CS)
# physeq <- prune_taxa(taxa_names(CS)[1:50], CS)
# plot_tree(physeq)
# plot_tree(physeq, "treeonly")
# plot_tree(physeq, nodelabf=nodeplotboot(), ladderize="left", color="travel_time")

##install.packages("agricolae")
library("agricolae")
##install.packages("haven")
