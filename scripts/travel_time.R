##TODO: Why do some calf_ids have multiple travel_times?

## Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
packages <- c("BiocManager",  "phyloseq", "tidyverse", 
              "microbiome", "DESeq2")
pacman::p_load(char = packages)

## Prep the phyloseq object
CS <- readRDS("results/dada/ps_filt.rds")
CS <- subset_samples(CS, !is.na(truck_nr))     # NOT SURE IF THIS SHOULD BE EXCLUDED?
CS <- subset_samples(CS, sampling_time != 48)
CS <- prune_taxa(taxa_sums(CS) > 0, CS)
sample_data(CS)$travel_time <-
  factor(sample_data(CS)$travel_time, levels = c(6, 12, 16))

## Lump sampling_times 6, 12, and 16 into 'SAT'
sample_data(CS)$sampling_time <-
  ifelse(sample_data(CS)$sampling_time %in% c(6, 12, 16),
         "SAT", sample_data(CS)$sampling_time)
sample_data(CS)$sampling_time <-
  factor(sample_data(CS)$sampling_time, levels = c(0, "SAT", 24, 72))

## Convert proportional abundance, agglomerate by phylum
CS_prop <- microbiome::transform(CS, "compositional")
CS_prop_phylum <- tax_glom(CS_prop, taxrank = "phylum")

## "Melt" the object to get all counts in long format
psm <- psmelt(CS_prop_phylum)

## Compute abundance differences relative to sampling_time 0:
ps_diff <- psm %>%
  select(OTU, phylum, calf_id, truck_nr, sampling_time, travel_time, Abundance) %>%
  ## Get abundances at different sampling_points side-by-side
  pivot_wider(id_cols = c(OTU, phylum, calf_id, truck_nr, travel_time),
              names_from = sampling_time, values_from = Abundance,
              names_prefix = "abundance_at_") %>%
  ## Get rid of rows with missing abundances
  filter(!if_any(starts_with("abundance_at_"), is.na)) %>%
  ## Put abundances at later timepoints (SAT, 24, 72) into 1 column
  pivot_longer(cols = c(abundance_at_SAT, abundance_at_24, abundance_at_72),
               names_to = "sampling_time",
               values_to = "abundance_at_later_sampling") %>% 
  ## Get abundance diff. between sampling_time 0 and each of the later sampling
  mutate(abund_diff = abundance_at_0 - abundance_at_later_sampling) %>%
  ## Get mean abundance diff. across calfs:
  group_by(OTU, phylum, travel_time, sampling_time) %>%
  summarize(abund_diff = mean(abund_diff)) %>%
  mutate(sampling_time = factor(sub("abundance_at_", "", sampling_time),
                                levels = c("SAT", 24, 72)))

## Filter out phyla with no big(ish) changes in abundance:
ps_diff <- ps_diff %>%
  group_by(OTU) %>%
  filter(max(abs(abund_diff)) > 0.01)

## Create the plot:
ggplot(ps_diff) +
  aes(x = phylum, y = abund_diff, fill = travel_time) +
  geom_col(position = position_dodge(preserve = "single"),
           color = "grey10") +
  facet_wrap(vars(sampling_time)) +
  scale_fill_brewer(palette = "Dark2") +
  labs(fill = "Travel time",
       y = "Difference in proprotional abundance\nrelative to Sampling Time 0") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("results/fig/samplingtime_abundance_diff.png", width = 8, height = 6)
