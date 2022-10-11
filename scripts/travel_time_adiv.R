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

## Estimate alpha diversity
div_df <- estimate_richness(CS, measures = "Shannon")%>%
  merge(., meta(CS) , by = "row.names") %>%
  drop_na()

## Compute abundance differences relative to sampling_time 0:
div_diff_df <- div_df %>%
  select(calf_id, truck_nr, sampling_time, travel_time, Shannon) %>%
  ## Get abundances at different sampling_points side-by-side
  pivot_wider(id_cols = c(calf_id, truck_nr, travel_time),
              names_from = sampling_time, values_from = Shannon,
              names_prefix = "div_at_") %>%
  ## Get rid of rows with missing abundances
  filter(!if_any(starts_with("div_at_"), is.na)) %>%
  ## Put abundances at later timepoints (SAT, 24, 72) into 1 column
  pivot_longer(cols = c(div_at_SAT, div_at_24, div_at_72),
               names_to = "sampling_time",
               values_to = "div_at_later_sampling") %>% 
  ## Get abundance diff. between sampling_time 0 and each of the later sampling
  mutate(div_diff = div_at_0 - div_at_later_sampling) %>%
  ## Get mean abundance diff. across calfs:
  #group_by(travel_time, sampling_time) %>%
  #summarize(div_diff = mean(div_diff)) %>%
  mutate(sampling_time = factor(sub("div_at_", "", sampling_time),
                                levels = c("SAT", 24, 72)))

## Create the plot:
ggplot(div_diff_df) +
  aes(x = sampling_time, y = div_diff, fill = travel_time) +
  geom_hline(yintercept = 0, color = "darkred") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2),
              color = "grey30") +
  scale_fill_brewer(palette = "Dark2") +
  labs(fill = "Travel time",
       x = "Sampling time",
       y = "Difference in diversity\nrelative to Sampling Time 0") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("results/fig/samplingtime_div_diff.png", width = 8, height = 6)
