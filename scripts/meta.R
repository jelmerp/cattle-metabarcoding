## Install and load packages
if (!require(pacman)) install.packages("pacman")
packages <- c('tidyverse', "readxl")
pacman::p_load(char = packages)

## Process the metadata
meta <- readxl::read_xlsx("metadata/Metadata.xlsx") %>%
  rename(calf_id = `Calf ID`,
         truck_nr = Group,
         sampling_time = Time,
         travel_time = txgp) %>%
  mutate(neg_control = ifelse(truck_nr == "200uL sterile H2O", TRUE, FALSE),
         mock_comm = ifelse(calf_id == "Mock Community", TRUE, FALSE),
         calf_id = sub("CONTROL 1", NA, calf_id),
         calf_id = sub("Mock Community", NA, calf_id),
         truck_nr = sub("200uL sterile H2O", NA, truck_nr)) %>% 
  as.data.frame()
colnames(meta)[1] <- 'sampleID'
rownames(meta) <- meta$sampleID

## Write output file
write_tsv(meta, "metadata/metadata.txt")
