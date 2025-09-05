################################################################################
## Script Name:        Summarise Seabird Data
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-09-04
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)
memory.limit(30000000)

############################### DEPENDENCIES ###################################
# Load required libraries
library(tidyverse)
library(magrittr)
source('./scripts/format_seabird_orders.R')

################################### DATA #######################################
# Read and inspect data
# Lu seabird data
seabird_taxa <- read_csv('./data/seabirds_allinOne_June_2025_AVONET.csv') %>%
  rename_with(., ~ str_to_lower(gsub(" |\\.", "_", .x))) 

# ebird taxonomy
ebird_taxa <- read_csv('./data/eBird_taxonomy_v2024.csv')

# H5 sequence metadata
sequence_meta <- read_csv('./data/gisaid_info_tbl-wg_infoGEO.2.wg_info_full.csv')

# Birdlife taxonomy
birdlife_taxa <- readxl::read_xlsx('~/Downloads/zip_file/HBW_BirdLife_List of Birds_v.9.xlsx') 


################################### MAIN #######################################
# Main analysis or transformation steps
seabird_orders_abrv <- c('sph', 
                         'pha', 
                         'sul',
                         'cha',
                         'pro')

birdlife_taxa <- birdlife_taxa %>%
  mutate(order_name = str_to_title(Order)) %>% 
  rename_with(., ~ str_to_lower(gsub(" ", "_", .x, fixed = TRUE))) %>%
  dplyr::select(ends_with('name')) %>%
  relocate(order_name, family_name, scientific_name, common_name) 

# Make generic names for ambiguous sequences
family_only_names <- ebird_taxa %>% 
  filter(CATEGORY == 'spuh') %>% 
  select(scientific_name = SCI_NAME,
         common_name = PRIMARY_COM_NAME,
         family_name = FAMILY,
         order_name = ORDER) %>%
  mutate(family_name = str_split_i(family_name, ' \\(', 1))


# Format sequence data
sequence_meta_with_full_taxa <- sequence_meta %>%
  # Filter out non-seabird/non-avian 
  filter(host == 'Avian') %>%
  filter(birdOrder %in% seabird_orders_abrv) %>%
  
  # Format remaining to get BirdLife scientific names
  rowwise() %>%
  mutate(Phaethontiformes = FormatPhaethontiformes(bird),
         Charadriiformes = FormatCharadriiformes(bird),
         Procellariiformes = FormatProcellariiformes(bird),
         Sphenisciformes = FormatSphenisciformes(bird),
         Suliformes = FormatSuliformes(bird)) %>%
  as_tibble()  %>%
  
  mutate(scientific_name = coalesce(Phaethontiformes,
                                    Charadriiformes,
                                    Procellariiformes, 
                                    Sphenisciformes,
                                    Suliformes)) %>%
  drop_na(scientific_name) %>%
  
  # Join BirdLife Binomial Nomenclature
  left_join(birdlife_taxa, by = 'scientific_name') %>%
  rows_patch(family_only_names, by = 'scientific_name', unmatched = 'ignore') %>%
  dplyr::select(-c(ends_with('formes'),
                   bird,
                   birdOrder))


# Summarise sequences across seabird orders/families
order_total <- sequence_meta_with_full_taxa %>%
  count(order_name)

family_total <- sequence_meta_with_full_taxa %>%
  count(family_name)

order_by_year <- sequence_meta_with_full_taxa %>%
  count(order_name, year)

family_by_year <- sequence_meta_with_full_taxa %>%
  count(family_name, year)

# Summary plots


################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################