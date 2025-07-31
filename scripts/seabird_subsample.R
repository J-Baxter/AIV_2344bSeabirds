################################################################################
## Script Name:        <INSERT_SCRIPT_NAME_HERE>
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-04-02
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6, # Avoid scientific notation
  digits = 7 # Set precision for numerical display
)
memory.limit(30000000)

############################### DEPENDENCIES ###################################
# Load required libraries
library(tidyverse)
library(magrittr)
library(ape)


# User functions
MakeTipNames <- function(data) {
  out <- data %>%
    unite(tipnames,
      virus_subtype,
      clade,
      isolate_id,
      host_order,
      collection_countryname,
      cluster_profile,
      collection_tipdate,
      sep = "|",
      remove = FALSE
    ) %>%
    mutate(tipnames = gsub(" ", "_", tipnames))

  return(out)
}


ReNameAlignment2 <- function(alignment, data) {
  z <- names(alignment)

  isolates <- regmatches(z, gregexpr("EPI_ISL_(china_){0,1}\\d+[^.|]*", z)) %>% unlist()

  new_seqnames <- sapply(isolates, function(x) data$new_tipnames[data$isolate_id %in% x]) %>%
    as.vector()

  names(alignment) <- new_seqnames

  return(alignment)
}


GroupSequences <- function(aln, snp_threshold = 0) {
  require(phangorn)
  require(igraph)

  # Ensure alignment is correctly formatted
  if (class(aln) != "PhyDat") {
    aln_formatted <- as.phyDat(aln)
  } else {
    aln_formatted <- aln
  }

  # Calculate hamming distance
  hd_normalised <- dist.hamming(aln_formatted) %>%
    as.matrix()
  hd_raw <- hd_normalised * ncol(aln)

  # Obtain groups of sequences for which HD < SNP threshold

  groups <- which(hd_raw <= snp_threshold,
    arr.ind = TRUE
  ) %>%
    dplyr::as_tibble(rownames = "tipnames") %>%
    filter(row != col) %>%
    dplyr::select(-tipnames) %>%
    # Infer network from HDs
    igraph::graph_from_data_frame(.,
      directed = F
    ) %>%
    components() %>%
    getElement("membership") %>%
    stack() %>%
    as_tibble() %>%
    mutate(ind = as.numeric(as.character(ind))) %>%
    mutate(tipnames = map_chr(ind, ~ rownames(aln)[.x])) %>%
    dplyr::select(c(tipnames, values)) %>%
    dplyr::distinct() %>%
    dplyr::rename(sequence_group = values)

  out <- tibble(tipnames = rownames(aln)) %>%
    left_join(groups) %>%
    mutate(
      sequence_group =
        ifelse(is.na(sequence_group),
          max(sequence_group, na.rm = T) + row_number() + 1,
          sequence_group
        )
    )


  return(out)
}


################################### DATA #######################################
# Read and inspect data
all_meta <- read_csv("./data/2024-09-09_meta.csv")

new_seabirds <- read_csv("./data/seabirds_tropicbirds_add_20250625.csv")

alignment_files <- c(
  list.files("./2024Aug18/region_alignments",
    pattern = "ha_",
    full.names = T
  ),
  list.files("./2024Aug18/reassortant_alignments",
    pattern = "ha_",
    full.names = T
  )
)


alignments <- lapply(alignment_files, read.dna, format = "fasta", as.matrix = F)

################################### MAIN #######################################
# Main analysis or transformation steps
new_seabirds_order <- new_seabirds %>%
  select(host_order = order) %>%
  distinct() %>%
  mutate(host_order = str_to_lower(host_order)) %>%
  mutate(is_seabird_order = 1)

new_seabirds_family <- new_seabirds %>%
  select(host_family = family) %>%
  distinct() %>%
  mutate(host_family = str_to_lower(host_family)) %>%
  mutate(is_seabird_family = 1)

new_seabirds_sciname <- new_seabirds %>%
  select(host_sciname = Scientific.name) %>%
  distinct() %>%
  mutate(host_sciname = str_to_lower(host_sciname)) %>%
  mutate(is_seabird_sciname = 1)

seabirds <- read_csv("./data/seabird_H5Nx_meta.csv") %>%
  mutate(is_seabird = 1) %>%
  select(tipnames, is_seabird)

temp <- do.call(c, alignments)
names(temp) <- gsub(".*\\.", "", names(temp))
temp <- temp[!duplicated(names(temp))]


groups <- GroupSequences(as.matrix(temp), snp_threshold = 1)

groups_with_seabirds <- groups %>%
  left_join(meta %>% select(tipnames, host_family)) %>%
  left_join(new_seabirds_family) %>%
  group_by(sequence_group) %>%
  # Replace non seabirds wth 0
  replace_na(list(is_seabird_family = 0)) %>%
  # only include groups with at least one seabird
  filter(sum(is_seabird_family) > 0) %>%
  ungroup()


ha_seabird_subsample <- all_meta %>%
  drop_na(collection_date) %>%
  filter(grepl("H5", virus_subtype)) %>%
  filter(isolate_id %in% str_extract(names(temp), "EPI_ISL_(china_){0,1}\\d+[^.|]*")) %>%
  filter(clade == "2344b") %>%
  select(
    tipnames,
    isolate_id,
    collection_date,
    collection_datemonth, 
    host_family, 
    host_order, 
    collection_regionname
  ) %>%
  left_join(groups_with_seabirds) %>%
  drop_na(sequence_group) %>%
  drop_na(host_family) %>%
  # drop_na(cluster_profile) %>%

  # mutate(new_host = case_when(is_seabird_family == 1 ~ 'seabird',
  # .default = host_family)) %>%
  mutate(host_cat = case_when(is_seabird_family == 1 ~ host_family,
    .default = host_order
  )) %>%
  mutate(collection_regionname = case_when(grepl("europe", collection_regionname) ~ "europe",
    grepl("africa", collection_regionname) ~ "africa",
    grepl("asia", collection_regionname) ~ "asia",
    grepl("(central|northern) america|caribbean", collection_regionname) ~ "central & northern america",
    grepl("south america|southern ocean", collection_regionname) ~ "south america",
    grepl("australia|melanesia", collection_regionname) ~ "australasia",
    .default = collection_regionname
  )) %>%
  group_by(is_seabird_family, collection_regionname, collection_datemonth, host_cat) %>%
  slice_sample(n = 1)


# select only genomes that are included
ha_selected_genomes <- temp[str_extract(names(temp), "EPI_ISL_(china_){0,1}\\d+[^.|]*") %in% ha_seabird_subsample$isolate_id]
ha_selected_genomes <- ha_selected_genomes[!duplicated(names(ha_selected_genomes))]


################################### OUTPUT #####################################
# Save output files, plots, or results
write.FASTA(ha_selected_genomes, "./data/alignments/ha_seabird_subsample.fasta")



write_delim(
  ha_seabird_subsample %>%
    ungroup() %>%
    # dplyr::select(tipnames, collection_date, collection_datemonth) %>%
    rename(name = tipnames) %>%
    mutate(
      collection_datemonth = paste(collection_datemonth, "-XX", sep = ""),
      date = coalesce(as.character(collection_date), collection_datemonth)
    ) %>%
    dplyr::select(name, date),
  delim = "\t",
  quote = "needed",
  "./data/alignments/ha_seabird_subsample_dates.tsv"
)


#################################### END #######################################
################################################################################