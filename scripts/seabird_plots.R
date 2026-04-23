################################################################################
## Script Name:        <INSERT_SCRIPT_NAME_HERE>
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-07-25
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
library(fitdistrplus)
library(scales)
library(rnaturalearth)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
FormatPhyloGeo <- function(mcc_file, posterior_file) {
  require(rnaturalearthdata)
  require(tidytree)
  require(sf)
  require(tidyverse)
  require(seraphim)
  require(treeio)

  # Import tree data
  mcc_tree <- read.beast(mcc_file)
  tree_tbl <- as_tibble(mcc_tree)

  # Guess most_recent_date
  most_recent_date <- tree_tbl %>%
    slice_min(height_median, n = 1, with_ties = FALSE) %>%
    pull(label) %>%
    str_extract(., "\\d{4}-\\d{2}-\\d{2}") %>%
    ymd() %>%
    decimal_date() %>%
    unique()


  # Extract TMRCA
  start_date <- tree_tbl %>%
    slice_max(height_median, n = 1, with_ties = FALSE) %>%
    pull(height_median) %>%
    as.numeric() %>%
    subtract(most_recent_date, .)
  print(start_date)

  # Scan posterior tree file
  allTrees <- scan(
    file = posterior_file,
    what = "",
    sep = "\n",
    quiet = T
  )


  localTreesDirectory <- "./phylo/temp_tree_dir/"
  do.call(file.remove, list(list.files("./phylo/temp_tree_dir/", full.names = TRUE)))

  burnIn <- 0
  randomSampling <- TRUE
  nberOfTreesToSample <- 100
  # mostRecentSamplingDatum <- most_recent_date
  coordinateAttributeName <- "location"
  treeExtractions(localTreesDirectory,
    allTrees,
    burnIn,
    randomSampling,
    nberOfTreesToSample,
    most_recent_date,
    coordinateAttributeName,
    nberOfCores = 8
  )

  # Step 4: Estimating the HPD region for each time slice ----

  # Format Polygons
  polygons <- suppressWarnings(spreadGraphic2(localTreesDirectory,
    nberOfExtractionFiles = nberOfTreesToSample,
    prob = 0.95,
    startDatum = start_date,
    precision = 0.08
  ))

  # Conver polgons to GEOMETRY, set coordinate system and add values
  polygons_sf <- lapply(polygons, st_as_sf) %>%
    bind_rows() %>%
    pivot_longer(cols = starts_with("2"), names_to = "year", values_to = "value") %>%
    filter(value == 1) %>%
    mutate(year = as.numeric(year)) %>%
    dplyr::select(-value) %>%
    st_set_crs(4326)


  # Format Nodes
  nodes_sf <- tree_tbl %>%
    dplyr::select(node, height_median, location1, location2, host_simplifiedhost, label) %>%
    mutate(height_median = as.numeric(height_median)) %>%
    replace_na(list(height_median = 0)) %>%
    mutate(year = most_recent_date - height_median) %>%
    # Convert to POINT & set coordinate system
    st_as_sf(
      coords = c("location2", "location1"),
      crs = 4326
    )

  # Format Arrows
  edges <- tree_tbl %>%
    dplyr::select(node, parent) %>%
    left_join(tree_tbl %>%
      dplyr::select(node, location1, location2)) %>%
    left_join(
      tree_tbl %>%
        dplyr::select(node, location1, location2),
      by = join_by(parent == node)
    ) %>%
    mutate(across(starts_with("location"), .fns = ~ as.numeric(.x))) %>%
    mutate(
      location1.y = case_when(location1.y == location1.x ~ location1.y + 0.00000001,
        .default = location1.y
      ),
      location2.y = case_when(location2.y == location2.x ~ location2.y + 0.00000001,
        .default = location2.y
      )
    ) %>%
    rename(
      start_lat = location1.x,
      start_lon = location2.x,
      end_lat = location1.y,
      end_lon = location2.y
    )

  edges_sf <- edges %>%
    rowid_to_column(var = "id") %>%
    pivot_longer(starts_with(c("start", "end")),
      names_to = c("type", ".value"),
      names_sep = "_"
    ) %>%
    # Convert coordinate data to sf POINT
    st_as_sf(
      coords = c("lon", "lat"),
      crs = 4326
    ) %>%
    # Convert POINT geometry to MULTIPOINT, then LINESTRING
    group_by(id) %>%
    summarise(do_union = FALSE) %>%
    st_cast("LINESTRING") %>%
    # Convert rhumb lines to great circles
    st_segmentize(units::set_units(20, km)) %>%
    # Wrap dateline correctly
    st_wrap_dateline()


  out <- list(
    polygons = polygons_sf,
    edges = edges_sf,
    nodes = nodes_sf
  )


  return(out)
}



################################### DATA #######################################
# Read and inspect data
seabird_outbreaks <- read_csv("./data/seabird_outbreak_combined_Sciname_updated_eventID_27Jan2026.csv") %>%
  mutate(across(ends_with("date"), .fns = ~ as_date(.x))) %>%
  mutate(observation.date = coalesce(observation.date, report.date)) %>%
  rename(observe_date = observation.date, 
         report_date = report.date) %>%
  drop_na(observe_date) 

seabird_outbreaks_old <- read_csv("./data/seabird_outbreak_16Dec2024.csv") %>%
  mutate(across(ends_with("date"), .fns = ~ as_date(format(as.Date(.x, '%d/%m/%Y'), '%Y-%m-%d')))) %>%
  mutate(observe_date = coalesce(observe_date, report_date)) %>%
  drop_na(observe_date) 

seabird_outbreaks_all <- seabird_outbreaks %>%
  bind_rows(seabird_outbreaks_old) %>%
  distinct(observe_date, Locality, Country, Region)

new_seabirds <- read_csv("./data/seabirds_tropicbirds_add_20250625.csv")

meta <- read_csv('./data/2026-04-20_meta.csv')
################################### MAIN #######################################
read_csv("./data/seabird_outbreak_combined_Sciname_updated_eventID_27Jan2026.csv") 
read_csv("./data/seabird_outbreak_16Dec2024.csv") 

# Main analysis or transformation steps
# 1a.
plt_1a <- seabird_outbreaks_all %>%
  dplyr::select(starts_with("observe")) %>%
  count(year(observe_date)) %>%
  ggplot(aes(x = `year(observe_date)`, y = n)) +
  geom_bar(stat = "identity", fill = "grey75") +
  geom_smooth(
    method = "glm",
    data = seabird_outbreaks_all %>%
      dplyr::select(starts_with("observe")) %>%
      count(year(observe_date)) %>%
      filter(`year(observe_date)` > 2015 & `year(observe_date)` < 2024),
    aes(
      x = `year(observe_date)`,
      y = n
    ),
    inherit.aes = F,
    formula = y ~ x,
    # se = F,
    method.args = list(family = gaussian(link = "log")),
    colour = "darkblue",
    fill = "blue",
    alpha = 0.2
  ) +
  scale_x_continuous("Year",
    expand = c(0, 0),
    limits = c(2005, 2026),
    breaks = pretty_breaks(n = 10)
  ) +
  scale_y_continuous("Number of Outbreak Events",
    limits = c(0, NA),
    expand = c(0, 0)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # Rotate x-axis labels
    axis.title.x = element_text(vjust = 0.5) # Adjust the x-axis title position if needed
  )


# 1b.
region_colours <- c(
  "europe" = "#1b9e77",
  "asia" = "#d95f02",
  "africa" = "#7570b3",
  "antarctica" = "#e7298a",
  "central & northern america" = "#66a61e",
  "south america" = "#e6ab02"
)

plt_1b <- seabird_outbreaks_all %>%
  mutate(new_region = case_when(Country %in% c("Canada", "Greenland", "United States of America", "Panama", "United States") ~ "central & northern america",
    Country %in% c("Argentina", "Ecuador", "Brazil", "Uruguay", "Chile", "Peru" , "Falkland Islands (Islas Malvinas)", "Falkland Islands (Malvinas)") ~ "south america",
    Country %in% c('Crozet Islands', "Prince Edward Islands" , "Kerguelen Islands") ~ "antarctica",
    .default = str_to_lower(Region)
  )) %>% 
  count(new_region, year(observe_date)) %>%
  ggplot(aes(
    x = `year(observe_date)`,
    y = n,
    colour = new_region
  )) +
  geom_line(linewidth = 1) +
  #scale_colour_manual(
    #values = region_colours,
    #"Region",
   # labels = str_to_title
  #) +
  scale_colour_brewer(palette = 'Blues', 
                      'Region',
                       labels = str_to_title) + 
  theme_minimal() +
  scale_x_continuous("Year",
                     limits = c(2005, 2026),
    expand = c(0, 0),
    breaks = pretty_breaks(n = 10)
  ) +
  scale_y_continuous("Number of Outbreak Events",
    expand = c(0, 0),
    limits = c(0, NA)
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(vjust = 0.5),
    legend.position = "inside",
    legend.position.inside = c(0, 1),
    legend.justification = c(0, 1)
  )

# 1c - proportion of seabird sequences (all birds vs seabird families)
seabird_meta <- meta %>%
  left_join(new_seabirds %>%
              dplyr::select(host_family = family) %>%
              mutate(host_family = str_to_lower(host_family)) %>%
              mutate(is_seabird = 1) %>%
              distinct()) %>%
  dplyr::select(
    host_family,
    collection_date,
    collection_datemonth,
    is_seabird
  ) %>%
  drop_na(collection_date) %>%
  mutate(collection_date = format(as.Date(collection_date, "%d/%m/%Y"), format = '%Y-%m-%d')) %>%
  filter(collection_date > as.Date("2017-01-01")) %>%
  replace_na(list(is_seabird = 0)) %>%
  mutate(fill = ifelse(is_seabird == 0, NA, host_family)) %>%
  count(collection_datemonth, fill) %>%
  mutate(fill = factor(fill, levels = c(
    "laridae",
    "sulidae",
    "phalacrocoracidae",
    "spheniscidae",
    "alcidae",
    "stercorariidae",
    "procellariidae",
    "fregatidae"
  ))) 

plt_1c <- ggplot(seabird_meta) +
  geom_area(
    aes(
      x = ymd(paste0(collection_datemonth, "-01")),
      y = n,
      fill = fill
    ),
    position = position_fill(reverse = TRUE)
  ) +
  #scale_fill_brewer("Host",
    #na.translate = F,
    #palette = "Dark2",
    #direction = 1,
    #labels = str_to_title
  #) +
  scale_fill_manual('Host',
                    values =  c("laridae" = '#11A579',
                                "sulidae" = '#7F3C8D',
                                "phalacrocoracidae" = '#3969AC',
                                "spheniscidae" = '#F2B701',
                                "alcidae" = '#E73F74',
                                "stercorariidae" ='#80BA5A',
                                "procellariidae" = '#E68310',
                                "fregatidae" = '#008695'),
                    labels = str_to_title,
                    na.translate = F) + 
  theme_minimal() +
  scale_x_date("Year",
    expand = c(0, 0),
    breaks = pretty_breaks(n = 10)
  ) +
  scale_y_continuous("Proportion of NFLGs",
    expand = c(0, 0),
    limits = c(0, NA)
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(vjust = 0.5),
    legend.position = "inside",
    legend.position.inside = c(0, 1),
    legend.justification = c(0, 1)
  )


# 1d. phylogeography of H5N1/2022/R10
#treefiles <- "./phylo/ha_11414114_subsampled_traits_1000.trees"
#mcc_treefiles <- "./phylo/ha_11414114_subsampled_traits_mcc.tree"

#phylogeo_list <- FormatPhyloGeo(mcc_treefiles, treefiles)

#polygons_sf <- phylogeo_list[["polygons"]]
#edges_sf <- phylogeo_list[["edges"]]
#nodes_sf <- phylogeo_list[["nodes"]]

# load base map
#map <- ne_countries(scale = "medium", returnclass = "sf")

# Plot in GGplot
#plt_1d <- ggplot() +
 # geom_sf(data = map) +

  # Plot HPD polygons
  #geom_sf(
   # data = polygons_sf,
   # aes(fill = year),
   # lwd = 0,
   # alpha = 0.07
  #) +

  # Plot Branches
  #geom_sf(
 #   data = edges_sf,
   # lwd = 0.2
  #) +

  # Plot nodes
 # geom_sf(
   # data = nodes_sf,
   # size = 1.5,
   # aes(fill = year, colour = year, shape = is.na(label))
  #) +
  #scale_shape_manual(values = c(19, 1), ) +

  # Set graphical scales - must be fixed
  #scale_fill_viridis_c("Year",
    #limits = c(2021.5, 2023.8),
    # breaks=c(2019, 2020,2021,2022,2023,2024),
   # direction = -1,
   # option = "C"
  #) +
  #scale_colour_viridis_c("Year",
    #limits = c(2021.5, 2023.8),
    ## breaks=c(2019, 2020,2021,2022,2023,2024),
   # direction = -1,
  #  option = "C"
  #) +
  #coord_sf(
  #  ylim = c(10, 75),
  #  xlim = c(-20, 50),
  #  expand = TRUE
  #) +
  #scale_x_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0)) +
  #guides(
    #colour = guide_colourbar(
     # theme = theme(
       # legend.key.height = unit(0.75, "lines"),
       # legend.key.width = unit(10, "lines")
     # ),
      # title.position = 'left',
     # title.vjust = 1,
     # position = "bottom"
    #),
    #fill = guide_colourbar(
     # theme = theme(
      #  legend.key.height = unit(0.75, "lines"),
       # legend.key.width = unit(10, "lines")
     # ),
     # title.vjust = 1,
      # title.position = 'left',
      #position = "bottom"
    #),
    #shape = "none"
  #) +
  #theme_void() +
  #theme(
    #plot.margin = grid::unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
    # legend.text = element_text(size = 8),
   # legend.position = "bottom",
    #panel.spacing = unit(2, "lines"),
    #strip.background = element_blank()
  #)

# 1e. Photo
#bird_photo <- readJPEG("./misc/gull_steals_food.jpg")



#plt_1e <- ggdraw() +
  #draw_image(bird_photo) +
  #theme_void() +
  #theme(plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "pt"))

# 1f. Global tree with seabird species (+ seabird only skygrid)
new_tree <- read.beast('./phylo/ha_seabird_srd06_hmc_constant_hipstr.tree')
new_names <- as_tibble(new_tree) %>% dplyr::select(label) %>%
  drop_na() %>%
  mutate(isolate_id = str_extract(label, "EPI_ISL_(china_){0,1}\\d+[^.|]*")) %>%
  rename(new = label) %>%
  left_join(meta %>% 
  dplyr::select(isolate_id, old = tipnames))

plt_1f <- new_tree %>%
  # drop.tip(drop) %>%
  rename_taxa(new_names, old, new) %>%
  left_join(
    meta %>%
      left_join(new_seabirds %>%
                  dplyr::select(host_family = family) %>%
                  mutate(host_family = str_to_lower(host_family)) %>%
                  mutate(is_seabird = 1) %>%
                  distinct()) %>%
      mutate(fill = factor(host_family, levels = c(
        "laridae",
        "sulidae",
        "phalacrocoracidae",
        "spheniscidae",
        "alcidae",
        "stercorariidae",
        "procellariidae",
        "fregatidae"
      ))) %>%
      mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                               grepl('africa', collection_regionname) ~ 'africa',
                                               grepl('asia', collection_regionname) ~ 'asia',
                                               collection_regionname %in% c("arctic", "western canada" , "south" , "eastern canada",
                                                                            "midwest", "northeast", "west", "pacific northwest",
                                                                            "northern america","central america" ) ~ 'central & northern america',
                                               grepl("south georgia and the south", collection_countryname) ~'antarctica',
                                               grepl('south america', collection_regionname) ~ 'south america',
                                               grepl('australia|melanesia', collection_regionname) ~ 'australasia',
                                               .default = collection_regionname)) %>%
      rename(label = tipnames),
    by = "label"
  ) %>%
  #left_join(seabirds) %>%
  ggtree(mrsd = "2024-03-28") +
  
  geom_tiplab(aes(colour = fill),
              align = TRUE,
              size = 0
  ) +
  
  scale_colour_manual('Host',
                      values =  c("laridae" = '#11A579',
                                  "sulidae" = '#7F3C8D',
                                  "phalacrocoracidae" = '#3969AC',
                                  "spheniscidae" = '#F2B701',
                                  "alcidae" = '#E73F74',
                                  "stercorariidae" ='#80BA5A',
                                  "procellariidae" = '#E68310',
                                  "fregatidae" = '#008695'),
                      guide = 'none',
                      labels = str_to_title,
                      na.translate = F
  )+
  
  # tip colour + shape = new sequences
  new_scale_colour() + 
  geom_tippoint(aes(fill = fill, colour = !is.na(is_seabird)),
                shape = 21,
                size = 2
  ) +
  scale_fill_manual('Host',
                    values =  c("laridae" = '#11A579',
                                "sulidae" = '#7F3C8D',
                                "phalacrocoracidae" = '#3969AC',
                                "spheniscidae" = '#F2B701',
                                "alcidae" = '#E73F74',
                                "stercorariidae" ='#80BA5A',
                                "procellariidae" = '#E68310',
                                "fregatidae" = '#008695'),
                    guide = guide_legend(
                      keywidth = 1.5,
                      keyheight = 1,
                      ncol = 1,
                      order = 1
                    ),
                    labels = str_to_title,
                    na.translate = F
  ) +
  
  
  scale_colour_manual(
    labels = c(
      "TRUE" = "Seabird",
      "FALSE" = "Other"
    ),
    values = c(
      "TRUE" = "black",
      "FALSE" = "#ff000000"
    ),
    guide = "none"
  ) +
  
  
  
  # node colour to show pp support
  new_scale_fill() +
  # geom_nodepoint(aes(colour = posterior), alpha = 0.7) +
  # scale_color_distiller(palette = 'YlOrRd', direction = 1, 'Posterior Support',
  # guide = guide_colourbar(order = 4)) +
  
  geom_fruit(
    geom = geom_tile,
    mapping = aes(fill = collection_regionname),
    width = 0.3,
    # colour = "white",
    # pwidth = 1.2,
    offset = 0.05
  ) +
  
  scale_fill_brewer('Region',
                    palette = 'Blues',
                    #values = region_colours,
                     labels = str_to_title,
                     na.translate = F,
  ) + 
  #scale_fill_discrete("Host",
  # na.value = NA,
  
  #guide = guide_legend(keywidth = 1.5, keyheight = 1, ncol = 1, order = 1)
  #) +
  
  
  theme_tree2(
    plot.margin = unit(c(1, 1, 1, 1), units = "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(vjust = 0.5),
    legend.position = "inside",
    legend.position.inside = c(0, 1),
    legend.justification = c(0, 1)
  ) +
  
  scale_x_continuous(
    # limits = c(2000, 2023),
    "Year",
    breaks = seq(2016, 2024, 1)
  )


################################### OUTPUT #####################################
# Save output files, plots, or results
rh <- align_plots(plt_1c, plt_1f, axis = "r", align = "v")
lh <- align_plots(plt_1a, plt_1f, axis = "l", align = "v")

top <- plot_grid(lh[[1]], plt_1b, rh[[1]], align = 'h', axis = 'tb', nrow = 1, scale = 0.95, labels = 'AUTO')

plt1 <-plot_grid(top, plt_1f, nrow = 2, rel_heights = c(0.33, 0.66), labels = c('', 'D'))


# Save output files, plots, or results
rh <- align_plots(plt_1a, plt_1b, plt_1c, axis = "r", align = "v")
bot <- align_plots(plt_1f, plt_1c, axis = "b", align = "h")

top <- plot_grid(lh[[1]], plt_1b, rh[[1]], align = 'h', axis = 'tb', nrow = 1, scale = 0.95, labels = 'AUTO')

plt1 <-plot_grid( bot[[1]], plot_grid(plt_1a, plt_1b, bot[[2]], axis = "r", align = "v", ncol = 1), ncol = 2, rel_widths  = c(0.66, 0.33))


ggsave( "~/Downloads/seabird_fig1.jpeg",
        plt1,
       height = 10,
       width = 10,
       dpi = 360
)

top <- plot_grid(lh[[1]], plt_1b, rh[[1]], align = 'h', axis = 'tb', nrow = 1, scale = 0.95, labels = 'AUTO')

plot_grid(top, plt_1f , nrow = 2, rel_heights = c(0.33, 0.66), labels = c('', 'D'))

lh <- plot_grid(top[[1]], plt_1b, bttm[[1]], axis = "lr", align = "v", nrow = 3, labels = "AUTO")
mid <- plot_grid(top[[2]], bttm[[3]],
  axis = "lr", align = "v", nrow = 2, labels = c("D", "E"),
  rel_heights = c(1, 1)
)



cowplot::plot_grid(lh,
  mid,
  bttm[[2]],
  ncol = 3,
  rel_widths = c(0.9, 0.9, 1),
  labels = c("", "", "F")
)


#################################### END #######################################
################################################################################
