library(sf)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(igraph)
library(ggplot2)
library(rnaturalearth)

# ----  SETUP  ----
#paths
INPUT_DIR <- "EGM_2019_SHP_20190312/DATA/Countries"
OUTPUT_DIR <- "results"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)
}
if (!dir.exists(INPUT_DIR)) stop(paste("Input directory not found:", INPUT_DIR))
country_dirs <- list.dirs(INPUT_DIR, full.names = FALSE, recursive = FALSE)
country_dirs <- country_dirs[country_dirs != ""]

#CRS (Coordinate Reference Systems)
CRS_METRIC <- 3035  # ETRS89-LAEA for distances
CRS_OUTPUT <- 4326  # WGS84 (Lat/Lon for output)

#data cleaning rules, from the EuroGlobalMap spec
MISSING_STRINGS <- c("N_A", "UNK", "N_P", "Not Applicable", "Unknown", "", "N/A", "NA")
VALID_TFC <- c(15, 31, 34) #"valid" station codes: 15: station, 31: joint station, 32: halt, 34: terminal
VALID_EXS <- c(28) # 'operational' lines
VALID_TUC <- c(26, 45) # 26: 'passenger', 45: 'general' lines

#topology rules, from the spec
MIN_STATION_DIST <- 200 #mrge stations closer than this
SNAP_TOLERANCE <- 50  #tolerance to match station to line endpoint

# ---- EXECUTION ----
#Loop through every folder found
for (target_country in country_dirs){
  
  message(paste0("\n--- Processing: ", target_country, " ---"))

  CURRENT_OUT_DIR <- file.path(OUTPUT_DIR, target_country) #output path
  if (!dir.exists(CURRENT_OUT_DIR)) dir.create(CURRENT_OUT_DIR)

  NODES_PATH <- file.path(INPUT_DIR, target_country, "RailrdC.shp") #input path
  LINES_PATH <- file.path(INPUT_DIR, target_country, "RailrdL.shp")  
  
  # ---- 1) CLEANING NODES ----
  message("Loading and cleaning station data")
  if (!file.exists(NODES_PATH)){
    warning(paste("Nodes file not found for", target_country))
    next
  }
  
  stations_raw <- st_read(NODES_PATH, quiet=TRUE)
  
  raw_tfc_counts <- table(stations_raw$TFC) #a specific check, due to HU bad data formatting (no proper stations)
  
  stations_prep <- stations_raw %>%
    mutate(
      TFC = ifelse(TFC<0, NA, TFC), # check for negative "missing" codes
      Label = coalesce(NAMN1, NAMN2),
      Label = ifelse(Label %in% MISSING_STRINGS | is.na(Label), "Unknown", Label),
      Country_ISO3 = ifelse(ICC %in% MISSING_STRINGS | ICC=="XX", "UNK", ICC)
    ) %>%
    filter(TFC %in% VALID_TFC) %>% 
    st_transform(CRS_METRIC)
  
  #still for the HU case: show that everything had to be filtered out
  if (nrow(stations_prep)==0 && nrow(stations_raw)>0){
    message(paste("  [!] WARNING: 0 stations kept. Available TFC codes were:", paste(names(raw_tfc_counts), collapse=", ")))
  }
  
  #merge stations that are too close
  if (nrow(stations_prep)>0){ #some countries may have none
    station_coords <- st_coordinates(stations_prep)
    dist_mat <- st_distance(stations_prep)
    hc <- hclust(as.dist(dist_mat), method="single")
    stations_prep$cluster_id <- cutree(hc, h=MIN_STATION_DIST)
    stations_clean <- stations_prep %>% group_by(cluster_id) %>%
                      summarise(
                        geometry = st_centroid(st_union(geometry)),
                        Label = first(Label),
                        Country_ISO3 = first(Country_ISO3),
                        Original_ID = first(RStationID)
                      ) %>% ungroup() %>% mutate(Type="Station")
  } else{
    stations_clean <- stations_prep # empty case
  }
  
  message(sprintf("  -> Reduced %d raw points to %d unique stations.", nrow(stations_raw), nrow(stations_clean)))
  
  # ---- 2) CLEANING EDGES ---
  message("Loading and cleaning segments")
  if (!file.exists(LINES_PATH)){
    warning(paste("Lines file not found for", target_country))
    next
  }
  
  lines_raw <- st_read(LINES_PATH, quiet=TRUE)
  lines_clean <- lines_raw %>% filter(EXS %in% VALID_EXS, TUC %in% VALID_TUC) %>%
                 st_transform(CRS_METRIC) %>% distinct(geometry, .keep_all=TRUE) #remove duplicates
  
  #other two cases with "improper formatting" ( different from the previous one) where found
  if (nrow(lines_clean)==0 && nrow(lines_raw)>0){
    message(paste("  [!] WARNING: 0 lines kept. Check EXS/TUC codes."))
    message(paste("      Available TUC:", paste(unique(lines_raw$TUC), collapse=", ")))
    message(paste("      Available EXS:", paste(unique(lines_raw$EXS), collapse=", ")))
  }
  
  message(sprintf("  -> Kept %d operational line segments.", nrow(lines_clean)))
  
  # ---- 3) TOPOLOGY BUILDING  ----
  #Let's now build the topology and create the 'edges.csv' file.
  #Since connectivity in 'RailrdL' is purely geometrical, there is no "Source ID" or "Target ID" column in the lines file. 
  #Instead, if a line ends at coordinate (X,Y) and another line starts at (X,Y), they are connected; we exploit this concept to identify the nodes directly from the tracks
  #(we're assuming, from the spec., that valid connections only happen at the ends of these segments, not in the middle).
  #We in fact extract all the endpoints first; if multiple lines meet at the same coordinate we treat them as a single unique node. This topological 'node', however, 
  #may be a station or a junction: we check its proximity to our cleaned station list. If a station (of 'nodes.csv') is within 50 meters, that node inherits the s
  #tation attributes, otherwise it is labeled as a junction. We can now determine the "source" and "target" for each segment stored in 'RailrdL'.
  message("Building topology")
  
  #safety check for empty lines
  if (nrow(lines_clean)==0){
    message(paste("No valid lines for", target_country, "- Skipping topology."))
    next
  }
  
  #extracting the coordinates of lines endpoints
  line_coords <- st_coordinates(lines_clean)
  starts <- line_coords[match(unique(line_coords[,"L1"]), line_coords[,"L1"]), 1:2]
  ends <- line_coords[nrow(line_coords) - match(unique(line_coords[,"L1"]), rev(line_coords[,"L1"]))+1, 1:2]
  all_endpoints <- rbind(starts, ends) #combining all endpoints
  all_endpoints_df <- as.data.frame(all_endpoints) %>% setNames(c("X", "Y"))
  
  #finding unique nodes (unique within 1m tolerance)
  nodes_sf <- st_as_sf(all_endpoints_df, coords=c("X", "Y"), crs=CRS_METRIC) %>% 
              group_by(geometry) %>% slice(1) %>%ungroup() %>%
              mutate(nodeID = row_number()) #for the naming requirement 
  
  message(sprintf("  -> Identified %d topological nodes", nrow(nodes_sf)))
  
  #matching the cleaned stations to the topological nodes
  if (nrow(stations_clean)>0){ #"if not HU"
    nearest_idx <- st_nearest_feature(stations_clean, nodes_sf)
    dists <- st_distance(stations_clean, nodes_sf[nearest_idx,], by_element=TRUE)
    
    station_mapping <- stations_clean %>%
                      mutate(
                        mapped_node_id = nodes_sf$nodeID[nearest_idx],
                        dist_m = as.numeric(dists) 
                      ) %>% filter(dist_m <= SNAP_TOLERANCE) %>% 
                      st_drop_geometry() %>%
                      select(mapped_node_id, Label, Country_ISO3, Type) %>%
                      distinct(mapped_node_id, .keep_all=TRUE)
    
    final_nodes <- nodes_sf %>%
                  left_join(station_mapping, by = c("nodeID" = "mapped_node_id")) %>%
                  mutate(
                    Label = replace_na(Label, "Junction"),
                    Type = replace_na(Type, "Junction"),
                    Country_ISO3 = replace_na(Country_ISO3, "UNK"))
  } else{ #"if HU"
    final_nodes <- nodes_sf %>% mutate(Label = "Junction", Type = "Junction", Country_ISO3 = "UNK")
  }
  
  #building the edge list:
  starts_sf <- st_as_sf(as.data.frame(starts), coords=1:2, crs=CRS_METRIC)
  ends_sf <- st_as_sf(as.data.frame(ends), coords=1:2, crs=CRS_METRIC)
  start_matches <- st_equals(starts_sf, final_nodes) #find which nodeID corresponds to start and end
  end_matches <- st_equals(ends_sf, final_nodes)
  
  get_id <- function(x) ifelse(length(x)> 0, final_nodes$nodeID[x[[1]]], NA) #just a helper function
  
  final_edges <- data.frame(nodeID_from=sapply(start_matches, get_id), nodeID_to=sapply(end_matches, get_id)) %>% 
                filter(!is.na(nodeID_from) & !is.na(nodeID_to)) %>% filter(nodeID_from != nodeID_to) # eemove self-loops

#---- 4) LOGICAL DEGREE ----
  #We can now move to the plotting phase. We will produce a plot that shows both stations and junctions, highligting only the former.
  #The dot size representing the stations will be proportional to its degree. The result, however, depends a lot on which we consider to be the nodes: 
  #stations alone or stations and junctions, but it also depends on the type of stations we consider.
  #Here we will consider stations only, and we will focus on people "flux", thus on passengers stations, therefore the degree will be due to 
  #station-station connections. To determine such connections we implement a "basic-search" algorithm:
  #start a "walker" at a station S. Move to all immediate neighbors. From those, keep moving to their neighbors. Stop if the walker: hits another 
  #station (and record this connection), hits a dead end, hits a node it has already visited in the specific search.
  #In this way we're also able to get the adjacency matrix of the network for station-station connections.
  message("Calculating station-to-station degree")
  
  nodes_df <- data.frame(name=as.character(final_nodes$nodeID), stringsAsFactors=FALSE) #ensuring a proper dataframe
  
  g_phys <- graph_from_data_frame(final_edges, directed=FALSE, vertices=nodes_df) #building the graph
  V(g_phys)$type <- final_nodes$Type
    
  station_ids <- final_nodes$nodeID[final_nodes$Type=="Station"]
  logical_degrees <- setNames(rep(0, nrow(final_nodes)), final_nodes$nodeID) #degree counter
  logical_edges_list <- list() #we store the logical egdes
  
  if (length(station_ids)>0){ #only run the crawler if there are actually stations
    adj_list <- as_adj_list(g_phys)
    is_station_vec <- V(g_phys)$type=="Station"
    names(is_station_vec) <- V(g_phys)$name
    
    for (root in station_ids){      
      root_char <- as.character(root)
      queue <- list(root_char)
      visited <- c(root_char)
      
      found_targets <- c() 
      
      head <- 1
      while(head<=length(queue)){
        curr <- queue[[head]]
        head <- head + 1
        
        neighbors <- adj_list[[curr]]
        
        for(nbr in neighbors){
          if(nbr %in% visited) next
          visited <- c(visited, nbr) # mark as visited to prevent loops
          if(is_station_vec[nbr]){ #found a station --> stop
            found_targets <- c(found_targets, nbr)
          } else { #it's a junction --> continue
            queue[[length(queue) + 1]] <- nbr
          }
        }
      }
      
      unique_targets <- unique(found_targets) #degree counting
      logical_degrees[root_char] <- length(unique_targets)
      
      if (length(unique_targets)>0){ #storing logical edges
        logical_edges_list[[length(logical_edges_list)+1]] <- data.frame(
          nodeID_from=root, nodeID_to=as.integer(unique_targets))
      }
    }
  }
  logical_edges_df <- bind_rows(logical_edges_list)#bind all the small dataframes
  final_nodes$Logical_Degree <- logical_degrees[as.character(final_nodes$nodeID)] #attach degree to nodes 

  
  # ---- 5) EXPORTING ----
  #country code -> country name
  country_lookup <- ne_countries(scale="medium", returnclass="sf") %>% st_drop_geometry() %>% select(iso_a2, admin) %>% distinct()
  
  final_nodes_wgs <- st_transform(final_nodes, CRS_OUTPUT)
  coords_wgs <- st_coordinates(final_nodes_wgs)
  
  nodes_export <- final_nodes_wgs %>% st_drop_geometry() %>%
                left_join(country_lookup, by=c("Country_ISO3"="iso_a2")) %>%
                mutate(
                  latitude = coords_wgs[,2],
                  longitude = coords_wgs[,1],
                  country_name = ifelse(is.na(admin), "Unknown", admin), #if lookup failed
                  nodeLabel = Label
                ) %>% select(nodeID, nodeLabel, latitude, longitude, country_name, 
                             country_ISO3=Country_ISO3, Logical_Degree, Type)
  
  #storing results
  file_nodes <- file.path(CURRENT_OUT_DIR, paste0("network_nodes_", target_country, ".csv"))
  file_edges <- file.path(CURRENT_OUT_DIR, paste0("network_edges_", target_country, ".csv"))
  file_logical <- file.path(CURRENT_OUT_DIR, paste0("network_logical_edges_", target_country, ".csv"))
    
  write_csv(nodes_export, file_nodes)
  write_csv(final_edges, file_edges)
  write_csv(logical_edges_df, file_logical)
  message(paste0("Results stored in folder ", CURRENT_OUT_DIR))
  
  # ---- 6) PLOTTING ----
  europe <- ne_countries(continent="Europe", returnclass="sf", scale="medium") #europe map
  #joining coordinates before plotting
  nodes_tbl <- nodes_export %>% select(nodeID, lat=latitude, lon=longitude)
  edges_plot <- final_edges %>%
                left_join(nodes_tbl, by=c("nodeID_from"="nodeID")) %>%
                rename(lat1=lat, lon1=lon) %>%
                left_join(nodes_tbl, by=c("nodeID_to"="nodeID")) %>%
                rename(lat2=lat, lon2=lon)
  
  p <- ggplot() +
    geom_sf(data=europe, fill="#f5f5f5", color="white") +
    #edges
    geom_segment(data=edges_plot, 
                 aes(x=lon1, y=lat1, xend=lon2, yend=lat2),
                 col="grey60", linewidth=0.1, alpha=0.5) +
    #junctions
    geom_point(data=nodes_export %>% filter(Type=="Junction"),
               aes(x=longitude, y=latitude),
               col="grey40", size=0.2, alpha=0.3) +
    #stations
    geom_point(data = nodes_export %>% filter(Type == "Station"), 
               aes(x = longitude, y = latitude, size = Logical_Degree),
               col = "firebrick", alpha = 0.9) +
    #styling
    scale_size_continuous(range=c(0.5, 4), name="Logical degree") +
    #zoom to specific country
    coord_sf(xlim = c(min(nodes_export$longitude)-1, max(nodes_export$longitude)+1), 
             ylim = c(min(nodes_export$latitude)-1, max(nodes_export$latitude)+1)) + 
    theme_void() + theme(legend.position="bottom") +
    labs(title=paste0("Rail network: ", target_country), 
         subtitle = "Red: Stations (sized by connectivity) | Grey: junctions & physical tracks")
  
  ggsave(file.path(CURRENT_OUT_DIR, paste0("network_map_", target_country, ".png")), p, width = 12, height = 10, bg = "white")
  message(paste0("  -> Plot saved as 'network_map_", target_country, ".png'."))
}


library(scales)

RESULTS_DIR <- "results"
SUMMARY_FILE <- file.path(RESULTS_DIR, "all_countries_stats_summary.csv")
if (!dir.exists(RESULTS_DIR)) stop("Results directory not found. Run data_proj_final.R first.")
country_dirs <- list.dirs(RESULTS_DIR, full.names = FALSE, recursive = FALSE)
country_dirs <- country_dirs[country_dirs != ""]

global_stats_log <- data.frame(
  Country = character(),
  Nodes = integer(),
  Edges = integer(),
  Giant_Comp_Size = integer(),
  Giant_Comp_Percent = numeric(),
  Is_Connected = logical(),
  Avg_Degree = numeric(),
  Global_Clustering = numeric(),
  Assortativity = numeric(),
  stringsAsFactors = FALSE
)

for (target_country in country_dirs){
  
  message(paste0("\n--- Analyzing: ", target_country, " ---"))
  CURRENT_DIR <- file.path(RESULTS_DIR, target_country)
  
  # we use the LOGICAL files (station-to-station)
  file_nodes <- file.path(CURRENT_DIR, paste0("network_nodes_", target_country, ".csv"))
  file_edges <- file.path(CURRENT_DIR, paste0("network_logical_edges_", target_country, ".csv"))
  if (!file.exists(file_nodes) || !file.exists(file_edges)){
    warning(paste("Missing logical network files for", target_country, "- Skipping."))
    next
  }
  
  nodes_df <- read_csv(file_nodes, show_col_types=FALSE)
  edges_df <- read_csv(file_edges, show_col_types=FALSE)

  stations_only <- nodes_df %>% filter(Type=="Station") # Only keep stations
  if (nrow(stations_only)==0){
    message("  -> No stations found (0 nodes). Skipping analysis.")
    next
  }
    
  g <- graph_from_data_frame(edges_df, directed=FALSE, vertices=stations_only) #undirected graph from logical edge list
  message(sprintf("  -> Nodes: %d | Edges: %d", vcount(g), ecount(g)))
  
  # 1)MACROSCOPIC ANALYSIS (robustness & components)
  #connected components
  comps <- components(g)
  giant_id <- which.max(comps$csize)
  giant_size <- comps$csize[giant_id]
  is_connected <- (comps$no == 1)
  giant_pct <- round((giant_size/vcount(g))*100, 2)
    
  message(sprintf("  -> Connected: %s | Giant Component: %d nodes (%s%%)", is_connected, giant_size, giant_pct))
  
  #2) MICROSCOPIC ANALYSIS (centrality & mixing)
  message("  -> Computing centralities")
  
  # - degree & assortativity
  deg <- degree(g)
  avg_k <- mean(deg)
  assortativity_val <- assortativity_degree(g, directed=FALSE)
  
  # - global clustering 
  clustering_val <- transitivity(g, type="global")
  
  # - centrality 
  bet_vals <- betweenness(g, normalized=TRUE)  #Betweenness
  close_vals <- harmonic_centrality(g, normalized=TRUE) #Harmonic centrality
  #Katz Centrality
  adj_mat <- as_adjacency_matrix(g, sparse=TRUE)
  eigen_res <- eigen(adj_mat, only.values=TRUE, symmetric=TRUE)
  lambda_max <- max(eigen_res$values)
  alpha_val <- 0.85 / lambda_max #std value
  katz_vals <- alpha_centrality(g, alpha = alpha_val)
  
  stations_only$Betweenness <- bet_vals[as.character(stations_only$nodeID)]
  stations_only$Closeness <- close_vals[as.character(stations_only$nodeID)]
  stations_only$Katz <- katz_vals[as.character(stations_only$nodeID)]
  
  #3) MESOSCALE ANALYSIS (communities)
  comm_det <- cluster_louvain(g) #Louvian method
  stations_only$Community <- membership(comm_det)[as.character(stations_only$nodeID)]
  modularity_score <- modularity(comm_det)
  message(sprintf("  -> Detected %d communities (Modularity Q: %.3f)", length(unique(stations_only$Community)), modularity_score))
  
  #4) EXPORTING
  row_stats <- data.frame(
    Country = target_country,
    Nodes = vcount(g),
    Edges = ecount(g),
    Giant_Comp_Size = giant_size,
    Giant_Comp_Percent = giant_pct,
    Is_Connected = is_connected,
    Avg_Degree = round(avg_k, 3),
    Global_Clustering = round(clustering_val, 4),
    Assortativity = round(assortativity_val, 4)
  )
  global_stats_log <- rbind(global_stats_log, row_stats)
  
  metrics_export <- stations_only %>% select(nodeID, Betweenness, Closeness, Katz, Community) #removing redundant info
  write_csv(metrics_export, file.path(CURRENT_DIR, paste0("station_metrics_", target_country, ".csv")))
  
  #5) PLOTTING
  #degree distribution
  p_deg <- ggplot(data.frame(Degree=deg), aes(x=Degree)) +
    geom_bar(fill="steelblue", color="white") +
    labs(title=paste0("Degree distribution: ", target_country), subtitle=paste0("Assortativity: ", round(assortativity_val, 3)),
         x="Degree (k)", y="Frequency") + theme_minimal()
  ggsave(file.path(CURRENT_DIR, paste0("plot_degree_", target_country, ".png")), p_deg, width=8, height=6)
  
  #communities on a map
  europe_map <- ne_countries(continent="Europe", returnclass="sf", scale="medium")
  p_comm <- ggplot() + geom_sf(data = europe_map, fill="#f0f0f0", color="white") +
    geom_point(data=stations_only, 
               aes(x=longitude, y=latitude, color=as.factor(Community)), size=1.5, alpha=0.8) +
    coord_sf(xlim=c(min(stations_only$longitude)-1, max(stations_only$longitude)+1), 
             ylim=c(min(stations_only$latitude)-1, max(stations_only$latitude)+1)) +
    labs(title=paste0("Rail communities: ", target_country), 
         subtitle=paste0("Modularity Q: ", round(modularity_score, 3)), color="Community") + 
    theme_void() + theme(legend.position="none")
  ggsave(file.path(CURRENT_DIR, paste0("plot_communities_", target_country, ".png")), p_comm, width = 10, height = 8)
  
  #centrality maps (Betweenness, Closeness, Katz)
  nodes_coord <- stations_only %>% select(nodeID, lon = longitude, lat = latitude)
  edges_plot <- edges_df %>% 
    left_join(nodes_coord, by=c("nodeID_from" = "nodeID")) %>% rename(lon1=lon, lat1=lat) %>%
    left_join(nodes_coord, by = c("nodeID_to" = "nodeID")) %>% rename(lon2=lon, lat2=lat)
  
  centrality_measures <- c("Betweenness", "Closeness", "Katz")
  for (measure in centrality_measures){
      limit_val <- quantile(stations_only[[measure]], 0.99, na.rm=TRUE) #kets guard from outliers
      if (limit_val==0) limit_val <- max(stations_only[[measure]], na.rm=TRUE)
      
      p_cent <- ggplot() + geom_sf(data=europe_map, fill="#f0f0f0", color="white") +
        geom_point(data=stations_only, aes(x=longitude, y=latitude, color=.data[[measure]]), 
                   size=0.6, alpha=0.9) + #stations
        scale_color_viridis_c(option="viridis", direction=1,
                              limits=c(0, limit_val), oob=scales::squish) + #clamps the outliers to the max limit
        coord_sf(xlim=c(min(stations_only$longitude)-1, max(stations_only$longitude)+1), 
                 ylim=c(min(stations_only$latitude)-1, max(stations_only$latitude)+1)) +
        labs(title=paste0(measure, " Centrality: ", target_country),
             subtitle="Color scale capped at 99th percentile to handle outliers",
             color = measure) +
        theme_void() + theme(legend.position = "right")
      ggsave(file.path(CURRENT_DIR, paste0("plot_", tolower(measure), "_", target_country, ".png")), 
             p_cent, width=10, height=8)
  }
}

write_csv(global_stats_log, SUMMARY_FILE)
message(paste("Global summary saved to:", SUMMARY_FILE))