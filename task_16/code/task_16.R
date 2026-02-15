# References:
# [1] A. Arenas, A. Diaz-Guilera, R. Guimerà, Communication in Networks with Hierarchical Branching, Phys. Rev. Lett., 86, 3196-3199, 2001
# [2] E. Pablo,  J. Gomez-Gardenes, Y. Moreno, Improved routing strategies for Internet traffic delivery, Phys. Rev. E, 70, 056105, 2004
# [3] E. Pablo,  J. Gomez-Gardenes, Y. Moreno, Dynamics of jamming transitions in complex networks, Europhys. Lett., 71, 325–331, 2005


# -----------------------------------------------------------------------------
# STEP 1: TOPOLOGY GENERATION
# -----------------------------------------------------------------------------
library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

generate_tree_topology <- function(branching_factor=3, levels=4){ # [1] uses z=3, m=4
  n_nodes <- (branching_factor^levels-1)/(branching_factor-1) # number of nodes
  g <- make_tree(n= n_nodes, children=branching_factor, mode= "undirected")
    
  V(g)$type <- "Tree"
  V(g)$id <- 1:vcount(g)
  return(g)
}

generate_scalefree_topology <- function(n_nodes=500, m_links=2){ #Barabasi-Albert model
  g <- sample_pa(n = n_nodes, power = 1, m = m_links, directed = FALSE)
    
  V(g)$type <- "Scale-Free"
  V(g)$id <- 1:vcount(g)
  return(g)
}

generate_random_topology <- function(n_nodes = 500, p_link = 0.02){ #Gilbert model
  g <- sample_gnp(n = n_nodes, p = p_link, directed = FALSE)
  
  comps <- components(g) #ensuring it is connected
  if (comps$no > 1){ #extract largest component if disconnected
    g <- induced_subgraph(g, which(comps$membership == which.max(comps$csize)))
  }

  V(g)$type <- "Random"
  V(g)$id <- 1:vcount(g)
  return(g)
}

get_distance_matrix <- function(g){ #pre-compute shortest paths
  distances(g, algorithm = "unweighted") #matrix D where D_i,j is the topological distance
}

# -----------------------------------------------------------------------------
# STEP 2: THE DYNAMICS ENGINE
# -----------------------------------------------------------------------------

#The routing function, from [2] and [3]
get_next_hop <- function(current_node, dest_node, neighbors_list, distance_matrix, queue_sizes, h){ #Decides which neighbor to send the packet to
  candidates <- neighbors_list[[current_node]] #get all neighbors of the current node
  if (dest_node %in% candidates){ #if destination is a neighbor, send directly
    return(dest_node)
  }
  
  #Calculate "Effective Distance" for all neighbors:  d_eff_i = h*d_i +(1−h)c_i
  #d_i=dist_matrix[i, dest_node], the distance from neighbor i to destination
  #c_i=queue_sizes[i], the congestion at neighbor i
  #h=1.0: pure shortest path, h<1.0: traffic aware 
  d_i <- distance_matrix[candidates, dest_node]
  c_i <- queue_sizes[candidates]
  costs <- (h*d_i)+((1 - h)*c_i)

  min_cost <- min(costs) #pick neighbor with minimum cost
  best_candidates <- candidates[which(costs == min_cost)] 
  if (length(best_candidates) > 1){ #break ties randomly
    return(sample(best_candidates, 1))
  } else {
    return(best_candidates)
  }
}

#Simulation Loop
run_simulation <- function(g, p_rate, h_param, n_steps=1000){
  # --- Initialization ---
  num_nodes <- vcount(g)
  dist_matrix <- get_distance_matrix(g)
  adj_list <- as_adj_list(g)
  
  #queues: a list of vectors --> queue[[i]] holds packet IDs at node i
  node_queues <- vector("list", num_nodes)
  for(i in 1:num_nodes) node_queues[[i]] <- integer(0)
  
  #'packet_registry': to track arrival times (to measure <T>, qunatificator used in [2]; here we'll actually use <T>=T_avg instead of <T>=<T_max>)
  #Instead of a df we use a pre-allocated matrix: lighter, allows for instant elements look-up, allows us to avoid increasing the size iteratively
  max_capacity <- ceiling(n_steps*p_rate*1.5) + 200 #etimation of the size needed
  packet_registry <- matrix(NA, nrow = max_capacity, ncol = 4) #Columns: 1=source, 2=dest, 3=t_create, 4=t_arrive. Row Index=packet ID.
  colnames(packet_registry) <- c("source", "dest", "t_create", "t_arrive")
  
  #'history': to track the order parameter
  history_N <- numeric(n_steps) # active packets
  packet_counter <- 0
  
  # --- Time evolution ---
  for (t in 1:n_steps){
    # --- Packet generation ---
    n_new <- rpois(1, lambda=p_rate) #p information packets are created, at a Poisson rate [1]
    if (n_new > 0){ #check if matrix expansion is needed
      if ((packet_counter + n_new) > nrow(packet_registry)) {
        cat("Matrix had been expanded\n")
        packet_registry <- rbind(packet_registry, matrix(NA, nrow=1000, ncol=4))
      }
      
      new_ids <- (packet_counter + 1):(packet_counter + n_new)
      
      #assigning, randomly, the source and destination of the newly created packets
      sources <- sample(1:num_nodes, n_new, replace=TRUE)
      dests <- sample(1:num_nodes, n_new, replace=TRUE)
      same_idx <- which(sources == dests) #to ensure source != dest
      while(length(same_idx) > 0){ 
        dests[same_idx] <- sample(1:num_nodes, length(same_idx), replace=TRUE)
        same_idx <- which(sources == dests)
      }

      #update registries  
      packet_registry[new_ids, 1] <- sources
      packet_registry[new_ids, 2] <- dests
      packet_registry[new_ids, 3] <- t
      
      packet_counter <- packet_counter + n_new

      #update queues: the new packets "pile up" in the queues of the corresponding sources
      for (k in 1:n_new){
        src <- sources[k]
        node_queues[[src]] <- c(node_queues[[src]], new_ids[k])
      }
    }
    
    # --- Actual routing  ---
    active_nodes <- which(sapply(node_queues, length) > 0) #nodes that have packets to send
    current_queue_sizes <- sapply(node_queues, length) #synchronous update approach --> we take a snapshot of the queue sizes before the "move" step
     
    transfers <- list() #list of (packet_id, from_node, to_node), see later
    for (u in active_nodes){ #routers with empty queues do nothing
      pkt_id <- node_queues[[u]][1] #the protocol is FIFO --> pick just 1 packet, the "oldest" one
      dest_node <- packet_registry[pkt_id, 2]
      next_node <- get_next_hop(current_node=u, dest_node=dest_node, neighbors_list=adj_list, distance_matrix=dist_matrix, 
                                queue_sizes=current_queue_sizes, h=h_param) #calculate the next hop using d_eff
      transfers[[length(transfers) + 1]] <- list(pkt=pkt_id, from=u, to=next_node) #storing router u's "plan of action"
    }
    
    # --- Executing the "plans of action" ---
    delivered_count <- 0
    for (tr in transfers){
      node_queues[[tr$from]] <- node_queues[[tr$from]][-1] #removing the packet from sources' queues 
      
      pkt_dest <- packet_registry[tr$pkt, 2] 
      if (tr$to == pkt_dest){ # if it is the destination: delivered --> remove from the network
        packet_registry[tr$pkt, 4] <- t #update t_arrive
        delivered_count <- delivered_count + 1
      } else{ #if not, append to destination queue
        node_queues[[tr$to]] <- c(node_queues[[tr$to]], tr$pkt)
      }
    }
    
    # --- Metric ---
    # N(t) = 'total active' = total_generated - total_delivered = sum of all queue lengths
    total_active <- sum(sapply(node_queues, length))
    history_N[t] <- total_active

    if (t %% 100 == 0) cat(".") #just to print progress every 100 steps (since computation is slow)
  }
  cat("\n")
  
  # --- Getting the results ---
  #1)microscopic analysis: betweenness
  node_betw <- betweenness(g, normalized=TRUE)
  final_queues <- sapply(node_queues, length)
  
  micro_stats <- data.frame(
    node_id = 1:num_nodes,
    queue_size = final_queues,
    betweenness = node_betw,
    topology_type = V(g)$type[1]
  )
  
  #2)an efficiency measure: average travel time, T_avg
  arrived_packets <- as.data.frame(packet_registry[1:packet_counter, , drop=FALSE]) #conversion to a df only of valid parts
  arrived_packets <- arrived_packets[!is.na(arrived_packets$t_arrive), ] #for packets that actually arrived
  
  if (nrow(arrived_packets) > 0){
    avg_time <- mean(arrived_packets$t_arrive - arrived_packets$t_create)
  } else {
    avg_time <- NA
  }
  
  #3)order parameter: slope of N(t)/p [1], [3]
  window <- floor(0.8 * n_steps):n_steps #we consider just the last part to ensure the system is not in a temporary regime
  #(equivalent to lim_t->+inf [3])
  if (length(window) > 1){
    fit <- lm(history_N[window] ~ window) #linear regression slope
    slope <- coef(fit)[["window"]]
    rho <- slope/p_rate
    rho <- max(0, min(1, rho))
  } else {
    rho <- 0
  }
  
  return(list(
    N_t = history_N,
    micro_stats = micro_stats,
    avg_travel_time = avg_time,
    rho = rho,
    p = p_rate,
    h = h_param
  ))
}

# -----------------------------------------------------------------------------
# STEP 3: RUNNING THE ANALYSIS
# -----------------------------------------------------------------------------
#Configuration
#simulation paramters
TIME_STEPS <- 1500
P_RANGE <- seq(1, 40, by=2) 
THRESHOLD <- 0.02 #threshold for detecting pc (rho > thresh)

#generating all the configurations we will test
N_std <- 200; m_std <- 2
N_high <- 500; m_high <- 4

#p for RN: p=2*m/(N-1)
p_std <- 2*m_std/(N_std-1)
p_incN <- 2*m_std/(N_high-1)
p_incM <- 2*m_high/(N_std-1)
p_incNM <- 2*m_high/(N_high-1)

h_values <- c(1.0, 0.8, 0.3)

n_runs <- 5

experiment_configs <- list()

#tree
tree_params <- list(
  list(z=3, m=4, label="Tree_z3m4"),
  list(z=3, m=5, label="Tree_z3m5"),
  list(z=3, m=3, label="Tree_z3m3"),
  list(z=4, m=4, label="Tree_z4m4")
)

for (tp in tree_params) {
  for (h in h_values) {
    experiment_configs[[length(experiment_configs) + 1]] <- list(
      name = paste0(tp$label, "_h", h),
      type = "Tree",
      params = list(z=tp$z, m=tp$m),
      h = h,
      n_runs = n_runs
    )
  }
}

#SF
sf_params <- list(
  list(n=N_std, m=m_std, label="SF_std"),
  list(n=N_high,m=m_std, label="SF_incN"),
  list(n=N_std, m=m_high, label="SF_incM"),
  list(n=N_high, m=m_high, label="SF_incNM")
)

for (sp in sf_params) {
  for (h in h_values) {
    experiment_configs[[length(experiment_configs)+1]] <- list(
      name = paste0(sp$label, "_h", h),
      type = "Scale-Free",
      params = list(n=sp$n, m=sp$m),
      h = h,
      n_runs = n_runs
    )
  }
}

#rnadom netw
rn_params <- list(
  list(n=N_std, p=p_std, label="Rnd_std"),
  list(n=N_high, p=p_incN, label="Rnd_incN"),
  list(n=N_std, p=p_incM, label="Rnd_incM"),
  list(n=N_high, p=p_incNM, label="Rnd_incNM")
)

for (rp in rn_params){
  for (h in h_values){
    experiment_configs[[length(experiment_configs)+1]] <- list(
      name = paste0(rp$label, "_h", h),
      type = "Random",
      params = list(n=rp$n, p_link=rp$p),
      h = h,
      n_runs = n_runs
    )
  }
}

#topology generator
get_topology_from_config <- function(type, params){
  if (type == "Tree") {
    return(generate_tree_topology(branching_factor=params$z, levels=params$m))
  } else if (type == "Scale-Free") {
    return(generate_scalefree_topology(n_nodes=params$n, m_links=params$m))
  } else if (type == "Random") {
    return(generate_random_topology(n_nodes=params$n))
  } else {
    stop(paste("Unknown topology type:", type))
  }
}

#Run
run_averaged_sweep <- function(config, p_seq){
  cat(paste0("\n--- Processing: ", config$name, " ---\n"))
  results_list <- list()
  pb <- txtProgressBar(min = 0, max = length(p_seq), style = 3)
  
  for (i in seq_along(p_seq)){
    p_val <- p_seq[i]
    rho_samples <- numeric(config$n_runs)
    time_samples <- numeric(config$n_runs)
    
    #multiple runs --> averaging the results
    for (r in 1:config$n_runs){
      g_instance <- get_topology_from_config(config$type, config$params) 
      sim_out <- run_simulation(g_instance, p_rate = p_val, h_param = config$h, n_steps = TIME_STEPS)
      
      rho_samples[r] <- sim_out$rho
      t_val <- sim_out$avg_travel_time
      time_samples[r] <- ifelse(is.na(t_val), NA, t_val) 
    }
    
    results_list[[i]] <- data.frame(
      p = p_val,
      rho_mean = mean(rho_samples),
      rho_sd = sd(rho_samples),
      avg_travel_time = mean(time_samples, na.rm=TRUE), 
      topology = config$type,
      config_name = config$name,
      h_param = config$h
    )
    setTxtProgressBar(pb, i)
  }
  close(pb)
  df_res <- do.call(rbind, results_list)

  #automatic p_c detection
  crit_row <- head(subset(df_res, rho_mean > THRESHOLD), 1) 
  if (nrow(crit_row) > 0){
    pc_est <- crit_row$p 
  } else{ 
    pc_est <- max(p_seq)
    warning(paste("Config", config$name, "did not jam!"))
  }
  df_res$p_critical <- pc_est
  
  return(df_res)
}

run_full_experiment <- function(){
  all_results <- data.frame()
  for (cfg in experiment_configs){
    res <- run_averaged_sweep(cfg, P_RANGE)
    all_results <- rbind(all_results, res)
  }
  save(all_results, file = "simulation_results_FULL_MERGED.RData")
  cat("\nExperiment complete. Data saved to 'simulation_results_FULL_MERGED.RData'.\n")
  return(all_results)
}

final_data <- run_full_experiment()


# -----------------------------------------------------------------------------
# STEP 4: PLOTTING
# -----------------------------------------------------------------------------

if(!file.exists("simulation_results_FULL_MERGED.RData")){
  stop("Data files not found! Run the merger script first.")
}
load("simulation_results_FULL_MERGED.RData")

final_data$base_config <- sub("_h[0-9.]+$", "", final_data$config_name) 
final_data$h_factor <- as.factor(final_data$h_param)

#handle cases where p_critical might be NA (non-jamming) (for p/pc)
final_data$p_critical[is.na(final_data$p_critical)] <- max(final_data$p)
final_data$p_normalized <- final_data$p / final_data$p_critical

big_text_theme <- theme(
  axis.text       = element_text(size = 16, color = "black"),
  axis.title      = element_text(size = 18),
  strip.text      = element_text(size = 14),
  legend.text     = element_text(size = 16),
  legend.title    = element_text(size = 16),
  plot.title      = element_text(size = 18),
  plot.subtitle   = element_text(size = 16)
)


#Fixed topology config, varying h

#trees
df_tree <- subset(final_data, topology=="Tree")
p1_topo_fixed <- ggplot(df_tree, aes(x=p, y=rho_mean, color=h_factor)) + 
  geom_line(linewidth=0.7, alpha=0.5) + geom_point(size=2) +
  geom_errorbar(aes(ymin=rho_mean-rho_sd, ymax=rho_mean+rho_sd), width=0.1) +
  facet_wrap(~base_config, scales="free_x") + 
  scale_x_log10() + theme_bw() + big_text_theme +
  labs(title="Tree Networks: Effect of Routing Strategy (h)", 
       subtitle="Panels = Topology Config | Color = h parameter",
       x="p/p_c", y="rho", color="h")

#SF
df_sf <- subset(final_data, topology=="Scale-Free")
p2_topo_fixed <- ggplot(df_sf, aes(x=p, y=rho_mean, color=h_factor)) + 
  geom_line(linewidth=0.7, alpha=0.5) + geom_point(size=2) +
  geom_errorbar(aes(ymin=rho_mean-rho_sd, ymax=rho_mean+rho_sd), width=0.1) +
  facet_wrap(~base_config, scales="free_x") + 
  scale_x_log10() + theme_bw() + big_text_theme +
  labs(title="Scale-Free Networks: Effect of Routing Strategy (h)", 
       subtitle="Panels = Topology Config | Color = h parameter",
       x="p", y="rho", color="h")

#random
df_rnd <- subset(final_data, topology=="Random")
p3_topo_fixed <- ggplot(df_rnd, aes(x=p, y=rho_mean, color=h_factor)) + 
  geom_line(linewidth=1) + geom_point(size=2) +
  geom_errorbar(aes(ymin=rho_mean-rho_sd, ymax=rho_mean+rho_sd), width=0.1) +
  facet_wrap(~base_config, scales="free_x") + 
  scale_x_log10() + theme_bw() + big_text_theme + 
  labs(title="Random Networks: Effect of Routing Strategy (h)", 
       subtitle="Panels = Topology Config | Color = h parameter",
       x="p/p_c", y="rho", color="h")



#Fixed h, varying topology
#trees (h fixed)
p1_h_fixed <- ggplot(df_tree, aes(x=p, y=rho_mean, color=base_config)) + 
  geom_line(linewidth=0.7,  alpha=0.5) + geom_point(size=2) + 
  geom_errorbar(aes(ymin=rho_mean-rho_sd, ymax=rho_mean+rho_sd), width=0.1) +
  facet_wrap(~h_factor, scales="free_x") + 
  scale_x_log10() + theme_bw() + big_text_theme +
  labs(title="Tree Networks: Effect of Topology", 
       subtitle="Panels = h parameter | Color = Topology Config",
       x="p/p_c", y="rho", y="rho", color="Config")

# SF (h fixed)
p2_h_fixed <- ggplot(df_sf, aes(x=p_normalized, y=rho_mean, color=base_config)) + 
  geom_line(linewidth=0.7, alpha=0.5) + geom_point(size=2) +
  geom_errorbar(aes(ymin=rho_mean-rho_sd, ymax=rho_mean+rho_sd), width=0.1) +
  facet_wrap(~h_factor, scales="free_x") + 
  scale_x_log10() +
  theme_bw() + big_text_theme + # <--- Added Custom Theme
  labs(title="Scale-Free Networks: Effect of Topology", 
       subtitle="Panels = h parameter | Color = Topology Config",
       x="p/p_c", y="rho", color="Config")

#random (h fixed)
p3_h_fixed <- ggplot(df_rnd, aes(x=p, y=rho_mean, color=base_config)) + 
  geom_line(linewidth=0.7, alpha=0.5) + geom_point(size=2) +
  geom_errorbar(aes(ymin=rho_mean-rho_sd, ymax=rho_mean+rho_sd), width=0.1) +
  facet_wrap(~h_factor, scales="free_x") + 
  theme_bw() + big_text_theme + scale_x_log10() +
  labs(title="Random Networks: Effect of Topology", 
       subtitle="Panels = h parameter | Color = Topology Config",
       x="p/p_c", y="rho", color="Config")


#Efficiency analysis
p_target <- 3
df_eff <- subset(final_data, p==p_target)
df_eff <- df_eff %>% group_by(topology) %>% 
          mutate(variant_id = as.factor(as.numeric(as.factor(base_config)))) %>% ungroup()

p4 <- ggplot(df_eff, aes(x=base_config, y=avg_travel_time, color=h_factor, shape=variant_id)) +
  geom_point(size=5, stroke=1.5) + theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  big_text_theme +
  labs(title = paste0("Efficiency analysis (at p=", p_target, ")"), 
       subtitle = "Average travel time (T) vs topology configuration",
       x="Topology Configuration", y="Avg Travel Time (T)", color="h", shape="Config ID")


ggsave("fig1_tree_topo_fixed.png", plot=p1_topo_fixed, width=10, height=6)
ggsave("fig2_sf_topo_fixed.png", plot=p2_topo_fixed, width=10, height=6)
ggsave("fig3_rnd_topo_fixed.png", plot=p3_topo_fixed, width=10, height=6)

ggsave("fig4_tree_h_fixed.png", plot=p1_h_fixed, width=10, height=6)
ggsave("fig5_sf_h_fixed.png", plot=p2_h_fixed, width=10, height=6)
ggsave("fig6_rnd_h_fixed.png", plot=p3_h_fixed, width=10, height=6)

ggsave("fig7_efficiency.png", plot=p4, width=10, height=6)

print(p1_h_fixed)
print(p2_h_fixed)
print(p3_h_fixed)
print(p1_topo_fixed)
print(p2_topo_fixed)
print(p3_topo_fixed)
print(p4)


#replace here the configuration, of each type of network, that turns out to be the best
get_microscopic_data <- function(
    tree_conf=list(z=3, m=3), 
    sf_conf=list(n=400, m=4), 
    rnd_conf=list(n=200, p_link=0.02)){
      cat("\nRunning microscopic analysis:\n")

      g_tree <- generate_tree_topology(branching_factor=tree_conf$z, levels=tree_conf$m)
      g_sf <- generate_scalefree_topology(n_nodes=sf_conf$n, m_links=sf_conf$m)
      g_rnd <- generate_random_topology(n_nodes=rnd_conf$n) 
      
      graphs <- list(Tree=g_tree, ScaleFree=g_sf, Random=g_rnd)
      p_rate_fixed <- 10 #to ensure conjestion
      h_values <- c(1.0, 0.8, 0.3)
      all_micro_stats <- data.frame()
      for (g_type in names(graphs)){
        g_curr <- graphs[[g_type]]
        for (h_val in h_values){
          cat(paste("  Simulating:", g_type, "h =", h_val, "\n"))
          res <- run_simulation(g_curr, p_rate=p_rate_fixed, h_param=h_val, n_steps=1000)

          df <- res$micro_stats
          df$network_type <- g_type
          df$h_value <- as.factor(h_val)
          df$strategy <- paste0(g_type, " (h=", h_val, ")")
          all_micro_stats <- rbind(all_micro_stats, df)
        }
      }
      save(all_micro_stats, file = "microscopic_results_7.RData")
      return(all_micro_stats)
}

micro_data <- get_microscopic_data()

# --- Plot 5: microscopic analysis ---
if(!file.exists("microscopic_results_7.RData")){
  stop("Data files not found! Run 'get_microscopic_data()' first.")
}
load("microscopic_results_7.RData")

p5 <- ggplot(micro_all, aes(x=betweenness, y=queue_size, color=h_value)) +
  geom_point(alpha=0.6, size=2) + scale_x_log10() + scale_y_log10() +
  theme_bw() + facet_grid(network_type ~ h_value) +
  labs(title="Microscopic congestion profile", subtitle="Queue size vs betweenness (Rows: network, Cols: h)",
       x="Betweenness", y="Queue size") + theme(legend.position="none")

# some statitics to better interpret plot 5
micro_stats_summary <- micro_all %>% group_by(network_type, h_value) %>%
  summarize(
    Total_Packets=sum(queue_size),
    Weighted_Avg_Queue=sum(queue_size*betweenness)/sum(betweenness),
    Congestion_Center_Mass=sum(queue_size*betweenness)/sum(queue_size),
    Correlation_Spearman=cor(queue_size, betweenness, method="spearman"),
    .groups = 'drop'
  )
cat("\nMICROSCOPIC STATISTICS SUMMARY:\n")
as.data.frame(micro_stats_summary)

ggsave("figure5_microscopic.png", plot=p5, width=10, height=8)
print(p5)
