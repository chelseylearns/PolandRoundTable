

#------------------------------------------------------------------------------
# Function for getting cross-betweenness for a given graph.
# Uses affiliation_dates (filtered to start and end) to create graph
# on_cols determines what categories create edges


get_cross_betweenness <- function(affils_by_date, on_cols, start, end) {
  # get edge list
  edgelist <- get_edgelist_members(
    affils_by_date,
    on_cols = on_cols,
    start = start,
    end = end,
    # weight_col = "weight"
  )
  if(is.null(edgelist) || nrow(edgelist) == 0) {
    # return 0 or NA
    return(tibble(Member.ID = character(), CrossBetweenness = numeric()))
  }

  # build the graph
  graph <- igraph::graph_from_data_frame(edgelist, directed = FALSE)

  # add affiliation as vertex attribute
  V(graph)$RT.Affiliation <- purrr::map_chr(
    V(graph)$name,
    ~ member_meta_info$RT.Affiliation[member_meta_info$Member.ID == .x]
  )

  # # create sets of Gov and Opp nodes for pairing
  # gov_idx <- which(V(graph)$RT.Affiliation == "Government")
  # opp_idx <- which(V(graph)$RT.Affiliation == "Opposition")
  #
  # if (length(gov_idx) == 0 || length(opp_idx) == 0) {
  #   cross_betweenness <- rep(0, igraph::vcount(graph))
  # } else {
  #   pairs_mat <- as.matrix(expand.grid(gov_idx, opp_idx))
  #   cross_betweenness <- igraph::betweenness(graph,
  #                                            directed = FALSE,
  #                                            pairs = pairs_mat,
  #                                            weights = E(graph)$weight)
  # }

  cross_betweenness <- cross_betweenness(graph)

  tibble(
    Member.ID = igraph::V(graph)$name,
    Cross.Betweenness = cross_betweenness
  )
}


#-----------------------------------------------------------------------------
# Calculate cross-betweenness
#
# Betweenness, but only considering shortest paths for node pairs from
# opposing political parties (Gov & Opp)


cross_betweenness <- function(graph) {
  # make groups of government nodes and opposition node
  gov_idx <- which(V(graph)$RT.Affiliation == "Government")
  opp_idx <- which(V(graph)$RT.Affiliation == "Opposition")

  # if none, then just return zeros
  if (length(gov_idx) == 0 || length(opp_idx) == 0) {
    return(rep(0, igraph::vcount(graph)))
  }

  scores <- numeric(igraph::vcount(graph))

  # loop over each pair gov-opp
  for (g in gov_idx) {
    for (o in opp_idx) {
      # get all shortest paths from g to o
      # ignore weights so that igraph doesn't do djikstra
      short_paths <- igraph::all_shortest_paths(graph, from = g, to = o, weights = NA)$res
      num_paths <- length(short_paths)
      if(num_paths == 0) {
        # skip if there is no shortest path from g to o
        next
      }

      # for each shortest path, add 1/num_paths to interior
      for (path_vec in short_paths) {
        if (length(path_vec) <= 2) {
          # no interior nodes
          next
        }
        # exclude first and last node
        interior <- path_vec[2:(length(path_vec) - 1)]
        scores[interior] <- scores[interior] + (1 / num_paths)
      }
    }
  }

  return(scores)

}

#-------------------------------------------------------------------------------
# Calculate cross-betweenness for a graph, but normalized by the cluster sizes
# of the starting and target nodes
# Uses Louvain method of community detection

cross_betweenness_normalized <- function(graph) {

  # detect communities
  louvain_communities <- igraph::cluster_louvain(graph)

  # create a vector for which community each node belongs to
  membership_vector <- igraph::membership(louvain_communities)

  # get sizes of communities
  community_sizes <- table(membership_vector)

  cluster_size_vector <- as.numeric(
    community_sizes[as.character(membership_vector)])

  # make groups of government nodes and opposition node
  gov_idx <- which(V(graph)$RT.Affiliation == "Government")
  opp_idx <- which(V(graph)$RT.Affiliation == "Opposition")

  # if none, then just return zeros
  if (length(gov_idx) == 0 || length(opp_idx) == 0) {
    return(rep(0, igraph::vcount(graph)))
  }

  scores <- numeric(igraph::vcount(graph))

  # loop over each pair gov-opp
  for (g in gov_idx) {
    for (o in opp_idx) {

      # normalizing factor
      factor <- 1 / (cluster_size_vector[g] * cluster_size_vector[o])

      # get all shortest paths from g to o
      # ignore weights so that igraph doesn't do djikstra
      short_paths <- igraph::all_shortest_paths(
        graph, from = g, to = o, weights = NA)$res
      num_paths <- length(short_paths)
      if(num_paths == 0) {
        # skip if there is no shortest path from g to o
        next
      }

      # for each shortest path, add 1/num_paths to interior
      for (path_vec in short_paths) {
        if (length(path_vec) <= 2) {
          # no interior nodes
          next
        }
        # exclude first and last node
        interior <- path_vec[2:(length(path_vec) - 1)]
        scores[interior] <- scores[interior] + (factor / num_paths)
      }
    }
  }

  return(scores)
}

#------------------------------------------------------------------------------
# Function for getting cross-betweenness for a given graph, but normalized
# by cluster size of the start and target nodes.
# Uses affiliation_dates (filtered to start and end) to create graph
# on_cols determines what categories create edges


get_cross_betweenness_norm <- function(affils_by_date, on_cols, start, end) {
  # get edge list
  edgelist <- get_edgelist_members(
    affils_by_date,
    on_cols = on_cols,
    start = start,
    end = end,
    # weight_col = "weight"
  )
  if(is.null(edgelist) || nrow(edgelist) == 0) {
    # return 0 or NA
    return(tibble(Member.ID = character(), CrossBetweenness = numeric()))
  }

  # build the graph
  graph <- igraph::graph_from_data_frame(edgelist, directed = FALSE)

  # add affiliation as vertex attribute
  V(graph)$RT.Affiliation <- purrr::map_chr(
    V(graph)$name,
    ~ member_meta_info$RT.Affiliation[member_meta_info$Member.ID == .x]
  )

  cross_betweenness <- cross_betweenness_normalized(graph)

  tibble(
    Member.ID = igraph::V(graph)$name,
    Cross.Betweenness.Norm = cross_betweenness
  )
}


#-------------------------------------------------------------------------------
# Get cross_betweenness for all graphs between start and end, currently
# set to produce yearly graphs -- need to change this at some point (after
# cross_betweenness is optimized).

cross_betweenness_all <- function(affils_by_date, on_cols, start, end) {
  # all_starts <- seq.Date(start, end, by = "month")
  # all_ends <- all_starts + 30
  all_starts <- seq.Date(from = start, to = end, by = "year")
  all_ends <- all_starts %m+% years(1) - days(1)


  results <- purrr::map2_dfr(all_starts, all_ends, function(s, e) {
    temp <- get_cross_betweenness(
      affils_by_date,
      on_cols,
      start = s,
      end = e
    )
    temp |>
      mutate(Start.Date = s, End.Date = e)
  })

  results
}


#-------------------------------------------------------------------------------
# Get cross_betweenness (normalized by cluster size of the start and
# target nodes) for all graphs between start and end, currently
# set to produce yearly graphs -- need to change this at some point (after
# cross_betweenness is optimized).

cross_betweenness_norm_all <- function(affils_by_date, on_cols, start, end) {
  # all_starts <- seq.Date(start, end, by = "month")
  # all_ends <- all_starts + 30
  all_starts <- seq.Date(from = start, to = end, by = "year")
  all_ends <- all_starts %m+% years(1) - days(1)


  results <- purrr::map2_dfr(all_starts, all_ends, function(s, e) {
    temp <- get_cross_betweenness_norm(
      affils_by_date,
      on_cols,
      start = s,
      end = e
    )
    temp |>
      mutate(Start.Date = s, End.Date = e)
  })

  results
}

