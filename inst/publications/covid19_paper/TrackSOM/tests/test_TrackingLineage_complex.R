library(testthat)

# how to run on rstudio
# setwd to the directory of this file
# then type test_dir(".")

source("../TrackingFunctions.R", chdir = TRUE)
source("test_helpers.R")

#' Test merge then split
test_that("merge then split", {
  
  n_codes <- c(4,rep(5,2))
  timesteps <- 3
  fsom.obj <- construct_dummy_metadata(n_codes, timesteps)
  
  # add code per cell
  dummy_codes_list <- lapply(seq_len(timesteps), function(i) {
    temp_codes <- create_codes(n_codes[i])
    return(temp_codes)
  })
  
  dummy_codes <- do.call(rbind, dummy_codes_list)
  fsom.obj[['FlowSOM']][['map']][['mapping']] <- dummy_codes
  
  # add metaclustering details
  meta_clusters_per_cell <- list()
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,2,3,3))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,2,2,2,3))
  meta_clusters_per_cell[['metacluster_3']] <- as.factor(c(1,2,2,3,3))
  
  fsom.obj[['metaclustering']][['metaclustering']] <- list()
  fsom.obj[['metaclustering']][['mcPerCell']] <- list()
  
  for (i in seq_len(timesteps)) {
    label <- paste('metacluster', i, sep = '_')
    fsom.obj[['metaclustering']][['metaclustering']][[label]] <- meta_clusters_per_cell[[label]]
    fsom.obj[['metaclustering']][['mcPerCell']][[label]] <- meta_clusters_per_cell[[label]]
  }
  
  # run tracking
  fsom.obj$tracking$lineage <- TrackingByLineage(fsom.obj)
  
  # test
  expected_lineage_ids <- list()
  expected_lineage_ids[['timestep_1']] <- c("A", "B", "C", "C")
  expected_lineage_ids[['timestep_2']] <- c("A", "(B,C)", "(B,C)", "(B,C)", "D")
  expected_lineage_ids[['timestep_3']] <- c("A", "(B,C)", "(B,C)", "((B,C)|1,D)", "((B,C)|1,D)")
  
  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = timesteps)
})

#' Test complicated scenario where there are merge, split, new cluster, code not filled the next day, etc.
test_that("complicated merge, split, empty code, new code", {
  
  n_codes <- c(3,7,13,18)
  timesteps <- 4
  fsom.obj <- construct_dummy_metadata(n_codes, timesteps)
  
  # add code per cell
  # can't use create_codes as we want some codes to not be filled with any cells
  # we prespecify what code is filled with cells
  codes_per_timestep <- list()
  codes_per_timestep[[1]] <- c(1,2,4)
  codes_per_timestep[[2]] <- c(1,2,4,6,7,9,15)
  codes_per_timestep[[3]] <- c(1,2,4,9,15,24,25,
                               6,26,27,21,
                               7,23)
  codes_per_timestep[[4]] <- c(1, 2, 33, 4, 6, 39, 9, 15, 21, 24, 25,
                               34, 38, 7, 40, 23,
                               26, 27)
  dummy_codes_list <- lapply(seq_len(timesteps), function(i) {
    num_codes <- n_codes[i]
    num_row = num_codes
    num_col = 2
    num_cells = num_row * num_col
    
    mat <- matrix(1:num_cells, nrow=num_row, ncol=num_col)
    
    rnd_2nd_col_vals <- runif(num_codes)
    codes_id <- codes_per_timestep[[i]]
    
    for (idx in seq_len(num_codes)) {
      mat[idx, 1] <- codes_id[idx]
      mat[idx, 2] <- rnd_2nd_col_vals[idx]
    }
    
    return(mat)
  })
  
  dummy_codes <- do.call(rbind, dummy_codes_list)
  fsom.obj[['FlowSOM']][['map']][['mapping']] <- dummy_codes
  
  ## The following is very repetitive.
  ## TODO make it better.
  meta_clusters_per_cell <- list()
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(sapply(1:num_codes_total, function(i) {
    if (i %in% c(1,2,4)) {
      return(1)
    } else {
      return(NA)
    }
  }))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(sapply(1:num_codes_total, function(i) {
    if (i %in% c(1,2,4,6,7,9,15)) {
      return(1)
    } else {
      return(NA)
    }
  }))
  meta_clusters_per_cell[['metacluster_3']] <- as.factor(sapply(1:num_codes_total, function(i) {
    if (i %in% c(1,2,4,9,15,24,25)) {
      return(1)
    } else if (i %in% c(6,26,27,21)) {
      return(2)
    } else if (i %in% c(7,23)) {
      return(3)
    } else {
      return(NA)
    }
  }))
  meta_clusters_per_cell[['metacluster_4']] <- as.factor(sapply(1:num_codes_total, function(i) {
    if (i %in% c(1, 2, 33, 4, 6, 39, 9, 15, 21, 24, 25)) {
      return(1)
    } else if (i %in% c(34, 38, 7, 40, 23)) {
      return(2)
    } else if (i %in% c(26, 27)) {
      return(3)
    } else {
      return(NA)
    }
  }))
  
  fsom.obj[['metaclustering']][['metaclustering']] <- list()
  fsom.obj[['metaclustering']][['mcPerCell']] <- list()
  
  for (i in seq_len(timesteps)) {
    label <- paste('metacluster', i, sep = '_')
    fsom.obj[['metaclustering']][['metaclustering']][[label]] <- meta_clusters_per_cell[[label]]
    fsom.obj[['metaclustering']][['mcPerCell']][[label]] <- meta_clusters_per_cell[[label]]
  }
  
  # run tracking
  fsom.obj$tracking$lineage <- TrackingByLineage(fsom.obj)
  
  # test
  expected_lineage_ids <- list()
  expected_lineage_ids[['timestep_1']] <- sapply(1:num_codes_total, function(i) {
    if (i %in% c(1,2,4)) {
      return("A")
    } else {
      return(NA)
    }
  })
  expected_lineage_ids[['timestep_2']] <- sapply(1:num_codes_total, function(i) {
    if (i %in% c(1,2,4,6,7,9,15)) {
      return("A")
    } else {
      return(NA)
    }
  })
  expected_lineage_ids[['timestep_3']] <- sapply(1:num_codes_total, function(i) {
    if (i %in% c(1,2,4,9,15,24,25)) {
      return("A")
    } else if (i %in% c(6,26,27,21)) {
      return("A|1")
    } else if (i %in% c(7,23)) {
      return("A|2")
    } else {
      return(NA)
    }
  })
  expected_lineage_ids[['timestep_4']] <- sapply(1:num_codes_total, function(i) {
    if (i %in% c(1, 2, 33, 4, 6, 39, 9, 15, 21, 24, 25)) {
      return("(A,A|1)")
    } else if (i %in% c(34, 38, 7, 40, 23)) {
      return("A|2")
    } else if (i %in% c(26, 27)) {
      return("A|1|1")
    } else {
      return(NA)
    }
  })
  
  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = timesteps)
})