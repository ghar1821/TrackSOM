library(testthat)

# how to run on rstudio
# setwd to the directory of this file
# then type test_dir(".")

source("../TrackingFunctions.R", chdir = TRUE)
source("test_helpers.R")


#' Test new meta clusters only experience slight movement in code
test_that("small movement", {
  
  n_codes <- c(4, 4)
  timesteps <- 2
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,1,2,2))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,1,2,2))
  
  fsom.obj[['metaclustering']][['metaclustering']] <- list()
  fsom.obj[['metaclustering']][['mcPerCell']] <- list()
  
  for (i in seq_len(timesteps)) {
    label <- paste('metacluster', i, sep = '_')
    fsom.obj[['metaclustering']][['metaclustering']][[label]] <- meta_clusters_per_cell[[label]]
    fsom.obj[['metaclustering']][['mcPerCell']][[label]] <- meta_clusters_per_cell[[label]]
  }
  
  # setup the code's centroid
  mat <- matrix(1:(n_codes[1]*2), nrow=n_codes[1], ncol=2)
  mat[1,1] <- 0.1
  mat[1,2] <- 0.1
  mat[2,1] <- 0.2
  mat[2,2] <- 0.2
  mat[3,1] <- 0.8
  mat[3,2] <- 0.8
  mat[4,1] <- 0.9
  mat[4,2] <- 0.9
  fsom.obj[['FlowSOM']][['map']][['coding']][[1]] <- mat
  
  mat <- matrix(1:(n_codes[2]*2), nrow=n_codes[2], ncol=2)
  mat[1,1] <- 0.11
  mat[1,2] <- 0.11
  mat[2,1] <- 0.21
  mat[2,2] <- 0.21
  mat[3,1] <- 0.81
  mat[3,2] <- 0.81
  mat[4,1] <- 0.91
  mat[4,2] <- 0.91
  fsom.obj[['FlowSOM']][['map']][['coding']][[2]] <- mat
  
  # run tracking
  fsom.obj$tracking$lineage <- TrackingByLineage(fsom.obj)
  fsom.obj$tracking$proximity <- TrackingByHistoricalProximity(fsom.obj, fsom.obj$tracking$lineage)
  
  # test
  expected_proximity_ids <- list()
  expected_proximity_ids[[2]] <- c("A", "A", "B", "B")
  test_proximities(proximity_ids = fsom.obj[['tracking']][['proximity']],
                   expected_proximity_ids = expected_proximity_ids,
                   timesteps = 2)
})

#' Test new meta clusters close to cluster A
test_that("new cluster close to A", {
  
  n_codes <- c(4, 6)
  timesteps <- 2
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,1,2,2,NA,NA))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,1,2,2,3,3))
  
  fsom.obj[['metaclustering']][['metaclustering']] <- list()
  fsom.obj[['metaclustering']][['mcPerCell']] <- list()
  
  for (i in seq_len(timesteps)) {
    label <- paste('metacluster', i, sep = '_')
    fsom.obj[['metaclustering']][['metaclustering']][[label]] <- meta_clusters_per_cell[[label]]
    fsom.obj[['metaclustering']][['mcPerCell']][[label]] <- meta_clusters_per_cell[[label]]
  }
  
  # setup the code's centroid
  nrow <- max(n_codes)
  ncol <- 2
  mat <- matrix(1:(nrow*ncol), nrow=nrow, ncol=ncol)
  mat[1,1] <- 0.1
  mat[1,2] <- 0.1
  mat[2,1] <- 0.2
  mat[2,2] <- 0.2
  mat[3,1] <- 0.8
  mat[3,2] <- 0.8
  mat[4,1] <- 0.9
  mat[4,2] <- 0.9
  mat[5,1] <- NA
  mat[5,2] <- NA
  mat[6,1] <- NA
  mat[6,2] <- NA
  fsom.obj[['FlowSOM']][['map']][['coding']][[1]] <- mat
  
  mat <- matrix(1:(nrow*ncol), nrow=nrow, ncol=ncol)
  mat[1,1] <- 0.11
  mat[1,2] <- 0.11
  mat[2,1] <- 0.21
  mat[2,2] <- 0.21
  mat[3,1] <- 0.81
  mat[3,2] <- 0.81
  mat[4,1] <- 0.91
  mat[4,2] <- 0.91
  mat[5,1] <- 0.12
  mat[5,2] <- 0.12
  mat[6,1] <- 0.22
  mat[6,2] <- 0.22
  fsom.obj[['FlowSOM']][['map']][['coding']][[2]] <- mat
  
  # run tracking
  fsom.obj$tracking$lineage <- TrackingByLineage(fsom.obj)
  fsom.obj$tracking$proximity <- TrackingByHistoricalProximity(fsom.obj, fsom.obj$tracking$lineage)
  
  # test
  expected_proximity_ids <- list()
  expected_proximity_ids[[2]] <- c("A", "A", "B", "B", "A", "A")
  test_proximities(proximity_ids = fsom.obj[['tracking']][['proximity']],
                   expected_proximity_ids = expected_proximity_ids,
                   timesteps = 2)
})

#' Test split cluster
test_that("split cluster", {
  
  n_codes <- c(4, 4)
  timesteps <- 2
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,1,2,2))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,1,2,3))
  
  fsom.obj[['metaclustering']][['metaclustering']] <- list()
  fsom.obj[['metaclustering']][['mcPerCell']] <- list()
  
  for (i in seq_len(timesteps)) {
    label <- paste('metacluster', i, sep = '_')
    fsom.obj[['metaclustering']][['metaclustering']][[label]] <- meta_clusters_per_cell[[label]]
    fsom.obj[['metaclustering']][['mcPerCell']][[label]] <- meta_clusters_per_cell[[label]]
  }
  
  # setup the code's centroid
  nrow <- max(n_codes)
  ncol <- 2
  mat <- matrix(1:(nrow*ncol), nrow=nrow, ncol=ncol)
  mat[1,1] <- 0.1
  mat[1,2] <- 0.1
  mat[2,1] <- 0.2
  mat[2,2] <- 0.2
  mat[3,1] <- 0.8
  mat[3,2] <- 0.8
  mat[4,1] <- 0.9
  mat[4,2] <- 0.9
  fsom.obj[['FlowSOM']][['map']][['coding']][[1]] <- mat
  
  mat <- matrix(1:(nrow*ncol), nrow=nrow, ncol=ncol)
  mat[1,1] <- 0.11
  mat[1,2] <- 0.11
  mat[2,1] <- 0.21
  mat[2,2] <- 0.21
  mat[3,1] <- 0.81
  mat[3,2] <- 0.81
  mat[4,1] <- 0.99
  mat[4,2] <- 0.99
  fsom.obj[['FlowSOM']][['map']][['coding']][[2]] <- mat
  
  # run tracking
  fsom.obj$tracking$lineage <- TrackingByLineage(fsom.obj)
  fsom.obj$tracking$proximity <- TrackingByHistoricalProximity(fsom.obj, fsom.obj$tracking$lineage)
  
  # test
  expected_proximity_ids <- list()
  expected_proximity_ids[[2]] <- c("A", "A", "B", "B")
  test_proximities(proximity_ids = fsom.obj[['tracking']][['proximity']],
                   expected_proximity_ids = expected_proximity_ids,
                   timesteps = 2)
})

#' Test merge cluster
test_that("merged cluster", {
  
  n_codes <- c(6, 6)
  timesteps <- 2
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,1,2,2,3,3))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,1,1,2,3,3))
  
  fsom.obj[['metaclustering']][['metaclustering']] <- list()
  fsom.obj[['metaclustering']][['mcPerCell']] <- list()
  
  for (i in seq_len(timesteps)) {
    label <- paste('metacluster', i, sep = '_')
    fsom.obj[['metaclustering']][['metaclustering']][[label]] <- meta_clusters_per_cell[[label]]
    fsom.obj[['metaclustering']][['mcPerCell']][[label]] <- meta_clusters_per_cell[[label]]
  }
  
  # setup the code's centroid
  nrow <- max(n_codes)
  ncol <- 2
  mat <- matrix(1:(nrow*ncol), nrow=nrow, ncol=ncol)
  mat[1,1] <- 0.1
  mat[1,2] <- 0.1
  mat[2,1] <- 0.2
  mat[2,2] <- 0.2
  mat[3,1] <- 0.3
  mat[3,2] <- 0.3
  mat[4,1] <- 0.4
  mat[4,2] <- 0.4
  mat[5,1] <- 0.8
  mat[5,2] <- 0.8
  mat[6,1] <- 0.9
  mat[6,2] <- 0.9
  fsom.obj[['FlowSOM']][['map']][['coding']][[1]] <- mat
  
  mat <- matrix(1:(nrow*ncol), nrow=nrow, ncol=ncol)
  mat[1,1] <- 0.11
  mat[1,2] <- 0.11
  mat[2,1] <- 0.21
  mat[2,2] <- 0.21
  mat[3,1] <- 0.31
  mat[3,2] <- 0.31
  mat[4,1] <- 0.49
  mat[4,2] <- 0.49
  mat[5,1] <- 0.81
  mat[5,2] <- 0.81
  mat[6,1] <- 0.91
  mat[6,2] <- 0.91
  fsom.obj[['FlowSOM']][['map']][['coding']][[2]] <- mat
  
  # run tracking
  fsom.obj$tracking$lineage <- TrackingByLineage(fsom.obj)
  fsom.obj$tracking$proximity <- TrackingByHistoricalProximity(fsom.obj, fsom.obj$tracking$lineage)
  
  # test
  expected_proximity_ids <- list()
  expected_proximity_ids[[2]] <- c("A & B", "A & B", "A & B", "B", "C", "C")
  test_proximities(proximity_ids = fsom.obj[['tracking']][['proximity']],
                   expected_proximity_ids = expected_proximity_ids,
                   timesteps = 2)
})
