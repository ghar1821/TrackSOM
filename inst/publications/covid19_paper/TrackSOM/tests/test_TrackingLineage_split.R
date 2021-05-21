library(testthat)

# how to run on rstudio
# setwd to the directory of this file
# then type test_dir(".")

source("../TrackingFunctions.R", chdir = TRUE)
source("test_helpers.R")


#' Test simple split
test_that("simple split", {
  
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,2,2,2))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,2,2,3))
  
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
  expected_lineage_ids[['timestep_1']] <- c("A", "B", "B", "B")
  expected_lineage_ids[['timestep_2']] <- c("A", "B", "B", "B|1")
  
  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = timesteps)
})

#' Test split clusters and new one as well
test_that("split and new", {
  
  n_codes <- c(4, 5)
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,2,2,2))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,2,2,3,4))
  
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
  expected_lineage_ids[['timestep_1']] <- c("A", "B", "B", "B")
  expected_lineage_ids[['timestep_2']] <- c("A", "B", "B", "B|1", "C")
  
  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = timesteps)
})

#' Test majority split is ordering insensitive
test_that("ids are allocated properly based on majority of codes regardless of the order in which codes/cells are stored", {
  
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,2,2,2))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,2,3,3))
  
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
  expected_lineage_ids[['timestep_1']] <- c("A", "B", "B", "B")
  expected_lineage_ids[['timestep_2']] <- c("A", "B|1", "B", "B")
  
  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = timesteps)
})

#' Test split into multipe clusters
test_that("split into more than 2 clusters", {
  
  n_codes <- c(5, 5)
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,2,2,2,2))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,2,3,3,4))
  
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
  expected_lineage_ids[['timestep_1']] <- c("A", "B", "B", "B","B")
  expected_lineage_ids[['timestep_2']] <- c("A", "B|1", "B", "B", "B|2")
  
  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = timesteps)
})

#' Test multi day split
test_that("split occurs over more than 2 days", {
  
  n_codes <- rep(4,3)
  timesteps <- 3
  # add metadata: number of rows/cells in each time point.
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,2,2,2))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,2,2,3))
  meta_clusters_per_cell[['metacluster_3']] <- as.factor(c(1,2,3,4))
  
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
  expected_lineage_ids[['timestep_1']] <- c("A", "B", "B", "B")
  expected_lineage_ids[['timestep_2']] <- c("A", "B", "B", "B|1")
  expected_lineage_ids[['timestep_3']] <- c("A", "B", "B|2", "B|1")
  
  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = timesteps)
})

#' Test complicated split
test_that("split occurs again and again", {
  
  n_codes <- rep(10,4)
  timesteps <- 4
  # add metadata: number of rows/cells in each time point.
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,2,rep(3,2),rep(4,6)))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,2,rep(3,2),rep(4,3),rep(5,3)))
  meta_clusters_per_cell[['metacluster_3']] <- as.factor(c(1,2,rep(3,2),rep(4,3),5,rep(6,2)))
  meta_clusters_per_cell[['metacluster_4']] <- as.factor(c(1,2,rep(3,2),rep(4,3),5,6,7))
  
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
  expected_lineage_ids[['timestep_1']] <- c("A", "B", rep("C",2), rep("D",6))
  expected_lineage_ids[['timestep_2']] <- c("A", "B", rep("C",2), rep("D",3), rep("D|1",3))
  expected_lineage_ids[['timestep_3']] <- c("A", "B", rep("C",2), rep("D",3), "D|1|1", rep("D|1",2))
  expected_lineage_ids[['timestep_4']] <- c("A", "B", rep("C",2), rep("D",3), "D|1|1", "D|1", "D|1|2")
  
  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = timesteps)
})

#' Test split only on last day
test_that("split only last day", {
  
  n_codes <- c(2,4,4)
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,2))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,rep(2,3)))
  meta_clusters_per_cell[['metacluster_3']] <- as.factor(c(1,2,3,4))
  
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
  expected_lineage_ids[['timestep_1']] <- c("A", "B")
  expected_lineage_ids[['timestep_2']] <- c("A", rep("B",3))
  expected_lineage_ids[['timestep_3']] <- c("A", "B", "B|1", "B|2")
  
  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = timesteps)
})