library(testthat)

# how to run on rstudio
# setwd to the directory of this file
# then type test_dir(".")

source("../TrackingFunctions.R", chdir = TRUE)
source("test_helpers.R")


#' Test to make sure ids are assigned alphabetically
#' 4 codes, 3 meta clusters.
#' 2 meta clusters contain 1 code each, last meta cluster contains 2 codes.
test_that("cluster id assignment", {
  
  n_codes <- rep(4,2)
  timesteps <- 2
  ## Data for timestep 2 is not going to be checked.
  ## TODO fix the tracking code not allowing only 1 time point. It will break.
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
  fsom.obj[['metaclustering']][['metaclustering']] <- list()
  fsom.obj[['metaclustering']][['metaclustering']][['metacluster_1']] <- as.factor(c(1,2,3,3))
  fsom.obj[['metaclustering']][['metaclustering']][['metacluster_2']] <- as.factor(c(1,2,3,3))
  fsom.obj[['metaclustering']][['mcPerCell']] <- list()
  fsom.obj[['metaclustering']][['mcPerCell']][['metacluster_1']] <- as.factor(c(1,2,3,3))
  fsom.obj[['metaclustering']][['mcPerCell']][['metacluster_2']] <- as.factor(c(1,2,3,3))
  
  fsom.obj$tracking$lineage <- TrackingByLineage(fsom.obj)
  
  expected_lineage_ids <- list()
  expected_lineage_ids[['timestep_1']] <- c('A', 'B', 'C', 'C')
  
  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = 1)
})

#' Test to make sure ids are assigned correctly if we have more than 26 clusters.
test_that("cluster id assignment beyond 26 clusters", {
  
  n_codes <- rep(27,2)
  timesteps <- 2
  ## Data for timestep 2 is not going to be checked.
  ## TODO fix the tracking code not allowing only 1 time point. It will break.
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
  fsom.obj[['metaclustering']][['metaclustering']] <- list()
  fsom.obj[['metaclustering']][['metaclustering']][['metacluster_1']] <- as.factor(c(1:27))
  fsom.obj[['metaclustering']][['metaclustering']][['metacluster_2']] <- as.factor(c(1:27))
  fsom.obj[['metaclustering']][['mcPerCell']] <- list()
  fsom.obj[['metaclustering']][['mcPerCell']][['metacluster_1']] <- as.factor(c(1:27))
  fsom.obj[['metaclustering']][['mcPerCell']][['metacluster_2']] <- as.factor(c(1:27))
  
  fsom.obj$tracking$lineage <- TrackingByLineage(fsom.obj)
  
  expected_lineage_ids <- list()
  expected_lineage_ids[['timestep_1']] <- c(toupper(letters), "AA")
  
  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = 1)
})

#' Test new clusters are given new id.
#' 4 codes in day 1, 7 codes in day 2. All 4 codes in day 1 are also filled in day 2 and thus carried over.
#' 3 meta clusters in day 1, 5 in day 2, which means some codes are placed within same meta cluster.
test_that("new clusters are given new id", {
  
  n_codes <- c(4, 7)
  timesteps <- 2
  ## Data for timestep 2 is not going to be checked.
  ## TODO fix the tracking code not allowing only 1 time point. It will break.
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,2,3,3))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,2,3,3,4,5,5))
  
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
  expected_lineage_ids[['timestep_2']] <- c("A", "B", "C", "C", "D", "E", "E")
  
  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = 1)
  
})

#' Test new clusters are given new id even though they are assigned to cells regardless of whether they are assigned to cells at the end or start of the data.
#' 4 codes in day 1, 6 codes in day 2. All 4 codes in day 1 are also filled in day 2 and thus carried over.
#' 3 meta clusters in day 1, 4 in day 2, which means some codes are placed within same meta cluster.
test_that("new clusters are given new id regardless of cells ordering", {
  
  n_codes <- c(4, 6)
  timesteps <- 2
  ## Data for timestep 2 is not going to be checked.
  ## TODO fix the tracking code not allowing only 1 time point. It will break.
  # add metadata: number of rows/cells in each time point.
  fsom.obj <- construct_dummy_metadata(n_codes, timesteps)
  
  # add code per cell
  dummy_codes_list <- lapply(seq_len(timesteps), function(i) {
    temp_codes <- create_codes(n_codes[i])
    return(temp_codes)
  })
  
  dummy_codes <- do.call(rbind, dummy_codes_list)
  # shift cells belonging to code 5 and 6 up such that they come before those in code 3 and 4
  # the function gives results ordered by codes per day, so we can just explicitly shift by using row number
  dummy_codes <- rbind(dummy_codes[1:6,], dummy_codes[9:10,], dummy_codes[7:8,])
  
  
  fsom.obj[['FlowSOM']][['map']][['mapping']] <- dummy_codes
  
  # add metaclustering details
  meta_clusters_per_cell <- list()
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,2,3,3))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,2,4,4,3,3))
  
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
  expected_lineage_ids[['timestep_2']] <- c("A", "B", "D", "D", "C", "C")
  
  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = 1)
  
})