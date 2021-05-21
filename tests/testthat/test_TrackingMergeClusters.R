context("Tracking merged clusters for merging allowed and disallowed")
library(TrackSOM)

#' Test simple merge
test_that("simple merge", {

  n_codes <- rep(5, 3)
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,2,3,3,4))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,2,2,2,3))
  meta_clusters_per_cell[['metacluster_3']] <- as.factor(c(1,2,2,2,2))

  fsom.obj[['metaclustering']][['metaclustering']] <- list()
  fsom.obj[['metaclustering']][['mcPerCell']] <- list()

  for (i in seq_len(timesteps)) {
    label <- paste('metacluster', i, sep = '_')
    fsom.obj[['metaclustering']][['metaclustering']][[label]] <- meta_clusters_per_cell[[label]]
    fsom.obj[['metaclustering']][['mcPerCell']][[label]] <- meta_clusters_per_cell[[label]]
  }

  # run tracking
  fsom.obj$tracking$lineage <- TrackingWithMerging(fsom.obj)

  # test
  expected_lineage_ids <- list()
  expected_lineage_ids[['timestep_1']] <- c("A", "B", "C", "C", "D")
  expected_lineage_ids[['timestep_2']] <- c("A", "(B,C)", "(B,C)", "(B,C)", "D")
  expected_lineage_ids[['timestep_3']] <- c("A", "((B,C),D)", "((B,C),D)", "((B,C),D)", "((B,C),D)")

  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = timesteps)
})

test_that("merge not allowed", {

  n_codes <- rep(5, 3)
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,2,3,3,4))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,2,2,2,3))
  meta_clusters_per_cell[['metacluster_3']] <- as.factor(c(1,2,2,2,2))

  fsom.obj[['metaclustering']][['metaclustering']] <- list()
  fsom.obj[['metaclustering']][['mcPerCell']] <- list()

  for (i in seq_len(timesteps)) {
    label <- paste('metacluster', i, sep = '_')
    fsom.obj[['metaclustering']][['metaclustering']][[label]] <- meta_clusters_per_cell[[label]]
    fsom.obj[['metaclustering']][['mcPerCell']][[label]] <- meta_clusters_per_cell[[label]]
  }

  # run tracking
  fsom.obj$tracking$lineage <- TrackingNoMerging(fsom.obj)

  # test
  expected_lineage_ids <- list()
  expected_lineage_ids[['timestep_1']] <- c("A", "B", "C", "C", "D")
  expected_lineage_ids[['timestep_2']] <- c("A", "B", "C", "C", "D")
  expected_lineage_ids[['timestep_3']] <- c("A", "B", "C", "C", "D")

  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = timesteps)
})

test_that("merge not allowed but split cluster", {

  n_codes <- rep(5, 3)
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,2,3,3,4))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,2,2,3,3))
  meta_clusters_per_cell[['metacluster_3']] <- as.factor(c(1,2,2,2,2))

  fsom.obj[['metaclustering']][['metaclustering']] <- list()
  fsom.obj[['metaclustering']][['mcPerCell']] <- list()

  for (i in seq_len(timesteps)) {
    label <- paste('metacluster', i, sep = '_')
    fsom.obj[['metaclustering']][['metaclustering']][[label]] <- meta_clusters_per_cell[[label]]
    fsom.obj[['metaclustering']][['mcPerCell']][[label]] <- meta_clusters_per_cell[[label]]
  }

  # run tracking
  fsom.obj$tracking$lineage <- TrackingNoMerging(fsom.obj)

  # test
  expected_lineage_ids <- list()
  expected_lineage_ids[['timestep_1']] <- c("A", "B", "C", "C", "D")
  expected_lineage_ids[['timestep_2']] <- c("A", "B", "C", "C|1", "D")
  expected_lineage_ids[['timestep_3']] <- c("A", "B", "C", "C|1", "D")

  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = timesteps)
})

test_that("merge not allowed but split and new cluster", {

  n_codes <- rep(6, 3)
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,2,3,3,4,NA))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,2,2,3,3,4))
  meta_clusters_per_cell[['metacluster_3']] <- as.factor(c(1,2,2,2,2,3))

  fsom.obj[['metaclustering']][['metaclustering']] <- list()
  fsom.obj[['metaclustering']][['mcPerCell']] <- list()

  for (i in seq_len(timesteps)) {
    label <- paste('metacluster', i, sep = '_')
    fsom.obj[['metaclustering']][['metaclustering']][[label]] <- meta_clusters_per_cell[[label]]
    fsom.obj[['metaclustering']][['mcPerCell']][[label]] <- meta_clusters_per_cell[[label]]
  }

  # run tracking
  fsom.obj$tracking$lineage <- TrackingNoMerging(fsom.obj)

  # test
  expected_lineage_ids <- list()
  expected_lineage_ids[['timestep_1']] <- c("A", "B", "C", "C", "D", NA)
  expected_lineage_ids[['timestep_2']] <- c("A", "B", "C", "C|1", "D", "E")
  expected_lineage_ids[['timestep_3']] <- c("A", "B", "C", "C|1", "D", "E")

  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = timesteps)
})

test_that("merge not allowed but more node in old cluster", {

  n_codes <- rep(6, 3)
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,2,3,3,4,NA))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,2,2,3,3,1))
  meta_clusters_per_cell[['metacluster_3']] <- as.factor(c(1,2,2,2,2,1))

  fsom.obj[['metaclustering']][['metaclustering']] <- list()
  fsom.obj[['metaclustering']][['mcPerCell']] <- list()

  for (i in seq_len(timesteps)) {
    label <- paste('metacluster', i, sep = '_')
    fsom.obj[['metaclustering']][['metaclustering']][[label]] <- meta_clusters_per_cell[[label]]
    fsom.obj[['metaclustering']][['mcPerCell']][[label]] <- meta_clusters_per_cell[[label]]
  }

  # run tracking
  fsom.obj$tracking$lineage <- TrackingNoMerging(fsom.obj)

  # test
  expected_lineage_ids <- list()
  expected_lineage_ids[['timestep_1']] <- c("A", "B", "C", "C", "D", NA)
  expected_lineage_ids[['timestep_2']] <- c("A", "B", "C", "C|1", "D", "A")
  expected_lineage_ids[['timestep_3']] <- c("A", "B", "C", "C|1", "D", "A")

  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = timesteps)
})

test_that("merge not allowed and nodes going their own way", {

  n_codes <- rep(4, 2)
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
  meta_clusters_per_cell[['metacluster_1']] <- as.factor(c(1,2,3,3))
  meta_clusters_per_cell[['metacluster_2']] <- as.factor(c(1,1,1,2))

  fsom.obj[['metaclustering']][['metaclustering']] <- list()
  fsom.obj[['metaclustering']][['mcPerCell']] <- list()

  for (i in seq_len(timesteps)) {
    label <- paste('metacluster', i, sep = '_')
    fsom.obj[['metaclustering']][['metaclustering']][[label]] <- meta_clusters_per_cell[[label]]
    fsom.obj[['metaclustering']][['mcPerCell']][[label]] <- meta_clusters_per_cell[[label]]
  }

  # run tracking
  fsom.obj$tracking$lineage <- TrackingNoMerging(fsom.obj)

  # test
  expected_lineage_ids <- list()
  expected_lineage_ids[['timestep_1']] <- c("A", "B", "C", "C")
  expected_lineage_ids[['timestep_2']] <- c("A", "B", "C", "C|1")

  test_lineages(lineage_ids = fsom.obj[['tracking']][['lineage']],
                expected_lineage_ids = expected_lineage_ids,
                timesteps = timesteps)
})

