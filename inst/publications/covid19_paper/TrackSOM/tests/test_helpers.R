library(testthat)

#' Create a dummy matrix containing cells and codes.
#' This function is written to make it easier to get a data structure
#' storing dummy cells for each code, the core data structure of flowsom's map[[mapping]] object.
#' To simply testing, each code will only have 1 cell.
#' The object have 2 columns, 1 column indicate the code a cell belongs to, 2nd column, no idea. Not enough documentation about this.
#' 
#' @param num_codes. NO DEFAULT. Numeric. Number of codes you want to create
#' 
#' @return matrix of size 2xnum_codes where each row is a cell, and 1st column is the code the cell belongs to.
create_codes <- function(num_codes) {
  
  num_row = num_codes
  num_col = 2
  num_cells = num_row * num_col
  
  mat <- matrix(1:num_cells, nrow=num_row, ncol=num_col)
  
  rnd_2nd_col_vals <- runif(num_codes)
  codes_id <- seq_len(num_codes)
  
  for (code_id in codes_id) {
    mat[code_id, 1] <- code_id
    mat[code_id, 2] <- rnd_2nd_col_vals[code_id]
  }
  
  return(mat)
}

#' Test to make sure that the meta cluster id assigned by tracking by lineage matches what is expected
#'
#'@param lineage_ids NO DEFAULT. A list of containing the id assigned TrackSOM's tracking by lineage determination. This is produced by TrackSOM's tracking[["lineage"]]
#'@param expect_lineage_ids NO DEFAULT. A list containing the expected id. Each element is a vector of id for a time point.
#'@param timesteps NO DEFAULT. Numeric value of number of time points to match the ids.
test_lineages <- function(lineage_ids, expected_lineage_ids, timesteps) {
  for (i in seq_len(timesteps)) {
    label <- paste("timestep", i, sep = "_")
    lineages <- as.vector(lineage_ids[[label]])
    expected_lineages <- as.vector(expected_lineage_ids[[label]])
    expect_equal(length(lineages), length(expected_lineages))
    for (l in 1:length(lineages)) {
      expect_equal(lineages[l], expected_lineages[l])
    }  
  }
}

#' Construct a metadata which describe the limit of the cells for each day.
#' This metadata is part of TrackSOM's data structure which describe for each day, how many cells are there.
#' So if we have 100 cells per day, for day 1, it will be cell 1-100th, day 2 will be cell 101th-201th, and so on.
#' 
#' @param n_codes NO DEFAULT. Number of codes per time point.
#' @param timesteps NO DEFAULT. Number of time points.
construct_dummy_metadata <- function(n_codes, timesteps) {
  fsom.obj <- list()
  for (i in seq_len(timesteps)) {
    label <- paste("timestep", i, sep = "_")
    if (i == 1) {
      fsom.obj[['FlowSOM']][['metaData']][[label]] <- c(1,n_codes[i])
    } else {
      label_prev <- paste("timestep", i-1, sep = "_")
      limit_prev_day <- fsom.obj[['FlowSOM']][['metaData']][[label_prev]][2]
      fsom.obj[['FlowSOM']][['metaData']][[label]] <- c(limit_prev_day + 1, limit_prev_day + n_codes[i])
    }
  }
  return(fsom.obj)
}

#' Test to make sure that the meta cluster id assigned by tracking by lineage matches what is expected
#'
#'@param proximity_ids NO DEFAULT. A list of containing the id assigned TrackSOM's tracking by proximity. This is produced by TrackSOM's tracking[["lineage"]]
#'@param expected_proximity_ids NO DEFAULT. A list containing the expected id. Each element is a vector of id for a time point.
#'@param timesteps NO DEFAULT. Numeric value of number of time points to match the ids.
test_proximities <- function(proximity_ids, expected_proximity_ids, timesteps) {
  expect_null(proximity_ids[[1]])
  for (i in 2:timesteps) {
    proximities <- as.vector(proximity_ids[[i]])
    expected_proximities <- as.vector(expected_proximity_ids[[i]])
    expect_equal(length(proximities), length(expected_proximities))
    for (l in 1:length(proximities)) {
      expect_equal(proximities[l], expected_proximities[l])
    }  
  }
}
