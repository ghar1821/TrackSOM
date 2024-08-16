test_that("collate results normal", {
  fsom.obj <- list()

  # code for each cell
  x <- matrix(1:12, nrow=6, ncol=2)
  x[1,1] <- 1
  x[1,2] <- 0.1
  x[2,1] <- 2
  x[2,2] <- 0.1
  x[3,1] <- 1
  x[3,2] <- 0.1

  x[4,1] <- 2
  x[4,2] <- 0.1
  x[5,1] <- 1
  x[5,2] <- 0.1
  x[6,1] <- 3
  x[6,2] <- 0.1
  fsom.obj[['FlowSOM']][['map']][['mapping']] <- x

  # meta cluster for each cell
  fsom.obj[['metaclustering']][['mcPerCell']][["metacluster_1"]] <- as.factor(c(1,2,1))
  fsom.obj[['metaclustering']][['mcPerCell']][["metacluster_2"]] <- as.factor(c(2,1,3))

  fsom.obj[['tracking']][['lineage']][['timestep_1']] <- as.factor(c("A", "B", NA))
  fsom.obj[['tracking']][['lineage']][['timestep_2']] <- as.factor(c("A", "B", "C"))

  cell.dat <- data.frame(x=c(1,2,3,4,5,6), y=c(7,8,5,3,6,3), day=c(1,1,1,2,2,2))
  res.dat <- ConcatenateClusteringDetails(tracksom.result = fsom.obj,
                                         dat = cell.dat,
                                         timepoint.col = "day",
                                         timepoints = c(1,2))
  # test here
  expect_equal(res.dat$x, c(1,2,3,4,5,6))
  expect_equal(res.dat$y, c(7,8,5,3,6,3))
  expect_equal(res.dat$day, c(1,1,1,2,2,2))
  expect_equal(res.dat$TrackSOM_metacluster, as.factor(c(1,2,1,2,1,3)))
  expect_equal(as.vector(res.dat$TrackSOM_metacluster_lineage_tracking),  c('A','B','A','B','A','C'))
})

test_that("collate results split", {
  fsom.obj <- list()

  # code for each cell
  x <- matrix(1:14, nrow=7, ncol=2)
  x[1,1] <- 1
  x[1,2] <- 0.1
  x[2,1] <- 2
  x[2,2] <- 0.1
  x[3,1] <- 3
  x[3,2] <- 0.1
  x[4,1] <- 4
  x[4,2] <- 0.1

  x[5,1] <- 1
  x[5,2] <- 0.1
  x[6,1] <- 3
  x[6,2] <- 0.1
  x[7,1] <- 4
  x[7,2] <- 0.1
  fsom.obj[['FlowSOM']][['map']][['mapping']] <- x

  # meta cluster for each cell
  fsom.obj[['metaclustering']][['mcPerCell']][["metacluster_1"]] <- as.factor(c(1,2,3,3))
  fsom.obj[['metaclustering']][['mcPerCell']][["metacluster_2"]] <- as.factor(c(1,2,3))

  # id per code
  fsom.obj[['tracking']][['lineage']][['timestep_1']] <- as.factor(c("A", "B", "C", "C"))
  # There is no cluster 2.
  fsom.obj[['tracking']][['lineage']][['timestep_2']] <- as.factor(c("A", NA, "C", "C|1"))

  cell.dat <- data.frame(x=c(1,2,3,4,5,6,7), y=c(7,8,5,3,6,3,7), day=c(1,1,1,1,2,2,2))
  res.dat <- ConcatenateClusteringDetails(tracksom.result = fsom.obj,
                                         dat = cell.dat,
                                         timepoint.col = "day",
                                         timepoints = c(1,2))
  # test here
  expect_equal(res.dat$x, c(1,2,3,4,5,6,7))
  expect_equal(res.dat$y, c(7,8,5,3,6,3,7))
  expect_equal(res.dat$day, c(1,1,1,1,2,2,2))
  expect_equal(res.dat$TrackSOM_metacluster, as.factor(c(1,2,3,3,1,2,3)))
  expect_equal(as.vector(res.dat$TrackSOM_metacluster_lineage_tracking),  c('A','B','C','C','A','C', 'C|1'))
})

test_that("collate results merge", {
  fsom.obj <- list()

  # code for each cell
  x <- matrix(1:14, nrow=7, ncol=2)
  x[1,1] <- 1
  x[1,2] <- 0.1
  x[2,1] <- 2
  x[2,2] <- 0.1
  x[3,1] <- 3
  x[3,2] <- 0.1
  x[4,1] <- 4
  x[4,2] <- 0.1

  x[5,1] <- 1
  x[5,2] <- 0.1
  x[6,1] <- 3
  x[6,2] <- 0.1
  x[7,1] <- 4
  x[7,2] <- 0.1
  fsom.obj[['FlowSOM']][['map']][['mapping']] <- x

  # meta cluster for each cell
  fsom.obj[['metaclustering']][['mcPerCell']][["metacluster_1"]] <- as.factor(c(1,2,3,4))
  fsom.obj[['metaclustering']][['mcPerCell']][["metacluster_2"]] <- as.factor(c(1,2,3))

  # id per code
  fsom.obj[['tracking']][['lineage']][['timestep_1']] <- as.factor(c("A", "B", "C", "D"))
  # There is no cluster 2.
  fsom.obj[['tracking']][['lineage']][['timestep_2']] <- as.factor(c("A", NA, "(C,D)", "(C,D)"))

  cell.dat <- data.frame(x=c(1,2,3,4,5,6,7), y=c(7,8,5,3,6,3,7), day=c(1,1,1,1,2,2,2))
  res.dat <- ConcatenateClusteringDetails(tracksom.result = fsom.obj,
                                         dat = cell.dat,
                                         timepoint.col = "day",
                                         timepoints = c(1,2))
  # test here
  expect_equal(res.dat$x, c(1,2,3,4,5,6,7))
  expect_equal(res.dat$y, c(7,8,5,3,6,3,7))
  expect_equal(res.dat$day, c(1,1,1,1,2,2,2))
  expect_equal(res.dat$TrackSOM_metacluster, as.factor(c(1,2,3,4,1,2,3)))
  expect_equal(as.vector(res.dat$TrackSOM_metacluster_lineage_tracking),  c('A','B','C','D','A','(C,D)','(C,D)'))
})
