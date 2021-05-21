context("Test running TrackSOM")
library(TrackSOM)
library(data.table)

test_that("Merging works for prescribed invariant", {
    data_files <- sapply(c(0:4), function(i) {
      return(paste0("~/Documents/GitHub/TrackSOM/inst/extdata/synthetic_d", i,".csv"))
    })


    data_files <- sapply(c(0:4), function(i) {
        system.file("extdata", paste0("synthetic_d", i, ".csv"), package = "TrackSOM")
    })
    fsom <- TrackSOM(
        inputFiles = data_files,
        dataFileType = '.csv',
        colsToUse = c("x", "y", "z"),
        nClus = 10,
        tracking = TRUE,
        noMerge = FALSE
    )
    dat <- lapply(data_files, function(d) fread(d))
    dat <- rbindlist(dat)

    timepoints <- c('Mock', 'SYN-1', 'SYN-2', 'SYN-3', 'SYN-4')
    dat <- ConcatenateClusteringDetails(tracksom.result = fsom,
                                        dat = dat,
                                        timepoint.col = 'timepoint',
                                        timepoints = timepoints)

    dat_meta <- unique(dat[, c("timepoint", "TrackSOM_metacluster_lineage_tracking")])

    for (t in timepoints) {
        dat_t <- dat_meta[dat_meta$timepoint == t]
        expect_equal(nrow(dat_t), 10)
    }
})

test_that("Merging works for prescribed variant", {

    data_files <- sapply(c(0:4), function(i) {
        system.file("extdata", paste0("synthetic_d", i, ".csv"), package = "TrackSOM")
    })

    nclust <- c(10, 15, 20, 21, 30)
    fsom <- TrackSOM(
        inputFiles = data_files,
        dataFileType = '.csv',
        colsToUse = c("x", "y", "z"),
        nClus = nclust,
        tracking = TRUE,
        noMerge = FALSE
    )
    dat <- lapply(data_files, function(d) fread(d))
    dat <- rbindlist(dat)

    timepoints <- c('Mock', 'SYN-1', 'SYN-2', 'SYN-3', 'SYN-4')
    dat <- ConcatenateClusteringDetails(tracksom.result = fsom,
                                        dat = dat,
                                        timepoint.col = 'timepoint',
                                        timepoints = timepoints)

    dat_meta <- unique(dat[, c("timepoint", "TrackSOM_metacluster_lineage_tracking")])

    for (i in c(1:5)) {
        t <- timepoints[i]
        nclus <- nclust[i]
        dat_t <- dat_meta[dat_meta$timepoint == t]
        expect_equal(nrow(dat_t), nclus)
    }
})

test_that("Merging works for autonomous adaptive", {
    data_files <- sapply(c(0:4), function(i) {
        return(paste0("~/Documents/GitHub/TrackSOM/inst/extdata/synthetic_d", i,".csv"))
    })


    data_files <- sapply(c(0:4), function(i) {
        system.file("extdata", paste0("synthetic_d", i, ".csv"), package = "TrackSOM")
    })
    fsom <- TrackSOM(
        inputFiles = data_files,
        dataFileType = '.csv',
        colsToUse = c("x", "y", "z"),
        maxMeta = 20,
        tracking = TRUE,
        noMerge = FALSE
    )
    dat <- lapply(data_files, function(d) fread(d))
    dat <- rbindlist(dat)

    timepoints <- c('Mock', 'SYN-1', 'SYN-2', 'SYN-3', 'SYN-4')
    dat <- ConcatenateClusteringDetails(tracksom.result = fsom,
                                        dat = dat,
                                        timepoint.col = 'timepoint',
                                        timepoints = timepoints)

    dat_meta <- unique(dat[, c("timepoint", "TrackSOM_metacluster_lineage_tracking")])

    for (t in timepoints) {
        dat_t <- dat_meta[dat_meta$timepoint == t]
        # can't predict how many it will produce
        expect_gt(nrow(dat_t), 2)
    }
})

test_that("No merging works for prescribed invariant", {
    data_files <- sapply(c(0:4), function(i) {
        return(paste0("~/Documents/GitHub/TrackSOM/inst/extdata/synthetic_d", i,".csv"))
    })


    data_files <- sapply(c(0:4), function(i) {
        system.file("extdata", paste0("synthetic_d", i, ".csv"), package = "TrackSOM")
    })
    fsom <- TrackSOM(
        inputFiles = data_files,
        dataFileType = '.csv',
        colsToUse = c("x", "y", "z"),
        nClus = 10,
        tracking = TRUE,
        noMerge = TRUE
    )
    dat <- lapply(data_files, function(d) fread(d))
    dat <- rbindlist(dat)

    timepoints <- c('Mock', 'SYN-1', 'SYN-2', 'SYN-3', 'SYN-4')
    dat <- ConcatenateClusteringDetails(tracksom.result = fsom,
                                        dat = dat,
                                        timepoint.col = 'timepoint',
                                        timepoints = timepoints)

    dat_meta <- unique(dat[, c("timepoint", "TrackSOM_metacluster_lineage_tracking")])

    for (t in timepoints) {
        dat_t <- dat_meta[dat_meta$timepoint == t]
        expect_gte(nrow(dat_t), 10)
    }
})

test_that("No merging works for prescribed variant", {

    data_files <- sapply(c(0:4), function(i) {
        system.file("extdata", paste0("synthetic_d", i, ".csv"), package = "TrackSOM")
    })

    nclust <- c(10, 15, 20, 21, 30)
    fsom <- TrackSOM(
        inputFiles = data_files,
        dataFileType = '.csv',
        colsToUse = c("x", "y", "z"),
        nClus = nclust,
        tracking = TRUE,
        noMerge = TRUE
    )
    dat <- lapply(data_files, function(d) fread(d))
    dat <- rbindlist(dat)

    timepoints <- c('Mock', 'SYN-1', 'SYN-2', 'SYN-3', 'SYN-4')
    dat <- ConcatenateClusteringDetails(tracksom.result = fsom,
                                        dat = dat,
                                        timepoint.col = 'timepoint',
                                        timepoints = timepoints)

    dat_meta <- unique(dat[, c("timepoint", "TrackSOM_metacluster_lineage_tracking")])

    for (i in c(1:5)) {
        t <- timepoints[i]
        nclus <- nclust[i]
        dat_t <- dat_meta[dat_meta$timepoint == t]
        expect_gte(nrow(dat_t), nclus)

    }
})

test_that("No merging works for autonomous adaptive", {
    data_files <- sapply(c(0:4), function(i) {
        return(paste0("~/Documents/GitHub/TrackSOM/inst/extdata/synthetic_d", i,".csv"))
    })


    data_files <- sapply(c(0:4), function(i) {
        system.file("extdata", paste0("synthetic_d", i, ".csv"), package = "TrackSOM")
    })
    fsom <- TrackSOM(
        inputFiles = data_files,
        dataFileType = '.csv',
        colsToUse = c("x", "y", "z"),
        maxMeta = 20,
        tracking = TRUE,
        noMerge = TRUE
    )
    dat <- lapply(data_files, function(d) fread(d))
    dat <- rbindlist(dat)

    timepoints <- c('Mock', 'SYN-1', 'SYN-2', 'SYN-3', 'SYN-4')
    dat <- ConcatenateClusteringDetails(tracksom.result = fsom,
                                        dat = dat,
                                        timepoint.col = 'timepoint',
                                        timepoints = timepoints)

    dat_meta <- unique(dat[, c("timepoint", "TrackSOM_metacluster_lineage_tracking")])

    for (t in timepoints) {
        dat_t <- dat_meta[dat_meta$timepoint == t]
        # can't predict how many it will produce
        expect_gt(nrow(dat_t), 2)
    }
})

