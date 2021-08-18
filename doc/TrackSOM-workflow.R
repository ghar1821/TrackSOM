## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(TrackSOM)

## -----------------------------------------------------------------------------
data.files.fullpath <- c(
  system.file("extdata", "synthetic_d0.csv", package = "TrackSOM"),
  system.file("extdata", "synthetic_d1.csv", package = "TrackSOM"),
  system.file("extdata", "synthetic_d2.csv", package = "TrackSOM"),
  system.file("extdata", "synthetic_d3.csv", package = "TrackSOM"),
  system.file("extdata", "synthetic_d4.csv", package = "TrackSOM")
)

## -----------------------------------------------------------------------------
print(data.files.fullpath)

## -----------------------------------------------------------------------------
data.files.fullpath.fcs <- c(
  system.file("extdata", "synthetic_d0.fcs", package = "TrackSOM"),
  system.file("extdata", "synthetic_d1.fcs", package = "TrackSOM"),
  system.file("extdata", "synthetic_d2.fcs", package = "TrackSOM"),
  system.file("extdata", "synthetic_d3.fcs", package = "TrackSOM"),
  system.file("extdata", "synthetic_d4.fcs", package = "TrackSOM")
)

## -----------------------------------------------------------------------------
print(data.files.fullpath.fcs)

## -----------------------------------------------------------------------------

library(data.table)
data.files <- c(
  system.file("extdata", "synthetic_d0.csv", package = "TrackSOM"),
  system.file("extdata", "synthetic_d1.csv", package = "TrackSOM"),
  system.file("extdata", "synthetic_d2.csv", package = "TrackSOM"),
  system.file("extdata", "synthetic_d3.csv", package = "TrackSOM"),
  system.file("extdata", "synthetic_d4.csv", package = "TrackSOM")
)
dat <- lapply(data.files, function(f) fread(f))

## -----------------------------------------------------------------------------
dat

## -----------------------------------------------------------------------------
tracksom.result <- TrackSOM(inputFiles = data.files.fullpath,
                            colsToUse = c('x', 'y', 'z'),
                            tracking = TRUE,
                            noMerge = TRUE,
                            nClus = c(3,3,9,7,15),
                            dataFileType = ".csv"
)

## -----------------------------------------------------------------------------
library(data.table)
InputDirectory <- "~/Documents/GitHub/TrackSOM/inst/extdata/"
data.files <- list.files(InputDirectory, ".csv")
dat <- lapply(data.files, function(f) fread(paste0(InputDirectory, f)))
dat <- rbindlist(dat)

## -----------------------------------------------------------------------------
head(dat)

## -----------------------------------------------------------------------------
tail(dat)

## -----------------------------------------------------------------------------
dat.clust <- ConcatenateClusteringDetails(
    tracksom.result = tracksom.result,
    dat = dat,
    timepoint.col = "timepoint",
    timepoints = c('Mock', 'SYN-1', 'SYN-2', 'SYN-3', 'SYN-4')
    )

## -----------------------------------------------------------------------------
head(dat.clust)

## -----------------------------------------------------------------------------
DrawNetworkPlot(dat = dat.clust,
                timepoint.col = "timepoint",
                timepoints = c('Mock', 'SYN-1', 'SYN-2', 'SYN-3', 'SYN-4'),
                cluster.col = 'TrackSOM_metacluster_lineage_tracking',
                marker.cols = c('x', 'y', 'z'))

## -----------------------------------------------------------------------------
list.files()

## -----------------------------------------------------------------------------
DrawTimeseriesHeatmap(dat = dat.clust,
                timepoint.col = "timepoint",
                timepoints = c('Mock', 'SYN-1', 'SYN-2', 'SYN-3', 'SYN-4'),
                cluster.col = 'TrackSOM_metacluster_lineage_tracking',
                marker.cols = c('x', 'y', 'z'))

## -----------------------------------------------------------------------------
list.files()

