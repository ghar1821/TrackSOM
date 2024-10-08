## ----message=FALSE, warning=FALSE---------------------------------------------
library(data.table)
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
dat <- lapply(data.files.fullpath, function(f) fread(f))
dat

## -----------------------------------------------------------------------------
timepoints <- seq(0, 4)

dat <- lapply(seq(length(data.files.fullpath)), function(data_file_i) {
    dt <- fread(data.files.fullpath[[data_file_i]])
    dt[['timepoint']] <- timepoints[data_file_i]
    return(dt)
})

dat <- rbindlist(dat)

head(dat)
tail(dat)

## -----------------------------------------------------------------------------
tracksom.result <- TrackSOM(inputFiles = data.files.fullpath,
                            colsToUse = c('x', 'y', 'z'),
                            tracking = TRUE,
                            noMerge = TRUE,
                            nClus = c(3,3,9,7,15),
                            dataFileType = ".csv"
)

## -----------------------------------------------------------------------------
data.files <- c(
  system.file("extdata", "synthetic_d0.csv", package = "TrackSOM"),
  system.file("extdata", "synthetic_d1.csv", package = "TrackSOM"),
  system.file("extdata", "synthetic_d2.csv", package = "TrackSOM"),
  system.file("extdata", "synthetic_d3.csv", package = "TrackSOM"),
  system.file("extdata", "synthetic_d4.csv", package = "TrackSOM")
)
dat <- lapply(data.files, function(f) fread(f))
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

