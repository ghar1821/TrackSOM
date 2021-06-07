# Sample TrackSOM workflow

# Import libraries
library(data.table)
library(TrackSOM)

# 1) Specify data files ----
## If using CSV file ----
InputDirectory <- "~/Documents/GitHub/TrackSOM/inst/extdata/"
data.files <- list.files(InputDirectory, ".csv")
# convert to absolute path
data.files.fullpath <- sapply(data.files, function(fname) {
    full.fname <- paste(InputDirectory, fname, sep="/")
})

## If using FCS file ----
InputDirectory <- "~/Documents/GitHub/TrackSOM/inst/extdata/"
data.files <- list.files(InputDirectory, ".fcs")
# convert to absolute path
data.files.fullpath <- sapply(data.files, function(fname) {
    full.fname <- paste(InputDirectory, fname, sep="/")
})

## If using data.frame or data.table ----
InputDirectory <- "~/Documents/GitHub/TrackSOM/inst/extdata/"
data.files <- list.files(InputDirectory, ".csv")
dat <- lapply(data.files, function(f) fread(paste0(InputDirectory, f)))

# 2) Specify more details and run TrackSOM ----
# column to be used by TrackSOM to do clustering
ClusteringCols <- c("x", "y", "z")

## Option 1: FlowSOM infer the optimal number of meta clusters. 
tracksom.result <- TrackSOM(inputFiles = data.files.fullpath,
                            colsToUse = ClusteringCols,
                            tracking = TRUE,
                            noMerge = TRUE,  # TODO change me if you want to allow merging
                            seed = 42,
                            xdim = 10,
                            ydim = 10,
                            maxMeta = 4,
                            dataFileType = ".fcs"  # TODO change me according to file type you have
)

## Option 2: FlowSOM creates same number of meta clusters per time point.
tracksom.result <- TrackSOM(inputFiles = data.files.fullpath,
                            colsToUse = ClusteringCols,
                            tracking = TRUE,
                            noMerge = TRUE,  # TODO change me if you want to allow merging
                            seed = 42,
                            xdim = 10,
                            ydim = 10,
                            nClus = 10,
                            dataFileType = ".csv"  # TODO change me according to file type you have
)

## Option 3: FlowSOM creates different number of meta clusters per time point.
meta.per.timepoint <- as.list(c(3, 3, 9, 7, 15))
tracksom.result <- TrackSOM(inputFiles = dat,
                            colsToUse = ClusteringCols,
                            tracking = TRUE,
                            noMerge = TRUE,  # TODO change me if you want to allow merging
                            seed = 42,
                            xdim = 8,
                            ydim = 8,
                            nClus = meta.per.timepoint,
                            dataFileType = "data.frame"  # TODO change me according to file type you have
)

# 3) Append the result to original data ----
## If you get TrackSOM to read the csv or fcs file, read those files in ----
cell.dat <- Spectre::read.files(file.loc = InputDirectory, file.type = '.csv')
cell.dat <- Spectre::do.merge.files(cell.dat)

## If using data.table, then just concatenate it into 1 giant data.table
cell.dat <- rbindlist(dat)

cell.dat <-
    ConcatenateClusteringDetails(
        tracksom.result = tracksom.result,
        dat = cell.dat,
        timepoint.col = "timepoint",
        timepoints = c('Mock', 'SYN-1', 'SYN-2', 'SYN-3', 'SYN-4')
    )
outdir <- "~/Documents/phd/tracksom/sample_result"
dir.create(outdir, recursive = TRUE)
setwd(outdir)
Spectre::write.files(cell.dat, "Result", divide.by = "timepoint")

# make some visualisation
TrackSOM::DrawNetworkPlot(dat = cell.dat,
                          timepoint.col = "timepoint",
                          timepoints = c('Mock', 'SYN-1', 'SYN-2', 'SYN-3', 'SYN-4'),
                          cluster.col = 'TrackSOM_metacluster_lineage_tracking',
                          marker.cols = ClusteringCols)
TrackSOM::DrawTimeseriesHeatmap(dat = cell.dat,
                                timepoint.col = "timepoint",
                                timepoints = c('Mock', 'SYN-1', 'SYN-2', 'SYN-3', 'SYN-4'),
                                cluster.col = 'TrackSOM_metacluster_lineage_tracking',
                                marker.cols = ClusteringCols)
