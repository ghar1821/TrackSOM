# Sample TrackSOM workflow

# Import libraries
library(pdist) # for distances
library(rlist) # Needed fot GetTransitions function
library(Spectre)

source('~/Documents/phd/code/FlowSOM-tracking/TrackSOM/TrackSOM.R')
source('~/Documents/phd/code/FlowSOM-tracking/TrackSOM/TrackingFunctions.R')
source('~/Documents/phd/code/FlowSOM-tracking/TrackSOM/UtilityFunctions.R')
source('~/Documents/phd/code/FlowSOM-tracking/TrackSOM/4_metaClustering.R')

setwd("/Users/givanna/Documents/phd/tracksom")
PrimaryDirectory <- getwd()

InputDirectory <- "~/Documents/phd/tracksom/disappearing_cluster/data_no_label/fcs_files"

data.files <- list.files(InputDirectory, ".fcs")
# convert to absolute path
data.files.fullpath <- sapply(data.files, function(fname) {
  full.fname <- paste(InputDirectory, fname, sep="/")
})

# column to be used by TrackSOM to do clustering
ClusteringCols <- c("x", "y", "z")


## Option 1: FlowSOM infer the optimal number of meta clusters. 
tracksom.result <- TrackSOM(inputFiles = data.files.fullpath,
                            colsToUse = ClusteringCols,
                            tracking = TRUE,
                            seed = 42,
                            xdim = 10,
                            ydim = 10,
                            rlen=8,
                            alpha=c(0.0796848513875157, 0.0264928490658524),
                            maxMeta = 4
)

## Option 2: FlowSOM creates same number of meta clusters per time point.
tracksom.result <- TrackSOM(inputFiles = data.files.fullpath,
                            colsToUse = ClusteringCols,
                            tracking = TRUE,
                            seed = 42,
                            xdim = 10,
                            ydim = 10,
                            rlen=8,
                            alpha=c(0.0796848513875157, 0.0264928490658524),
                            nClus = 3
)

## Option 3: FlowSOM creates different number of meta clusters per time point.
meta.per.timepoint <- as.list(c(3, 3, 9, 7, 15))
tracksom.result <- TrackSOM(inputFiles = data.files.fullpath,
                            colsToUse = ClusteringCols,
                            tracking = TRUE,
                            seed = 42,
                            xdim = 8,
                            ydim = 8,
                            rlen=8,
                            alpha=c(0.0796848513875157, 0.0264928490658524),
                            nClus = meta.per.timepoint
)

cell.dat <- Spectre::read.files(file.loc = InputDirectory,
                                file.type = ".fcs")
cell.dat <- Spectre::do.merge.files(cell.dat)
head(cell.dat)

cell.dat <- consolidate.tracksom.result(tracksom.result = tracksom.result,
                                        dat = cell.dat,
                                        divide.by = "FileName")

setwd("~/Documents/phd/tracksom")
Spectre::write.files(cell.dat, "Result", divide.by = "FileName")
