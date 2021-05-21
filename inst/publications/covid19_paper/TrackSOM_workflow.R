# Sample TrackSOM workflow

# Import libraries
library(pdist) # for distances
library(rlist) # Needed fot GetTransitions function
library(Spectre)

code.dir <- '~/Documents/phd/code/TrackSOM/TrackSOM'
source(paste0(code.dir, '/TrackSOM.R'))
source(paste0(code.dir, '/TrackingFunctions.R'))
source(paste0(code.dir, '/UtilityFunctions.R'))

setwd("~/Documents/phd/code/TrackSOM/covid19_paper/binned.data")
PrimaryDirectory <- getwd()


data.files <- list.files(PrimaryDirectory, ".fcs")
# convert to absolute path
data.files.fullpath <- sapply(data.files, function(fname) {
  full.fname <- paste(PrimaryDirectory, fname, sep="/")
})

# column to be used by TrackSOM to do clustering
ClusteringCols <- c("CD19_asinh_noiseRed_aligned",
                    "HLA-DR_asinh_noiseRed_aligned",
                    "CD8_asinh_noiseRed_aligned",
                    "CD16_asinh_noiseRed_aligned",
                    "CD45RA_asinh_noiseRed_aligned",
                    "TCRgd_asinh_noiseRed_aligned",
                    "CD3_asinh_noiseRed_aligned",
                    "CD56_asinh_noiseRed_aligned",
                    "CD4_asinh_noiseRed_aligned",
                    "CD27_asinh_noiseRed_aligned",
                    "CD14_asinh_noiseRed_aligned")

tracksom.result <- TrackSOM(inputFiles = data.files.fullpath,
                            colsToUse = ClusteringCols,
                            tracking = TRUE,
                            seed = 42,
                            xdim = 10,
                            ydim = 10,
                            nClus = 40
)

dat.dir <- '~/Documents/phd/code/TrackSOM/covid19_paper/binned.data.csv'
dat.list <- Spectre::read.files(file.loc = dat.dir,
                                file.type = ".csv")
cell.dat <- data.table::rbindlist(dat.list, fill = TRUE)
head(cell.dat)

cell.dat <- consolidate.tracksom.result(tracksom.result = tracksom.result,
                                        dat = cell.dat,
                                        divide.by = "FileName")
setwd("~/Documents/phd/code/TrackSOM/covid19_paper")
dir.create("TrackSOM_result")
setwd("TrackSOM_result")
Spectre::write.files(cell.dat, "Tracked", divide.by = "FileName")
