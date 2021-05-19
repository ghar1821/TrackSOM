# All FlowSOM Series functions here

# Import libraries
library(FlowSOM)
library(flowCore)

#' Run TrackSOM algorithm
#' TrackSOM builds on FlowSOM to accomodate the tracking of cluster evolution using time-series input data
#'
#' OUTPUTS: a list of 2-3 lists (3 if tracking is set to TRUE)
#' 1) FlowSOM - a list which is very similar to output of FlowSOM function. Adds "coding" for each time step under map
#' 2) Metaclustering - a list of two lists (Metaclustering and mcPerCell)
#'     - Metaclustering: contains the clustering levels for each node in grid for each time step
#'     - mcPerCell: contains the clusering levels for each datapoint for each time step
#' 3) Tracking (optional) - a list of two lists (lineage and proximity)
#'     - lineage: contains the unique IDs as determined by lineage determination for each time step
#'     - proximity: contains the unique IDs as determined by historical proximity for each time step
#' 
#' @usage TrackSOM (inputFiles, colsToUse, maxMeta, ...)
#' 
#' @param inputFiles NO DEFAULT. Vector of filenames, each corresponds to the file containing data for a time point
#' @param colsToUse NO DEFAULT. Vector of column names to be used to train SOM. One list which must be same for all time periods.
#' @param maxMeta DEFAULT NULL. Numeric. Maximum number of meta clusters for all time point. TrackSOM use this to limit the number of meta clusters allowed when it is given the freedom to decide the optimum number of meta clusters.
#' @param nClus DEFAULT NULL. Numeric or Vector of numbers. If single number, TrackSOM will produce same number of meta clusters per time point. Otherwise, it will use the entries in the vector.
#' Number of meta clusters must be at least 3 for each time point.
#' @param tracking DEFAULT NULL. Whether to track cluster changes (TRUE), or just cluster (NULL)
#' @param ... other parameters for FlowSOM functions ReadInput, BuildSOM and BuildMST. See FlowSOM vignette for specific parameter information
#'
#'
#' @export

TrackSOM <- function(inputFiles, colsToUse, 
                     tracking=TRUE, 
                     maxMeta=NULL, 
                     nClus=NULL,
                     seed=NULL, 
                     xdim=10, 
                     ydim=10, 
                     rlen=10,
                     pattern = '.fcs', 
                     compensate=FALSE, 
                     spillover=NULL,
                     transform=FALSE, 
                     toTransform=NULL, 
                     transformFunction=flowCore::logicleTransform(),
                     scale=TRUE,
                     scaled.center=TRUE, 
                     scaled.scale=TRUE, 
                     silent=FALSE, 
                     tSNE=FALSE,
                     importance=NULL, 
                     ...) {
  
  ### Some checks to make sure that FlowSOM can run properly ###
  # 1. Meta clusters are greater than grid size and that it's greater than 2
  if (is.null(maxMeta)) {
    grid.size <- xdim * ydim
    for (num.clust in c(nClus)) {
      if (grid.size < num.clust) {
        stop("Cannot have grid size ", 
             xdim, " x ", ydim,
             " (xdim x ydim)  smaller than number of metaclusters ",
             nClus)
      }
      if (num.clust <= 2) {
        stop("Please set number of meta clusters to be more than 2. Otherwise
             ConsensusClusteringPlus will fail.")
      }
    } 
  }
  # 2. Check the maxMeta is greater than 3, otherwise you are running the risk of 
  # setting number of metaclusters to 2, which will cause ConsensusClusterPlus to fail.
  else {
    if (maxMeta <= 3) {
      stop("Please set maxMeta to be more than 3. Otherwise you are running risk of
           optimal metaclusters to be set to 2, which will cause ConsensusClusterPlus to fail.")
    }
  }
  
  
  
  # files <- c("input_1.fcs", "input_2.fcs", ... , "input_n.fcs")
  fsom <- ReadInput(inputFiles, pattern=pattern, 
                    compensate=compensate, 
                    spillover=spillover, 
                    transform=transform, 
                    toTransform=toTransform, 
                    transformFunction = transformFunction, 
                    scale=scale,
                    scaled.center=scaled.center, 
                    scaled.scale=scaled.scale, 
                    silent=TRUE)
  
  ## STEP 1: Build SOM on data from all time points
  
  # Build SOM on all data (all time periods concatenated)
  fsom <- BuildSOM(fsom, 
                   colsToUse=colsToUse, 
                   xdim=xdim, 
                   ydim=ydim, 
                   rlen=rlen,
                   silent=silent,
                   importance=importance,
                   ...)

  # Build MST on all data, this builds using codes from final grid only (as of now)
  fsom <- BuildMST(fsom, 
                   silent=silent,
                   tSNE=tSNE)

  ## STEP 2: Reconstruct the SOM for individual time point. From the big SOM above, for each time point, go through each node
  ## then and see if there is any data in the node belonging to that time point.
  ## if yes, collate them, and compute the node's coordinate based on the data therein (for that time point only)
  
  # Create list for building metaclustering from each time period
  cl <- list()
  
  # Create empty list for metaclustering per cell in to be added for each time period
  metaPerCell <- list()
  
  codes <- vector("list", length(inputFiles))
  cluster_codes <- vector("list", length(inputFiles))
  
  # Loop through data subsets and create new codes (a matrix of averages for each column) based on data in a time point only
  for (i in 1:length(inputFiles)) {
    # obtain the data just for a single time point (whichever is pointed by variable i)
    subset_ <- FlowSOMSubset(fsom, (fsom$metaData[[i]][1]):(fsom$metaData[[i]][2]))
    
    # 1st apply, will store the new centroid of each SOM node in a matrix.
    # the 2nd apply: get the centroid of each SOM's node based on the mean of only data points for this time point.
    # apply(..,2,mean): ... subset data so only include data for a SOM node. 2 means compute mean over columns, mean is well mean.
    code <- as.matrix(t(sapply(seq_len(fsom$map$nNodes), function(x) {
      apply(subset(subset_$data, subset_$map$mapping[,1] == x),2,mean)
    })))
    
    # only keep the columns to be used for clustering
    code <- code[, colsToUse]
    code <- ifelse(is.nan(code), NA, code)
    codes[[i]] <- code
    
  }
  # Loop through files and create a metaclustering for each file, based on codes (mean of each datapoint in given node at time point)
  fsom$map$coding <- codes
  
  ## STEP 3: for each time point, run meta clustering on the reconstructed SOM
  ## TODO: make it flexible such that user can define number of meta clusters per time point rather than leaving it to TrackSOM to work out.
  
  # Run clustering - this will perform clustering only on nodes with datapoints in a given time period
  for (i in 1:(length(inputFiles))) {
    map <- fsom$map$coding[[i]]
    map2 <- map[,1]
    cluster_codes <- map[rowSums(!is.na(map)) > 0,] # only cluster on nodes with datapoints in the given time period
    
    if (is.null(nClus)){
      mc <- MetaClustering(cluster_codes,method="metaClustering_consensus", max=maxMeta, seed=seed) # run clustering
      mc <- as.factor(clustWithNAs(map2, mc))
      
    } else {
      if (length(nClus) == 1) { # check if list or single value for number of clusters
        mc <- metaClustering_consensus(cluster_codes, k = nClus, seed = seed)
        mc <- as.factor(clustWithNAs(map2, mc))
        
      } else {
        mc <- metaClustering_consensus(cluster_codes, k = nClus[[i]], seed = seed)
        mc <- as.factor(clustWithNAs(map2, mc))
      }
      
    }
    
    cl <- c(cl, setNames(list(mc), paste("metacluster", i, sep="_")))
    
    # create named list of metacluster for each cell
    mcl <- (mc[fsom$map$mapping[(fsom$metaData[[i]][1]):(fsom$metaData[[i]][2]),1]])
    metaPerCell <- c(metaPerCell, setNames(list(mcl), paste("metacluster", i, sep="_")))
  }
  
  
  meta <- list()
  meta$metaclustering <- cl
  meta$mcPerCell <- metaPerCell
  fsom <- list("FlowSOM"=fsom, "metaclustering"=meta)
  
  if (isTRUE(tracking)) {
    fsom$tracking$lineage <- TrackingByLineage(fsom)
    fsom$tracking$proximity <- TrackingByHistoricalProximity(fsom, fsom$tracking$lineage)
  }
  
  fsom
}


# EXAMPLE:
# fsom <- TrackSOM(data_files, 
#                       colsToUse = c(1:14), nClus = 20, xdim=10, ydim=10, rlen=10, tracking=TRUE, silent=FALSE)'

# =========================================================================================

# Function used to deal with clustering since run on only non-missing data

clustWithNAs <- function(values, mcs) {
  j <- 1
  for (i in 1:length(values)) {
    if (!is.na(values[i])) {
      values[i] <- mcs[j]
      j <- j + 1
    }
  }
  return(values)
  
}

# =========================================================================================
