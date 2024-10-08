# All TrackSOM functions ----

#' TrackSOM
#'
#' Run the TrackSOM algorithm.
#' Will run Tracking by default with no merging allowed.
#'
#' @param inputFiles NO DEFAULT. Either a vector of filenames or a list of `data.table` objects
#' or a `data.table` object.
#' If a vector, each element should correspond to a file belonging to a timepoint.
#' If a list, each element should correspond to a `data.table` object belonging to a timepoint.
#' If a `data.table` object, make sure `timepoints` and `timepointCol` is not NULL.
#' Make sure the vector or list is ordered such that the 1st element corresponds to
#' the 1st time-point, the 2nd element corresponds to the 2nd time-point, etc.
#' @param colsToUse NO DEFAULT. Vector of column names to be used to train SOM.
#' One list which must be same for all time periods.
#' @param timepoints DEFAULT NULL. A vector of time-points in order.
#' Only used if inputFiles is a `data.table` object
#' @param timepointCol DEFAULT NULL. A character denoting the column in inputFile that
#' corresponds to the time-point each cell belongs to.
#' Only used if inputFiles is a `data.table` object
#' @param maxMeta DEFAULT NULL. Numeric. Maximum number of meta clusters for all
#' time point. TrackSOM use this to limit the number of meta clusters allowed
#' when it is given the freedom to decide the optimum number of meta clusters.
#' @param nClus DEFAULT NULL. Numeric or Vector of numbers. If single number,
#' TrackSOM will produce same number of meta clusters per time point.
#' Otherwise, it will use the entries in the vector. Number of meta clusters
#' must be at least 3 for each time point.
#' @param tracking DEFAULT TRUE. Whether to track cluster changes (TRUE), or
#' just cluster (NULL).
#' @param noMerge DEFAULT FALSE. Whether to allow meta-cluster merging (FALSE) 
#' or not (TRUE).
#' @param xdim DEFAULT 10. SOM grid size.
#' @param ydim DEFAULT 10. SOM grid size.
#' @param seed DEFAULT 42. Random seed number.
#' @param ... other parameters for FlowSOM functions ReadInput, BuildSOM and
#'   BuildMST. See FlowSOM vignette for specific parameter information
#'
#' @return A list containing the following elements:
#' 1) FlowSOM: a list which is very similar to output of FlowSOM function.
#' 2) Metaclustering: a list of two lists:
#'   2.1) Metaclustering: contains the clustering levels for each node in grid
#'   for each time step.
#'   2.2) mcPerCell: contains the clusering levels for each datapoint for each
#'   time step.
#' 3) Tracking: a list of containing the unique IDs as determined by the
#' tracking function for each time step
#'
#' @usage TrackSOM (inputFiles, colsToUse, maxMeta, ...)
#'
#' @importFrom flowCore logicleTransform
#' @import data.table
#' @import FlowSOM
#'
#' @examples
#'
#' data_files <- sapply(c(0:4), function(i) {
#' system.file("extdata", paste0("synthetic_d", i, ".fcs"), package="TrackSOM")
#' })
#' use_cols <- c("x", "y", "z")
#' tracksom_result <- TrackSOM(inputFiles = data_files,
#'                             colsToUse = use_cols,
#'                             nClus = 10,
#'                             dataFileType = ".fcs"
#' )
#' tracksom_res <- ExportClusteringDetailsOnly(tracksom_result)
#'
#' @export

TrackSOM <- function(inputFiles, colsToUse,
                     timepoints = NULL,
                     timepointCol = NULL,
                     tracking = TRUE,
                     noMerge = FALSE,
                     maxMeta = NULL,
                     nClus = NULL,
                     seed = 42,
                     xdim = 10,
                     ydim = 10,
                     rlen = 10,
                     dataFileType = c("data.frame", ".csv", ".fcs"),
                     compensate = FALSE,
                     spillover = NULL,
                     transform = FALSE,
                     toTransform = NULL,
                     transformFunction = flowCore::logicleTransform(),
                     scale = FALSE,
                     scaled.center = TRUE,
                     scaled.scale = TRUE,
                     silent = FALSE,
                     ...) {
  set.seed(seed)

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
  # 2. Check the maxMeta is greater than 3, otherwise you are running the risk
  # of setting number of metaclusters to 2, which will cause
  # ConsensusClusterPlus to fail.
  else {
    if (maxMeta <= 3) {
      stop("Please set maxMeta to be more than 3. Otherwise you are running risk of
           optimal metaclusters to be set to 2, which will cause ConsensusClusterPlus to fail.")
    }
  }

  # Read files


  fsom <- ParseInput(
      inputFiles = inputFiles, 
        compensate = compensate, 
        spillover = spillover,
        transform = transform,
        toTransform = toTransform,
        transformFunction = transformFunction,
        scale = scale,
        scaled.center = scaled.center,
        scaled.scale = scaled.scale,
        silent = silent,
        colsToUse = colsToUse,
        timepoints = timepoints,
        timepointCol = timepointCol
  )



  ## STEP 1: Build SOM on data from all time points

  # Build SOM on all data (all time periods concatenated)
  fsom <-
    BuildSOM(
      fsom,
      colsToUse = colsToUse,
      xdim = xdim,
      ydim = ydim,
      silent = silent,
      ...
    )

  # Build MST on all data, this builds using codes from final grid only (as of
  # now)
  fsom <- BuildMST(fsom, silent = silent)

  # STEP 2: Reconstruct the SOM for individual time point. From the big SOM
  # above, for each time point, go through each node then and see if there is
  # any data in the node belonging to that time point. if yes, collate them, and
  # compute the node's coordinate based on the data therein (for that time point
  # only)

  # Create list for building metaclustering from each time period
  cl <- list()

  # Create empty list for metaclustering per cell in to be added for each time
  # period
  metaPerCell <- list()

  codes <- vector("list", length(inputFiles))
  cluster_codes <- vector("list", length(inputFiles))

  if (!silent) {
      message("Extracting SOM nodes for each time point")
  }
  
  # Loop through data subsets and create new codes (a matrix of averages for
  # each column) based on data in a time point only.
  for (i in 1:length(inputFiles)) {
    # obtain the data just for a single time point (whichever is pointed by
    # variable i)
    subset_ <-
      FlowSOMSubset(fsom, (fsom$metaData[[i]][1]):(fsom$metaData[[i]][2]))

    # 1st apply, will store the new centroid of each SOM node in a matrix. the
    # 2nd apply: get the centroid of each SOM's node based on the mean of only
    # data points for this time point. apply(..,2,mean): ... subset data so only
    # include data for a SOM node. 2 means compute mean over columns, mean is
    # well mean.
    code <- as.matrix(t(sapply(seq_len(fsom$map$nNodes), function(x) {
      apply(subset(subset_$data, subset_$map$mapping[,1] == x),2,mean)
    })))

    # only keep the columns to be used for clustering
    code <- code[, colsToUse]
    code <- ifelse(is.nan(code), NA, code)
    codes[[i]] <- code

  }
  # Loop through files and create a metaclustering for each file, based on codes
  # (mean of each datapoint in given node at time point)
  fsom$map$coding <- codes

  # STEP 3: for each time point, run meta clustering on the reconstructed SOM
  # TODO: make it flexible such that user can define number of meta clusters
  # per time point rather than leaving it to TrackSOM to work out.

  # Run clustering - this will perform clustering only on nodes with datapoints
  # in a given time period

  if (!silent) {
      message("Running meta clustering")
  }
  
  for (i in 1:(length(inputFiles))) {
    map <- fsom$map$coding[[i]]
    map2 <- map[,1]
    # only cluster on nodes with datapoints in the given time period
    cluster_codes <- map[rowSums(!is.na(map)) > 0,]

    num_codes <- nrow(cluster_codes)
    
    if (!silent) {
        message(paste("Meta clustering time point", i, "with", num_codes,
                      "SOM nodes"))
    }
    
    
    if (is.null(nClus)){
      # if the requested max meta is larger than number of nodes available, the
      # following will fail. Thus need to just exist here and fail
      if (maxMeta > num_codes) {
        stop(
          paste(
            "Time point",
            i,
            "only contains",
            num_codes,
            "SOM nodes, less than requested number of meta clusters",
            maxMeta,
            '. Please reduce it.'
          )
        )
      }

      # run meta-clustering
      mc <-
        MetaClustering(cluster_codes,
                       method = "metaClustering_consensus",
                       max = maxMeta,
                       seed = seed)
      mc <- as.factor(clustWithNAs(map2, mc))

    } else {
      # check if list or single value for number of clusters
      if (length(nClus) == 1) {
        if (nClus > num_codes) {
          stop(
            paste(
              "Time point",
              i,
              "only contains",
              num_codes,
              "SOM nodes, less than requested number of meta clusters",
              nClus,
              '. Please reduce it.'
            )
          )
        }
        mc <- metaClustering_consensus(cluster_codes, k = nClus, seed = seed)
        mc <- as.factor(clustWithNAs(map2, mc))

      } else {
        if (nClus[[i]] > num_codes) {
          stop(
            paste(
              "Time point",
              i,
              "only contains",
              num_codes,
              "SOM nodes, less than requested number of meta clusters",
              nClus[[i]],
              '. Please reduce it.'
            )
          )
        }
        mc <- metaClustering_consensus(cluster_codes, k = nClus[[i]], seed = seed)
        mc <- as.factor(clustWithNAs(map2, mc))
      }

    }

    cl <- c(cl, setNames(list(mc), paste("metacluster", i, sep="_")))

    # create named list of metacluster for each cell
    mcl <-
      (mc[fsom$map$mapping[(fsom$metaData[[i]][1]):(fsom$metaData[[i]][2]), 1]])
    metaPerCell <-
      c(metaPerCell, setNames(list(mcl), paste("metacluster", i, sep = "_")))
  }


  meta <- list()
  meta$metaclustering <- cl
  meta$mcPerCell <- metaPerCell
  fsom <- list("FlowSOM"=fsom, "metaclustering"=meta)

  if (isTRUE(tracking)) {
    if (isTRUE(noMerge)) {
      fsom$tracking$lineage <- TrackingNoMerging(fsom)
    } else {
      fsom$tracking$lineage <- TrackingWithMerging(fsom)
    }
  }

  fsom
}


#' clustWithNAs
#'
#' Function used to tease apart with the meta-clustering result for the SOM nodes.
#'
#'
#' @param values a list of SOM nodes
#' @param mcs a list of meta-clusters
#'
#' @return an updated list of SOM nodes
#'
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


