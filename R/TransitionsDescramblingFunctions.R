# Helper functions to de-scramble meta-cluster ID ----

#' Get round bracket
#'
#' Helper method to find round brackets in meta-cluster ID
#'
#' @param cl Character cluster ID
#'
get_idx_roundclsbracket <- function(cl) {
  tail(gregexpr("\\)", cl)[[1]], n=1)
}

#' Get pipe character
#'
#' Helper method to detect pipe character in meta-cluster ID
#'
#' @param cl Character cluster ID
#'
#'
get_idx_pipe <- function(cl) {
  tail(gregexpr("\\|", cl)[[1]], n=1)
}

#' GetTransitionsAsEdges
#'
#' Given a TrackSOM clustering and tracking, descramble the meta-cluster ID.
#' Handy to work out the meta-cluster transitions.
#'
#' @param dat Data.table storing the clustered gated dataset
#' @param timepoints Vector containing the time points in order
#' @param timepoint_col Character Column denoting time points in dat
#' @param cluster_col Character Column denoting cluster id produced by tracking
#'   by lineage determination in dat
#'
#' @import gtools
#' @import stringr
#' @import data.table
#'
#' @examples
#'
#' library(data.table)
#' data_files <- sapply(c(0:4), function(i) {
#'   system.file("extdata", paste0("synthetic_d", i, ".fcs"), package="TrackSOM")
#' })
#' use_cols <- c("x", "y", "z")
#' tracksom_result <- TrackSOM(inputFiles = data_files,
#'                             colsToUse = use_cols,
#'                             nClus = 10,
#'                             dataFileType = ".fcs"
#' )
#'
#' data <- lapply(c(0:4), function(i) {
#'   fread(system.file("extdata", paste0("synthetic_d", i, ".csv"), package="TrackSOM"))
#' })
#' data <- rbindlist(data)
#' tracksom_res <- ConcatenateClusteringDetails(
#'   dat = data,
#'   tracksom.result = tracksom_result,
#'   timepoint.col = 'timepoint',
#'   timepoints = c("Mock", "SYN-1", "SYN-2", "SYN-3", "SYN-4")
#' )
#' GetTransitionsAsEdges(tracksom_res,
#'                       'timepoint',
#'                       c("Mock", "SYN-1", "SYN-2", "SYN-3", "SYN-4"),
#'                       "TrackSOM_metacluster_lineage_tracking")
#'
#' @export
#'
GetTransitionsAsEdges <- function(dat,
                                  timepoints,
                                  timepoint_col,
                                  cluster_col) {

  # Remove cells not assigned clusters
  dat_bk <- dat[!is.na(dat[[cluster_col]]),]

  edge_df <- data.frame(from=character(),
                        to=character())

  #### Find all the transitions ####
  for (tp_idx in c(2:length(timepoints))) {
    prev_tp_dat <- dat_bk[dat_bk[[timepoint_col]] == timepoints[tp_idx-1], ]
    prev_tp_clust <- mixedsort(unique(prev_tp_dat[[cluster_col]]))

    complex_clusters <-
      as.vector(prev_tp_clust[lapply(prev_tp_clust, get_idx_roundclsbracket) > -1])
    simple_n_pipe_clusters <- setdiff(prev_tp_clust, complex_clusters)
    pipe_clusters <-
      simple_n_pipe_clusters[sapply(simple_n_pipe_clusters, get_idx_pipe) > -1]
    simple_clusters <- setdiff(simple_n_pipe_clusters, pipe_clusters)

    # any of the above could have empty character in it which we don't want.
    # this can happen when no merging is allowed. Complex_clusters is empty
    complex_clusters <- complex_clusters[!complex_clusters %in% ""]
    pipe_clusters <- pipe_clusters[!pipe_clusters %in% ""]
    simple_clusters <- simple_clusters[!simple_clusters %in% ""]

    # Order the vector based on the length of each element
    simple_clusters <-
      simple_clusters[order(nchar(simple_clusters), simple_clusters, decreasing = TRUE)]
    pipe_clusters <-
      pipe_clusters[order(nchar(pipe_clusters), pipe_clusters, decreasing = TRUE)]
    complex_clusters <-
      complex_clusters[order(nchar(complex_clusters), complex_clusters, decreasing = TRUE)]

    curr_tp_dat <- dat_bk[dat_bk[[timepoint_col]] == timepoints[tp_idx], ]
    curr_tp_clust <- mixedsort(unique(curr_tp_dat[[cluster_col]]))

    # this filter out clusters that only exist in current time point
    # these are new clusters and have no predecessors
    curr_tp_clust_existing <-
      sapply(as.vector(curr_tp_clust), function(cl) {
        if (grepl(",", cl, fixed = TRUE) ||
            grepl("|", cl, fixed = TRUE) || cl %in% simple_clusters) {
          return(cl)
        }
      })
    curr_tp_clust_existing <- unlist(curr_tp_clust_existing)

    for (cl in curr_tp_clust_existing) {
      clean_cl <- cl
      ## Match the complex clusters first
      for (cls in complex_clusters) {
        cls_found <- grepl(cls, clean_cl, fixed = TRUE)
        if (cls_found) {
          df <- data.frame(paste0(tp_idx-1,'_', cls), paste0(tp_idx, '_', cl))
          names(df) <- c("from","to")
          edge_df <- rbind(edge_df, df)
          clean_cl <- gsub(cls, "", clean_cl, fixed = TRUE)
        }
      }

      ## Match the split simple clusters first
      for (cls in pipe_clusters) {
        cls_found <- grepl(cls, clean_cl, fixed = TRUE)
        if (cls_found) {

          # To prevent matching something like RR|1 to R|1, we replace R|1 from
          # RR|1 and if there is alphabet, then it's not the right match
          clean_cl_without_cls <- gsub(cls, "", clean_cl, fixed = TRUE)
          character_in_clean_cl_without_cls <- grepl("^[A-Za-z]+$", clean_cl_without_cls, perl = T)

          if (!character_in_clean_cl_without_cls) {
            df <- data.frame(paste0(tp_idx-1,'_', cls), paste0(tp_idx, '_', cl))
            names(df) <- c("from","to")
            edge_df <- rbind(edge_df, df)
            clean_cl <- gsub(cls, "", clean_cl, fixed = TRUE)
          }

          df <- data.frame(paste0(tp_idx-1,'_', cls), paste0(tp_idx, '_', cl))
          names(df) <- c("from","to")
          edge_df <- rbind(edge_df, df)
          clean_cl <- gsub(cls, "", clean_cl, fixed = TRUE)
        }
      }

      ## Then simple clusters
      for (cls in simple_clusters) {
        cls_found <- grepl(cls, clean_cl, fixed = TRUE)
        if (cls_found) {
          df <- data.frame(paste0(tp_idx-1,'_',cls), paste0(tp_idx,'_',cl))
          names(df) <- c("from","to")
          edge_df <- rbind(edge_df, df)
          clean_cl <- gsub(cls, "", clean_cl, fixed = TRUE)
        }
      }
    }
  }
  return(edge_df)
}
