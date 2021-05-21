# All the tracking functions ----

#' TrackingNoMerging
#'
#' Perform tracking of meta-clusters where merging is NOT allowed.
#' Need to be run after FlowSOM's BuildSOM which contorts the SOM grid and
#' MetaClustering which perform meta-cluster of the SOM node.
#'
#' @param fsom list containing a FlowSOM object and meta-clustering lists.
#' @return A list containing N lists of unique IDs as determined by lineage,
#'   where N is no. of time steps.
#'
#' @import pdist
#' @import rlist
#'
#' @seealso \code{\link{TrackSOM}}.
#'
TrackingNoMerging <- function(fsom) {
  # Make lists that are needed for all time periods (clustering,
  # iterations and number)
  new_mcs <- list()
  id_params <- list("n" = 1, "iter" = 1)
  id_split_list <- list()

  # There is no tracking needed for first time step, but we need to get IDs
  # from this time step
  get_ids_1 <- GetFirstID(fsom$metaclustering$metaclustering[[1]],
                          fsom$metaclustering$mcPerCell[[1]])
  new_mcs[[1]] <- get_ids_1$new_ids
  id_params$n <- get_ids_1$n
  id_params$iter <- get_ids_1$iter


  # Loop through remaining time steps and perform tracking
  for (timestep in 2:length(fsom$FlowSOM$metaData)) {
    # For the second time point, we want to base the IDs off the
    # first time point
    mc_prev <- new_mcs[[timestep - 1]]

    # flowsom meta cluster per day - the numeric meta cluster
    # a vector, grid size length. An element is the meta cluster for the node.
    mc <- fsom$metaclustering$metaclustering[[timestep]]
    # meta cluster for each data point.
    mcPerCell <-
      (mc[fsom$FlowSOM$map$mapping[
        (fsom$FlowSOM$metaData[[timestep]][1]):
          (fsom$FlowSOM$metaData[[timestep]][2]),
        1]])

    # List of nodes per metacluster for this time step
    # The name of the list is the metacluster. The content is the node
    # for that metacluster
    nodes_per_clust <- GetNodesPerCluster(mc)
    # List of nodes per metacluster from previous time step
    # same as above, but the name of the list is the tracksom meta cluster id
    # (A, A|1, etc.).
    nodes_per_clust_prev <- GetNodesPerCluster(mc_prev)
    new_ids <-
      vector(mode = "list", length = length(nodes_per_clust)) # Empty list to store IDs
    mc_new <- mc

    # Creates data matrix with number of overlapping nodes between metaclusters
    # from the two time steps.
    # Each row is metacluster for current time point
    # Each column is metacluster for previous time point
    # Each cell contains number of nodes common between meta cluster in current
    # and previous time point
    data <-
      matrix(
        data = NA,
        nrow = length(nodes_per_clust),
        ncol = length(nodes_per_clust_prev)
      )
    for (m in 1:length(nodes_per_clust)) {
      for (n in 1:length(nodes_per_clust_prev)) {
        c <-
          length(which(nodes_per_clust_prev[[n]] %in% nodes_per_clust[[m]]))
        if (c > 0) {
          data[m, n] = c
        }
      }
    }

    # Create new matrix to hold the ID labels for easy merging and splitting of
    # ID labels later
    new_data <-
      matrix(
        data = NA,
        nrow = length(nodes_per_clust),
        ncol = length(nodes_per_clust_prev)
      )

    # Find all rows and columns which is not NA
    # Meaning getting the meta clusters which have shared nodes
    nonzero <- which(data != 0, arr.ind = T)

    # So this bit identify the meta clusters that are brought forward from
    # previous day.
    # row = for each meta in current day
    for (row in 1:dim(data)[1]) {
      # sum(!is.na(data[row,])) = count number of previous day meta cluster this
      # current day meta cluster share with
      # if one, it means it's either a brought forward meta or split
      if (sum(!is.na(data[row,])) == 1) {
        # find the meta cluster from previous day it shares with
        column <- which(!is.na(data[row,]))
        # if equal to 1, it means that previous day meta cluster is not split.
        # otherwise there will be multiple current meta clusters i.e. sum is > 1.
        if (sum(!is.na(data[, column])) == 1) {
          new_data[row, column] <- names(nodes_per_clust_prev)[column]
        }
      } else if (sum(!is.na(data[row,])) > 1) {
        # find the meta cluster from previous day it shares with
        columns <- which(!is.na(data[row,]))
        for (column in columns) {
          new_data[row, column] <- names(nodes_per_clust_prev)[column]
        }
      }
    }

    # If row is empty, this is a new cluster (no affiliation with previous time step clusters), assign a new ID
    for (row in 1:length(nodes_per_clust)) {
      if (sum(is.na(data[row,])) == length(data[row,])) {
        get_new_id <- GetNewID(id_params$n, id_params$iter)
        new_ids[[row]] <- get_new_id$new_id
        id_params$n <- get_new_id$n
        id_params$iter <- get_new_id$iter
      }
    }

    # For each meta cluster in previous (split) and current (merge), count
    # the number of times it appears in the nonzero, i.e. number of meta
    # clusters which it shares with either previous or current.
    # For instance, merge_table:
    # 1 2 3 4
    # 0 0 0 2
    # it means meta cluster 1,2,3 in current time point does not share any nodes
    # with any meta in previous day, and meta cluster 4 share nodes with 2 meta
    # clusters in previous day.
    split_table <-
      table(factor(nonzero[, 2], levels = 1:length(nodes_per_clust_prev))) # columns


    # Creates list of all new meta clusters which are split
    # which is basically the row which have more than 2 columns which are not empty.
    split_list <- c()
    for (x in 1:length(split_table)) {
      if (split_table[x] > 1) {
        split_list <- c(split_list, x)
      }
    }

    # Create IDs for splits
    ranked_splits <- c()
    split_ids <- c()

    for (n in split_list) {
      # Sort by number of nodes
      sorted_split <-
        sort(
          data[, n],
          decreasing = T,
          index.return = T,
          na.last = T
        ) # Sort by number of nodes
      count <- sum(!is.na(sorted_split$x))
      ranked_splits[[n]] <- sorted_split$ix[1:count]

      # Check here if this ID has already split or not
      # If it has split, we will need to see how many times it has and add one
      # to it for the next split.
      # If it has not yet, we will need to split it and add it to the global
      # split list

      name <- names(nodes_per_clust_prev)[n]

      i <- CheckSplits(name, id_split_list)
      id <- c()
      for (j in 1:count) {
        if (j == 1) {
          id <- c(id, name)
        } else {
          id <- c(id, paste(name, i + 1, sep = "|"))
          i <- i + 1
          id_split_list[[name]] <- i
        }

      }
      split_ids[[n]] <- id
    }

    # Assign this ID to new data matrix
    for (i in 1:length(ranked_splits)) {
      if (!is.null(ranked_splits[[i]])) {
        n <- length(ranked_splits[[i]])
        for (m in 1:n) {
          new_data[ranked_splits[[i]][m], i] <- split_ids[[i]][m]
        }
      }
    }

    # Creates a list meta clusters which are either new, carried forward, or split
    for (i in 1:length(new_ids)) {
      if (is.null(new_ids[[i]])) {
        result <- c()
        merge_list <- which(!is.na(new_data[i,]))
        last_ids <- c()
        if (length(merge_list) == 1) {
          new_ids[[i]] <- new_data[i, merge_list[[1]]]

        }
      }
    }
    mc_new_noMerge <- list()
    for (node in c(1:length(mc_new))) {
      meta_cluster <- mc_new[node]
      # This mean the node is not included in this day
      if (is.na(meta_cluster)) {
        new_id <- NA
      } else {
        new_id <- new_ids[[meta_cluster]]
        if (is.null(new_id)) {
          # this means this is one of the merged meta cluster.

          # check if this meta split then merge.
          # if this is a split meta, we need to just get previously worked out id
          # do this by looking at the previous day's meta cluster for the node.
          # work out the index,
          # then check if there's id assigned to it
          prev_meta <- mc_prev[node]
          # then work out which index is this in the data (the numeric)
          idx_prev_meta <-
            which(names(nodes_per_clust_prev) == prev_meta)
          new_id <- new_data[meta_cluster, idx_prev_meta]

          # this is not a split cluster that merged
          # TODO maybe not needed
          if (is.null(new_id)) {
            # this is not a split node that merge, we assume previous day
            new_id <- as.character(mc_prev[node])

            # if this new_id is NA, it means the node wasn't part of previous
            # day. We need to assign new id.
            if (is.na(new_id)) {
              get_new_id <- GetNewID(id_params$n, id_params$iter)
              id_params$n <- get_new_id$n
              id_params$iter <- get_new_id$iter
              new_id <- get_new_id$new_id
            }
          } else if (length(new_id) == 0) {
            new_id <- NA
          }
        }
      }
      mc_new_noMerge[[node]] <- new_id
    }

    # if we have a meta cluster which receive new node, what happen is that the
    # ID for that node will be NA as this meta cluster may contain nodes from
    # various meta cluster in previous day,
    # we just assign it new cluster
    new_nodes <- which(is.na(mc_new_noMerge))

    # but if the node is not assigned a meta-cluster in this time-point,
    # the mc_new_noMerge entry will also be NA
    # so need to filter them out
    unassigned_nodes <- which(is.na(mc_new))
    new_nodes <- setdiff(new_nodes, unassigned_nodes)

    for (new_node in new_nodes) {
      get_new_id <- GetNewID(id_params$n, id_params$iter)
      mc_new_noMerge[[new_node]] <- get_new_id$new_id
      id_params$n <- get_new_id$n
      id_params$iter <- get_new_id$iter
    }

    mc_new_noMerge <- unlist(mc_new_noMerge)
    mc_new_noMerge <- as.factor(mc_new_noMerge)
    new_mcs[[timestep]] <- mc_new_noMerge
  }

  tracking <- list()

  # Create names for each timestep and return a list
  for (j in 1:length(new_mcs)) {
    tracking <-
      c(tracking, setNames(list(new_mcs[[j]]), paste("timestep", j, sep = "_")))

  }
  tracking

}

#' TrackingWithMerging
#'
#' Perform tracking of meta-clusters where merging is allowed.
#' Need to be run after FlowSOM's BuildSOM which contorts the SOM grid and
#' MetaClustering which perform meta-cluster of the SOM node.
#'
#' @param fsom list containing a FlowSOM object and meta-clustering lists.
#' @return A list containing N lists of unique IDs as determined by lineage,
#'   where N is no. of time steps.
#'
#' @import pdist
#' @import rlist
#'
#' @seealso \code{\link{TrackSOM}}.
#'
TrackingWithMerging <- function(fsom) {


  # Make lists that are needed for all time periods (clustering, iterations and
  # number)
  new_mcs <- list()
  id_params <- list("n"=1, "iter"=1)
  id_split_list <- list()


  # There is no tracking needed for first time step, but we need to get IDs from
  # this time step
  get_ids_1 <-
    GetFirstID(fsom$metaclustering$metaclustering[[1]],
               fsom$metaclustering$mcPerCell[[1]])
  new_mcs[[1]] <- get_ids_1$new_ids
  id_params$n <- get_ids_1$n
  id_params$iter <- get_ids_1$iter


  # Loop through remaining time steps and perform tracking
  for (timestep in 2:length(fsom$FlowSOM$metaData)) {

    # For the second time point, we want to base the IDs off the first time point
    if (timestep == 2) {
      mc_prev <- new_mcs[[1]]
    } else {
      mc_prev <- new_mcs[[timestep-1]]
    }

    mc <- fsom$metaclustering$metaclustering[[timestep]]
    mcPerCell <- (mc[fsom$FlowSOM$map$mapping[(fsom$FlowSOM$metaData[[timestep]][1]):(fsom$FlowSOM$metaData[[timestep]][2]),1]])

    # List of nodes per metacluster for this time step
    nodes_per_clust <- GetNodesPerCluster(mc)
    # List of nodes per metacluster from previous time step
    nodes_per_clust_prev <- GetNodesPerCluster(mc_prev)
    # Empty list to store IDs
    new_ids <- vector(mode = "list", length = length(nodes_per_clust))
    mc_new <- mc


    # Creates data matrix with number of overlapping nodes between metaclusters
    # from the two time steps
    data <-
      matrix(
        data = NA,
        nrow = length(nodes_per_clust),
        ncol = length(nodes_per_clust_prev)
      )
    for (m in 1:length(nodes_per_clust)) {
      for (n in 1:length(nodes_per_clust_prev)) {
        c <- length(which(nodes_per_clust_prev[[n]] %in% nodes_per_clust[[m]]))
        if (c > 0) {
          data[m,n] = c
        }
      }
    }

    # Create new matrix to hold the ID labels for easy merging and splitting of
    # ID labels later
    new_data <-
      matrix(
        data = NA,
        nrow = length(nodes_per_clust),
        ncol = length(nodes_per_clust_prev)
      )

    # Find all rows and columns with observations
    nonzero <- which(data !=0, arr.ind = T)

    # These tables count the number of occurences in a row (not number of nodes)
    # columns
    split_table <- table(factor(nonzero[,2], levels= 1:length(nodes_per_clust_prev)))
    # rows
    merge_table <- table(factor(nonzero[,1], levels= 1:length(nodes_per_clust)))

    # If both the columns and corresponding row have only 1 nonzero observation
    # In this case, the metacluster is same as last period, no splits or merges
    # so can assign it now.
    for (row in 1:dim(data)[1]) {
      if (sum(!is.na(data[row,]))== 1) {
        column <- which(!is.na(data[row,]))
        if (sum(!is.na(data[,column])) == 1) {
          new_data[row,column] <- names(nodes_per_clust_prev)[column]
        }

      }
    }

    # If row is empty, this is a new cluster (no affiliation with previous time
    # step clusters), assign a new ID
    for (row in 1:length(nodes_per_clust)) {
      if (sum(is.na(data[row,])) == length(data[row,])) {
        get_new_id <- GetNewID(id_params$n, id_params$iter)
        new_ids[[row]] <- get_new_id$new_id
        id_params$n <- get_new_id$n
        id_params$iter <- get_new_id$iter
      }
    }

    # Creates list of all columns that need to be split
    split_list <- c()
    for (x in 1:length(split_table)) {
      if (split_table[x] > 1) {
        split_list <- c(split_list, x)
      }
    }


    # Creates a list of all rows that need to be merged
    merge_list <- c()
    for (x in 1:length(merge_table)) {
      if (merge_table[x] > 1) {
        merge_list <- c(merge_list, x)
      }
    }


    # Create IDs for splits
    ranked_splits <- c()
    split_ids <- c()

    for (n in split_list) {
      # Sort by number of nodes
      sorted_split <- sort(data[,n], decreasing=T, index.return=T, na.last=T)
      count <- sum(!is.na(sorted_split$x))
      ranked_splits[[n]] <- sorted_split$ix[1:count]

      # Check here if this ID has already split or not If it has split, we will
      # need to see how many times it has and add one to it for the next split
      # If it has not yet, we will need to split it and add it to the global
      # split list

      name <- names(nodes_per_clust_prev)[n]

      i <- CheckSplits(name, id_split_list)
      id <- c()
      for (j in 1:count) {
        if (j ==1) {
          id <- c(id, name)
        } else {
          id <- c(id, paste(name, i+1, sep = "|"))
          i <- i+1
          id_split_list[[name]] <- i
        }

      }
      split_ids[[n]] <- id
    }

    # Assign this ID to new data matrix
    for (i in 1:length(ranked_splits)) {
      if (!is.null(ranked_splits[[i]])) {
        n <- length(ranked_splits[[i]])
        for (m in 1:n) {
          new_data[ranked_splits[[i]][m],i] <- split_ids[[i]][m]
        }
      }
    }

    # Find merges and rank them

    ranked_merges <- c()

    merge_ids <- c()

    for (n in merge_list) {
      sorted_merge <- sort(data[n,], decreasing=T, index.return=T, na.last=T)
      count <- sum(!is.na(sorted_merge$x))
      if (count > 0) {
        ranked_merges[[n]] <- sorted_merge$ix[1:count]
      }
    }



    # Now add the merge ids that are not already accounted for to new matrix

    for (i in 1:length(ranked_merges)) {
      if (!is.null(ranked_merges[[i]])) {
        n <- length(ranked_merges[[i]])
        for (m in 1:n) {
          if (is.na(new_data[i, ranked_merges[[i]][m]]))
            new_data[i, ranked_merges[[i]][m]] <-
              names(nodes_per_clust_prev)[ranked_merges[[i]][m]]
        }
      }
    }



    # Here we will fill in missing ids in the ID list with merges of the non-na
    # row IDs

    for (i in 1:length(new_ids)) {

      if (is.null(new_ids[[i]])) {
        result <- c()
        merge_list <- which(!is.na(new_data[i,]))
        last_ids <- c()
        if (length(merge_list) == 1) {
          new_ids[[i]] <- new_data[i,merge_list[[1]]]
        } else {
          for (m in merge_list) {
            last_ids <- c(last_ids, new_data[i,m])
          }

          last_ids_ <- unique(last_ids)

          if (length(last_ids_) == 1) {
            new_ids[[i]] <- last_ids_
          } else {
            ids <- paste(last_ids_, collapse= ",")

            new_ids[[i]] <- (paste("(", ids , sep="", ")"))
          }

        }

      }

    }



    # Assign the new ids to levels for output
    levels(mc_new)[1:length(nodes_per_clust)] <-
      c(unlist(new_ids))[1:length(nodes_per_clust)]
    new_mcs[[timestep]] <- mc_new


  }

  tracking <- list()

  # Create names for each timestep and return a list
  for (j in 1:length(new_mcs)) {
    tracking <-
      c(tracking, setNames(list(new_mcs[[j]]), paste("timestep", j, sep = "_")))

  }
  tracking

}

#' GetNewID
#'
#' Function which gets a new ID given the number in the alphabet is next and the
#' iteration in the alphabet.
#'
#' @param number current number of letters already issued.
#' @param iteration the current iteration, i.e. the number of times the letter
#'   has been repeated (AA then iteration = 2)
#' @return A list containing three elements:
#'   1) new_id: the new ID
#'   2) n: new number to be used as input for the next GetNewID function
#'   3) iter: new iteration to be used as input for the next GetNewID function
#'
GetNewID <- function(number, iteration) {

  let <- c() # create empty letters list
  n_iteration <- iteration

  # Fill in letters list
  for (alphabet in LETTERS) {
    let <- c(let, strrep(alphabet, n_iteration))
  }

  new_id <- let[number]

  n <- number +1
  if (n > 26) {
    n_iteration <- iteration +1
    n <- 1
  }

  new_list <- list("new_id"=new_id, "n"=n, "iter"=n_iteration)
  return(new_list)
}

#' GetFirstID
#'
#' Function which obtains the IDs for the first time point in tracking.
#'
#' @param mc the meta-cluster ID for each node
#' @param mcPerCell metaclustering per datapoint
#'
#' @return A list containing three elements:
#'  1) new_ids: the new IDs
#'  2) n: new number to be used as input for GetNewID function
#'  3) iter: new iteration to be used as unput to GetNewID function
#'
GetFirstID <- function(mc, mcPerCell) {

  let <- LETTERS
  iteration <- 1
  j<-1


  for (i in 1:length(levels(mc))) {
    levels(mc)[[i]] <- let[j]
    j<- j+1

    if (j == 27) {
      iteration <- iteration +1

      for (alphabet in LETTERS) {
        let <- c(strrep(LETTERS, iteration))
        j <- 1
      }
    }
  }

  list("new_ids"=mc, "n"=j, "iter"=iteration)
}


#' GetNodesPerCluster
#'
#' Function which generates a list of the nodes for each metacluster, used to
#' generate data matrix for tracking.
#'
#' @param mc A factor of the meta-cluster per node
#'
#' @return A list of node IDs for each metacluster ID
#'
GetNodesPerCluster <- function(mc) {
  data <- list()

  for (id in unique(unlist(mc))) {
    if (!is.na(id)) {
      data[[id]] <- which(mc %in% id)
    }

  }
  data
}



#' CheckSplits
#'
#' A function to check number of splits and return number of splits
#'
#' @param id The meta-cluster id to check
#' @param id_split_list The list of meta-clusters currently exist
#'
#' @return The number of splits for the meta-cluster
#'
CheckSplits <- function(id, id_split_list) {

  if (id %in% names(id_split_list)) {
    n_splits <- id_split_list[[id]]
  } else {
    n_splits <- 0
  }
  n_splits
}


