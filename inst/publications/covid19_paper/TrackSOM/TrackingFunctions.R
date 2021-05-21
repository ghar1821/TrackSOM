
# All Tracking functions here

# Import libraries
library(pdist) # for distances
library(rlist) # Needed fot GetTransitions function

# --------------------------------------------------------------------------------------------------------------------

# Input: a TrackSOM object containing the FlowSOM and metaclustering lists
# Output: a list containing N lists of unique IDs as determined by lineage, where N is no. of time steps

TrackingByLineage <- function(fsom) {
  
  
  # Make lists that are needed for all time periods (clustering, iterations and number)
  new_mcs <- list()
  id_params <- list("n"=1, "iter"=1)
  id_split_list <- list()
  
  
  # There is no tracking needed for first time step, but we need to get IDs from this time step 
  get_ids_1 <- GetFirstID(fsom$metaclustering$metaclustering[[1]], fsom$metaclustering$mcPerCell[[1]])
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
    
    nodes_per_clust <- GetNodesPerCluster(mc) # List of nodes per metacluster for this time step
    nodes_per_clust_prev <- GetNodesPerCluster(mc_prev) # List of nodes per metacluster from previous time step
    new_ids <- vector(mode = "list", length = length(nodes_per_clust)) # Empty list to store IDs
    mc_new <- mc 
    
    
    # Creates data matrix with number of overlapping nodes between metaclusters from the two time steps
    data <- matrix(data=NA, nrow=length(nodes_per_clust), ncol=length(nodes_per_clust_prev))
    for (m in 1:length(nodes_per_clust)) {
      for (n in 1:length(nodes_per_clust_prev)) {
        c <- length(which(nodes_per_clust_prev[[n]] %in% nodes_per_clust[[m]]))
        if (c > 0) {
          data[m,n] = c
        }
      }
    }
    
    # Create new matrix to hold the ID labels for easy merging and splitting of ID labels later
    new_data <- matrix(data=NA, nrow=length(nodes_per_clust), ncol=length(nodes_per_clust_prev))
    
    # Find all rows and columns with observations
    nonzero <- which(data !=0, arr.ind = T)
    
    # These tables count the number of occurences in a row (not number of nodes)
    split_table <- table(factor(nonzero[,2], levels= 1:length(nodes_per_clust_prev))) # columns
    merge_table <- table(factor(nonzero[,1], levels= 1:length(nodes_per_clust))) # rows
    
    # If both the columns and corresponding row have only 1 nonzero observation
    # In this case, the metacluster is same as last period, no splits or merges so can assign it now.
    for (row in 1:dim(data)[1]) {
      if (sum(!is.na(data[row,]))== 1) {
        column <- which(!is.na(data[row,]))
        if (sum(!is.na(data[,column])) == 1) {
          new_data[row,column] <- names(nodes_per_clust_prev)[column]
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
      sorted_split <- sort(data[,n], decreasing=T, index.return=T, na.last=T) # Sort by number of nodes
      count <- sum(!is.na(sorted_split$x))
      ranked_splits[[n]] <- sorted_split$ix[1:count]
      
      # Check here if this ID has already split or not
      # If it has split, we will need to see how many times it has and add one to it for the next split
      # If it has not yet, we will need to split it and add it to the global split list
      
      name <- names(nodes_per_clust_prev)[n]
      
      i <- check_splits(name, id_split_list)
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
            new_data[i, ranked_merges[[i]][m]] <- names(nodes_per_clust_prev)[ranked_merges[[i]][m]]
        }
      }
    }
    
    
    
    # Here we will fill in missing ids in the ID list with merges of the non-na row IDs
    
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
    levels(mc_new)[1:length(nodes_per_clust)] <- c(unlist(new_ids))[1:length(nodes_per_clust)]
    new_mcs[[timestep]] <- mc_new
    
    
  }
  
  tracking <- list()
  
  # Create names for each timestep and return a list
  for (j in 1:length(new_mcs)) {
    tracking <- c(tracking, setNames(list(new_mcs[[j]]), paste("timestep", j, sep="_")))
    
  }
  tracking
  
}

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------------

# Inputs: 
#   1) A TrackSOM object containing FlowSOM and metaclustering lists
#   2) Output from TrackingByLineage function
# Output: A list of N lists containing the unique IDs as determined by historical proximity for each of N time steps

TrackingByHistoricalProximity <- function(fsom, tracking_by_lineage) {
  
  n_nodes <- length(c(tracking_by_lineage[[1]]))
  vector <- numeric(length = n_nodes)
  new_mcs <- list()
  
  
  # For each time step beginning with timestep 2:
  for (timestep in 2:length(fsom$FlowSOM$metaData)) {
    
    clust_list <- list()
    mc_new <- fsom$metaclustering$metaclustering[[timestep]]
    
    # Create distance matrix between pairs of codes from
    matrix <- as.matrix(pdist(fsom$FlowSOM$map$coding[[timestep]], fsom$FlowSOM$map$coding[[timestep-1]]))
    matrix[matrix == 0] <- NA # if == 0, this means that there were no data points
    
    
    nodes_per_clust <- GetNodesPerCluster(tracking_by_lineage[[timestep]])
    
    
    for (i in 1:n_nodes) {
      min <- which.min(matrix[i,])
      if (length(min) > 0) {
        vector[i] <- min
      } else {
        vector[i] <- NA
      }
      
    }
    
    cluster_per_node <- as.factor(tracking_by_lineage[[timestep-1]][vector])
    
    
    for (clust in 1:length(nodes_per_clust)) {
      names <- c()
      for (node in nodes_per_clust[[clust]]) {
        names <- c(names, as.character(cluster_per_node[[node]]))
      }
      
      unique_names <- unique(names)
      if (length(unique_names) > 1) {
        names_ <- paste(unique_names, collapse= " & ")
      } else {
        names_ <- unique_names
      }
      clust_list[[clust]] <- names_
    }
    levels(mc_new)[1:length(nodes_per_clust)] <- c(unlist(clust_list))[1:length(nodes_per_clust)]
    
    new_mcs[[timestep]] <- mc_new
    
    
  }
  
  new_mcs
  
}


# --------------------------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------------------------




# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------------

# Function which gets a new ID given the number in the alphabet is next and the iteration in the alphabet
# Inputs:
#   1) number - current number of letter in alphabet
#   2) iteration - current iteration
# 
# Output: a list containing three lists (new_id, n and iter)
#   1) new_id - the new ID
#   2) n - new number to be used as input for GetNewID function
#   3) iter - new iteration to be used as unput to GetNewID function

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


# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------------

# Function which obtains the IDs for the first time point in tracking
# Inputs:
#   1) mc - metaclustering per node
#   2) mcPerCell - metaclustering per datapoint
# 
# Output: a list containing three lists (new_id, n and iter)
#   1) new_ids - the new IDs
#   2) n - new number to be used as input for GetNewID function
#   3) iter - new iteration to be used as unput to GetNewID function

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


# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------------

# Function which generates a list of the nodes for each metacluster, used to generate data matrix for tracking
#   Input:    factor of metaclustering per node
#   Output:   list of node IDs for each metacluster ID

GetNodesPerCluster <- function(mc) {
  data <- list()
  
  for (id in unique(unlist(mc))) {
    if (!is.na(id)) {
      data[[id]] <- which(mc %in% id)
    }
    
  }
  data
}

# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------------

## Function to check number of splits and return number of splits

check_splits <- function(id, id_split_list) {
  
  if (id %in% names(id_split_list)) {
    n_splits <- id_split_list[[id]]
  } else {
    n_splits <- 0
  }
  n_splits
}



