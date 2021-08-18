#' DrawNetworkPlot
#'
#' Draw TrackSOM network plot where each node represents a meta-cluster,
#' coloured by the time-point it exists or median/mean marker values.
#' Meta-clusters are connected with edges, indicating a transition over time.
#'
#' @param dat Data table. The dataset
#' @param timepoint.col Column denoting the time point column
#' @param timepoints Vector of time-points
#' @param cluster.col Column denoting the TrackSOM meta-cluster ID
#' @param marker.cols Columns denoting the markers to plot
#' @param calculation.type The aggregation (mean/median) type for colouring the
#'   nodes
#' @param node.size Whether the nodes are to be sized based on the proportion of
#'   cells it captured ('auto') or just constant (a number denoting the node
#'   size)
#' @param min.node.size Smallest node size. Only used if node.size is 'auto'
#' @param max.node.size Largest node size. Only used if node.size is 'auto'
#' @param arrow.length Length of the arrow
#' @param arrow.head.gap Gap between the node and the head of the arrow
#' @param standard.colours Colour scheme for colouring the nodes. Can only be
#'   either one of the colours specified by the viridis package
#'   \url{https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html}
#'   or spectral, bupu, rdylbu.
#' @param img.height Height of the image for network plot export
#' @param img.width Width of the image for network plot export
#' @param load.temp.data Whether to load some intermediate calculations or not.
#'   These intermediate files contain the details of each node and edge, and
#'   produced by the function after the very first run.
#' @param mapping A data.table containing the annotation (cellular population)
#'   each meta-cluster in each time-point represents. The data.table must
#'   contain a column meta-cluster ID column named as per cluster.col,
#'   time-point column named as per timepoint.col, and a cell population column
#'   named as per population.col.
#' @param population.col Column denoting the cell population name column. Only
#'   used if mapping is not NULL.
#' @param graph.layout The layout of the network plot. Could be anything from
#'   the ggraph, graphlayouts, and igraph packages.
#' @param legend.position Position of the legend, as per ggplot2.
#' @param no_merge Whether meta-cluster merging was allowed. This only affect
#'   the network plot coloured by the meta-cluster origin. If merging was
#'   permitted, any merged meta-clusters will be coloured grey. If merging was
#'   not permitted, any split meta-clusters will be coloured the same as their
#'   parent.
#' @param file.format Image export file format
#' @param line.width Width of each edge
#'
#' @usage DrawNetworkPlot(dat, timepoint.col, timepoints, cluster.col,
#'  marker.cols)
#'
#' @import gtools
#' @import tidygraph
#' @import ggraph
#' @import igraph
#' @import RColorBrewer
#' @import viridis
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
#' DrawNetworkPlot(tracksom_res,
#'                 'timepoint',
#'                 c("Mock", "SYN-1", "SYN-2", "SYN-3", "SYN-4"),
#'                 "TrackSOM_metacluster_lineage_tracking",
#'                 use_cols)
#'
#' @export
#'
DrawNetworkPlot <- function(dat,
                            timepoint.col,
                            timepoints,
                            cluster.col,
                            marker.cols,
                            calculation.type = c('mean', 'median'),
                            node.size = 'auto',
                            min.node.size = 6,
                            max.node.size = 20,
                            arrow.length = 3,
                            arrow.head.gap = 4,
                            standard.colours = 'spectral',
                            img.height = 9,
                            img.width = 9,
                            load.temp.data = FALSE,
                            mapping = NULL,
                            population.col = NULL,
                            graph.layout = 'kk',
                            legend.position = 'right',
                            no_merge = FALSE,
                            file.format = 'pdf',
                            line.width = 1) {

  set.seed(42)

  # TODO change graph.layout to vector to restrict the "acceptable" layout

  calculation.type <- match.arg(calculation.type)

  ## Have to make sure column header have no hyphen
  if (TRUE %in% stringr::str_detect(colnames(dat), '-')) {
    message("Some column headers have hyphen (-) in it.
            Please rename them first!")
    message("No plots are created.")
    return()
  }

  if (isFALSE(load.temp.data)) {
    message("Calculating edges")
    # To store transitions
    edge.df <-
      GetTransitionsAsEdges(dat, timepoints, timepoint.col, cluster.col)

    message("Computing node details")
    # Then get extra details on the nodes
    # Get all cluster ID and concatenate them with the time-point they belong to
    # E.g. 1_A, 2_(A,C), etc.
    all.clust <- lapply(c(1:length(timepoints)), function(tp.idx) {
      tp.dat <- dat[dat[[timepoint.col]] == timepoints[tp.idx], ]
      tp.clust <- mixedsort(unique(tp.dat[[cluster.col]]))
      tp.clust.name <- sapply(tp.clust, function(c) {
        return(paste0(tp.idx, '_', c))
      })
      return(tp.clust.name)
    })
    all.clust <- unlist(all.clust)

    cluster.ids <- sapply(all.clust, function(c) {
      strsplit(c, '_')[[1]][2]
    })
    timepoints.idx <- sapply(all.clust, function(c) {
      strsplit(c, '_')[[1]][1]
    })
    timepoints.actual <- sapply(timepoints.idx, function(c) {
      timepoints[as.numeric(c)]
    })
    node.dat <- data.table(nodeId = all.clust,
                           clusterId = cluster.ids,
                           timepointIdx = timepoints.idx,
                           timepoint = timepoints.actual)

    message("Calculating marker's average per node")
    grp_dat <- c(cluster.col, timepoint.col)

    # Compute each cluster's mean or median expression
    if (calculation.type == 'mean') {
      average_marker_expressions <-
        dat[, lapply(.SD, mean), by = grp_dat, .SDcols = marker.cols]
    } else if (calculation.type == 'median') {
      average_marker_expressions <-
        dat[, lapply(.SD, median), by = grp_dat, .SDcols = marker.cols]
    }

    node.dat <- merge(node.dat, average_marker_expressions,
                      by.x = c("clusterId", "timepoint"),
                      by.y = grp_dat)

    # Add an extra option which will colour just the node using the very first
    # cluster.We only want clusters which is alphabetical, the very first
    # cluster
    cluster.origins <- GetClusterOrigins(node.dat, no_merge)
    node.dat$origin <- as.factor(cluster.origins)

    # compute node sizes
    if (node.size == 'auto') {
      node.sizes <- GetClusterProportions(node.dat = node.dat,
                                            dat = dat,
                                            timepoint.col = timepoint.col,
                                            timepoints = timepoints,
                                            cluster.col = cluster.col)

      node.dat$ProportionOfCells <- node.sizes
    }

    # if population mapping is given, then we map
    if (!is.null(mapping)) {
      node.dat <-
        merge(node.dat,
              mapping,
              by.x = c("clusterId", "timepoint"),
              by.y = grp_dat)
    }

    message("Saving node and edge details")
    fwrite(node.dat, "network_plot_node_details.csv")
    fwrite(edge.df, "network_plot_edge_details.csv")

  } else {
    message("Reading node and edge details")
    node.dat <- fread("network_plot_node_details.csv")
    edge.df <- fread("network_plot_edge_details.csv")

    # handle na
    indx <- which(sapply(node.dat, is.character))
    for (j in indx)
      set(
        node.dat,
        i = grep("^$|^ $", node.dat[[j]]),
        j = j,
        value = NA_character_
      )

  }
  # order node dat based on the timepoints
  node.dat$timepoint <- factor(node.dat$timepoint, levels=timepoints)

  setcolorder(node.dat, c("nodeId", setdiff(names(node.dat), "nodeId")))

  # Start plotting

  message("Start drawing plots")

  # Setup layout
  g <- graph_from_data_frame(edge.df, directed = TRUE, node.dat)
  lay <- create_layout(g, layout = graph.layout)

  # Draw by time point
  message("Drawing plots coloured by time point")

  geompoint.prop <- GetGeomPoint(node.size, 'timepoint')

  colour.palette <- GetColourPallete(standard.colours,
                                       n.col = length(timepoints),
                                       factor = timepoints)

  plt <- ggraph(lay) +
    geom_edge_link(arrow = arrow(length = unit(arrow.length, 'mm')),
                   end_cap = circle(arrow.head.gap, 'mm'),
                   width = line.width) +
    geompoint.prop +
    colour.palette

  if (node.size == 'auto') {
    plt <-
      plt + scale_size_continuous(range = c(min.node.size, max.node.size)) +
      guides(colour = guide_legend(override.aes = list(size = min.node.size))) +
      labs(color = 'Time-point', size = 'Proportion')
  } else {
    plt <- plt + labs(color='Time-point')
  }

  plt <- plt + theme(text = element_text(size=20),
                     panel.background = element_rect(fill = 'white'),
                     aspect.ratio = 1,
                     legend.position = legend.position,
                     plot.margin = grid::unit(c(0,0,0,0), "mm"))
  ggsave(paste0("network_colBy_timepoints.", file.format),
         plot = plt,
         width = img.width,
         height = img.height,
         dpi = 1000,
         limitsize = FALSE)

  # Draw by the origin cluster
  # +1 for the colour brewer because we have NA

  message("Drawing plots coloured by origin")

  unique.origins <- as.character(mixedsort(unique(node.dat$origin)))
  colour.palette <- GetColourPallete(standard.colours,
                                       n.col = length(unique.origins),
                                       factor = unique.origins)
  geompoint.prop <- GetGeomPoint(node.size, 'origin')

  plt <- ggraph(lay) +
    geom_edge_link(arrow = arrow(length = unit(arrow.length, 'mm')),
                   end_cap = circle(arrow.head.gap, 'mm'),
                   width = line.width) +
    geompoint.prop +
    colour.palette

  if (node.size == 'auto') {
    plt <- plt + scale_size_continuous(range = c(min.node.size, max.node.size)) +
      guides(colour = guide_legend(override.aes = list(size=min.node.size))) +
      labs(color='Cluster origin', size='Proportion')
  } else {
    plt <- plt + labs(color='Cluster origin')
  }

  plt <- plt + theme(text = element_text(size=20),
                     panel.background = element_rect(fill = 'white'),
                     aspect.ratio = 1,
                     legend.position = legend.position,
                     plot.margin = grid::unit(c(0,0,0,0), "mm"))
  ggsave(paste0("network_colBy_origin.", file.format),
         plot = plt,
         width = img.width,
         height = img.height,
         dpi = 1000,
         limitsize = FALSE)

  # Draw by population if given
  if (!is.null(mapping)) {
    message(paste("Drawing plots coloured by", population.col))

    unique.populations <-
      as.character(mixedsort(unique(node.dat[[population.col]])))
    colour.palette <- GetColourPallete(standard.colours,
                                         n.col = length(unique.populations),
                                         factor = unique.populations)
    geompoint.prop <- GetGeomPoint(node.size, population.col)

    plt <- ggraph(lay) +
      geom_edge_link(arrow = arrow(length = unit(arrow.length, 'mm')),
                     end_cap = circle(arrow.head.gap, 'mm'),
                     width = line.width) +
      geompoint.prop +
      colour.palette

    if (node.size == 'auto') {
      plt <- plt + scale_size_continuous(range = c(min.node.size, max.node.size)) +
        guides(colour = guide_legend(override.aes = list(size=min.node.size))) +
        labs(color=population.col, size='Proportion')
    } else {
      plt <- plt + labs(color=population.col)
    }

    plt <- plt + theme(text = element_text(size=20),
                       panel.background = element_rect(fill = 'white'),
                       aspect.ratio = 1,
                       legend.position = legend.position,
                       plot.margin = grid::unit(c(0,0,0,0), "mm"))
    ggsave(paste0("network_colBy_", population.col, ".", file.format),
           plot = plt,
           width = img.width,
           height = img.height,
           dpi = 1000,
           limitsize = FALSE)
  }

  # Colour plot by markers
  colour.palette <- GetColourPallete(standard.colours)

  for (marker in marker.cols) {

    message(paste0("Drawing plots coloured by ", marker))

    geompoint.prop <- GetGeomPoint(node.size, marker)

    plt <- ggraph(lay) +
      geom_edge_link(arrow = arrow(length = unit(arrow.length, 'mm')),
                     end_cap = circle(arrow.head.gap, 'mm'),
                     width = line.width) +
      geompoint.prop +
      colour.palette +
      labs(size='Proportion')

    if (node.size == 'auto') {
      plt <- plt + scale_size_continuous(range = c(min.node.size, max.node.size))
    }

    plt <- plt + theme(text = element_text(size=20),
                       panel.background = element_rect(fill = 'white'),
                       aspect.ratio = 1,
                       legend.position = legend.position,
                       plot.margin = grid::unit(c(0,0,0,0), "mm"))
    ggsave(paste0("network_colBy_", marker, '.', file.format),
           plot = plt,
           width = img.width,
           height = img.height,
           dpi = 1000,
           limitsize = FALSE)
  }

}

#' GetColourPallete
#'
#' Function to get colour pallete for network plot's nodes
#'
#' @param standard.colours the colour scheme
#' @param n.col number of colours to get
#' @param factor discrete values for each colour
#'
#' @return colour mapping
#'
GetColourPallete <- function(standard.colours,
                             n.col = 50, factor = NULL) {

  if(tolower(standard.colours) == "spectral"){
    colour.list <-
      rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(n.col))
    if (is.null(factor)) {
      colour.palette <- scale_colour_gradientn(colours = colour.list)
    } else {
      colour.palette <- scale_colour_manual(values = colour.list,
                                            breaks = factor,
                                            na.value='grey')
    }
  }
  else if (tolower(standard.colours) %in% c('bupu', 'rdylbu')) {
    if(tolower(standard.colours) == "bupu"){
      colour.list <-
        rev(colorRampPalette(RColorBrewer::brewer.pal(9, "BuPu"))(n.col))
    } else {
      colour.list <-
        rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu"))(n.col))

    }
    if (is.null(factor)) {
      colour.palette <- scale_colour_gradientn(colours = colour.list)
    } else {
      colour.palette <- scale_colour_manual(values = colour.list,
                                            breaks = factor,
                                            na.value='grey')
    }
  }
  else {
    if (is.null(factor)) {
      colour.palette <- scale_color_viridis(option=tolower(standard.colours))
    } else {
      colour.palette <-
        scale_color_viridis(
          discrete = TRUE,
          option = tolower(standard.colours),
          na.value = 'grey'
        )
    }
  }
  return(colour.palette)
}

#' GetClusterOrigins
#'
#' Function to infer the origins of a meta-cluster
#'
#' @param node.dat data.table containing the nodes details
#' @param no_merge whether meta-cluster merging was enabled or not
#'
#' @return a vector of meta-cluster origins
#'
GetClusterOrigins <- function(node.dat,
                              no_merge) {

  clusters <- as.vector(unique(node.dat$clusterId))
  cluster.origins <- sapply(clusters, function(cl.id) {
    alphabet.only <- grepl('^[A-Z]+$', cl.id)
    if (alphabet.only) {
      return(cl.id)
    } else {
      if (!no_merge) {
        return('NotOrigin')
      } else {
        # if the user decided to use this with merging, then it'll just be
        # coloured by one of the meta cluster
        return(strsplit(cl.id, '|')[[1]][1])
      }
    }
  })

  # remove the cluster which is not alphabet only (as above function assign it
  # to null)
  cluster.origins <- cluster.origins[cluster.origins != 'NotOrigin']

  cluster.colours <- sapply(c(1:nrow(node.dat)), function(i) {
    row <- node.dat[i,]
    cluster.remain <- row$clusterId %in% names(cluster.origins)

    if (cluster.remain) {
      return(as.character(cluster.origins[row$clusterId]))
    } else {
      return(NA)
    }
  })
}

#' GetClusterProportions
#'
#' Function to count the proportin of cells in each meta-cluster
#'
#' @param node.dat data.table containing the node details
#' @param dat the raw clustered and tracked data
#' @param timepoint.col the name of the column denoting the time-point
#' @param timepoints a vector of time-points available in the data
#' @param cluster.col the name of the column denoting the meta-cluster ID
#'
#' @return a vector of meta-clusters' proportion
#'
GetClusterProportions <- function(node.dat,
                                  dat,
                                  timepoint.col,
                                  timepoints,
                                  cluster.col) {


  # count cell per time point
  cell.count.per.tp <- dat[, .N, by=dat[[timepoint.col]]]

  proportions <- sapply(c(1:nrow(node.dat)), function(i) {
    node.row <- node.dat[i,]
    sub.dat <- dat[dat[[timepoint.col]] == node.row$timepoint &
                     dat[[cluster.col]] == node.row$clusterId, ]
    cell.cnt <- cell.count.per.tp[cell.count.per.tp$dat == node.row$timepoint, N]
    return(nrow(sub.dat)/cell.cnt)
  })


  return(proportions)
}

#' GetGeomPoint
#'
#' Function to get node object from
#'
#' @param node.size whether the nodes are sized based on proportion of cells in
#'   the meta-clusters or constant number
#' @param col.by what to colour the node based on
#'
#' @return a geom node point object
#'
GetGeomPoint <- function(node.size, col.by) {
  if (node.size == 'auto') {
    geompoint <-
      geom_node_point(aes_string(colour = col.by, size = 'ProportionOfCells'))
  }
  else {
    geompoint <- geom_node_point(aes_string(colour = col.by), size=node.size)
  }
  return(geompoint)
}
