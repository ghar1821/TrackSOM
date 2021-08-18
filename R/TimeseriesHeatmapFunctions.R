#' DrawTimeseriesHeatmap
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
#' @param colours Colour scheme for colouring the nodes. Can only be
#'   either one of the colours specified by the viridis package
#'   \url{https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html}
#'   or spectral, bupu, rdylbu.
#' @param plot.height Height of the image for network plot export
#' @param plot.width Width of the image for network plot export
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
#' @param map.cluster.id Whether to map the meta-cluster ID to numeric. Handy if
#'   the meta-cluster IDs get too long.
#' @param cluster.ordering A custom ordering for the meta-clusters (y-axis)
#' @param legend.position Position of the legend, as per ggplot2.
#' @param no_merge Whether meta-cluster merging was allowed. This only affect
#'   the network plot coloured by the meta-cluster origin. If merging was
#'   permitted, any merged meta-clusters will be coloured grey. If merging was
#'   not permitted, any split meta-clusters will be coloured the same as their
#'   parent.
#' @param file.format Image export file format
#' @param line.width Width of each edge
#' @param font.size Size of the font
#'
#' @usage DrawTimeseriesHeatmap(dat, timepoint.col, timepoints, cluster.col,
#'  marker.cols)
#'
#' @import ggplot2
#' @import viridis
#' @import data.table
#' @import stringr
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
#' DrawTimeseriesHeatmap(tracksom_res,
#'                       'timepoint',
#'                       c("Mock", "SYN-1", "SYN-2", "SYN-3", "SYN-4"),
#'                       "TrackSOM_metacluster_lineage_tracking",
#'                       use_cols)
#'
#' @export
#'
DrawTimeseriesHeatmap <- function(dat,
                                    timepoints,
                                    timepoint.col,
                                    cluster.col,
                                    marker.cols,
                                    colours = 'viridis',
                                    plot.width = 12,
                                    plot.height = 15,
                                    font.size = 16,
                                    file.format = 'pdf',
                                    min.node.size = 6,
                                    max.node.size = 20,
                                    plot.title = NULL,
                                    legend.position = 'right',
                                    load.temp.data = FALSE,
                                    arrow.length = 2,
                                    mapping = NULL,
                                    population.col = NULL,
                                    cluster.ordering = NULL,
                                    map.cluster.id = FALSE,
                                    calculation.type = c('mean', 'median'),
                                    line.width = 1) {



  if (!load.temp.data) {
    message("Computing node details")

    # So we need a data table with day, cluster id, centroid, proportion column
    # for just the bubble plot
    grp <- c(timepoint.col, cluster.col)

    # Get the centroids and cell proportions ====
    # https://stackoverflow.com/questions/50710266/r-data-table-sum-all-columns-whose-names-are-stored-in-a-vector
    calculation.type = match.arg(calculation.type)

    if (calculation.type == 'mean') {
      centroids <- dat[, lapply(.SD, mean), by = grp, .SDcols = marker.cols]
    } else if (calculation.type == 'median') {
      centroids <- dat[, lapply(.SD, median), by = grp, .SDcols = marker.cols]
    }

    cellCount <- dat[, .(count = .N), by = grp]
    cellCount_perDay <- dat[, .(countPerDay = .N), by = timepoint.col]
    cellCount <- merge(cellCount, cellCount_perDay, by=timepoint.col, all.x = TRUE)
    cellCount[, proportion:=(count/countPerDay)]

    dat_chart <- merge(centroids, cellCount, by=grp)

    # append the mapping of cluster and population name.
    if(!is.null(mapping)) {
      dat_chart <- merge(dat_chart, mapping, by=grp)
    }

    message("Computing edge details")

    # Get the edges to link the clusters ====
    # the edges need network plot to work out edges
    edges <- GetTransitionsAsEdges(dat = dat,
                                      timepoints = timepoints,
                                      timepoint_col = timepoint.col,
                                      cluster_col = cluster.col)
    edges <- data.table(edges)
    # have to convert this to data table with columns: id, timepoint?
    from_tp <- sapply(edges$from, function(fr) {
      tp <- as.numeric(str_split_fixed(fr, "_", 2)[1,1])
      timepoints[tp]
    })
    to_tp <- sapply(edges$to, function(to) {
      tp <- as.numeric(str_split_fixed(to, "_", 2)[1,1])
      timepoints[tp]
    })
    from_cl <- sapply(edges$from, function(fr) {
      str_split_fixed(fr, "_", 2)[1,2]
    })
    to_cl <- sapply(edges$to, function(to) {
      str_split_fixed(to, "_", 2)[1,2]
    })
    edges_chart <- data.table(timepoint=c(from_tp, to_tp))
    edges_chart$cluster <- c(from_cl,to_cl)
    edges_chart$group <- rep(c(1:length(from_tp)),2)

    message("Saving node and edge details")
    fwrite(dat_chart, "timeseries_heatmap_node_details.csv")
    fwrite(edges_chart, "timeseries_heatmap_edges_details.csv")
  } else {
    message("Reading node and edge details")
    dat_chart <- fread("timeseries_heatmap_node_details.csv")
    edges_chart <- fread("timeseries_heatmap_edges_details.csv")
  }

  # order the data based on the timepoints given
  # convert ID to numeric if need be
  if (map.cluster.id) {
    if (!is.null(cluster.ordering)) {
      sorted_unique_clust <- cluster.ordering
    } else {
      sorted_unique_clust <- sort(unique(dat[[cluster.col]]))
    }
    unique_clust <- sapply(c(1:length(sorted_unique_clust)), function(x)
      as.character(x))
    names(unique_clust) <- sorted_unique_clust
    dat_chart$mapped_cluster_id <-
      sapply(dat_chart[[cluster.col]], function(cl) {
        unique_clust[cl]
      })
    dat_chart$mapped_cluster_id <-
      factor(dat_chart$mapped_cluster_id, levels = unique_clust)
    edges_chart$mapped_cluster_id <- sapply(edges_chart$cluster, function(cl) {
      unique_clust[cl]
    })

    to_plot_node <- 'mapped_cluster_id'
    to_plot_edge <- 'mapped_cluster_id'

    mapped_dat<- data.table(cluster=names(unique_clust),
                            numeric_id=unique_clust)
    setnames(mapped_dat, 'cluster', cluster.col)
    fwrite(mapped_dat, paste0(cluster.col, '_numericMapping.csv'))
  } else {
    # if don't need to convert then order the data based on the cluster ordering
    # if given
    if (!is.null(cluster.ordering)) {
      dat_chart[[cluster.col]] <- factor(dat_chart[[cluster.col]], levels=cluster.ordering)
    } else {
      # well if cluster ordering is not given, then we just do lexicographical
      # sort
      dat_chart[[cluster.col]] <-
        factor(dat_chart[[cluster.col]], levels = sort(unique(dat_chart[[cluster.col]])))
    }
    to_plot_node <- cluster.col
    to_plot_edge <- 'cluster'
  }

  # Draw plots ====
  if (colours == 'spectral') {
    col.scheme <- colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))
  } else {
    col.scheme <- colorRampPalette(c(viridis_pal(option = colours)(50)))
  }


  for (marker in marker.cols) {
    message(paste("Drawing timeseries heatmap coloured by", marker))
    plt <-
      ggplot(dat_chart, aes_string(x = timepoint.col, y = to_plot_node)) +
      geom_point(aes_string(size = "proportion", fill = marker),
                 alpha = 0.75,
                 shape = 21) +
      geom_line(
        aes_string(x = "timepoint",
                   y = to_plot_edge, group = "group"),
        edges_chart,
        arrow = arrow(length = unit(as.numeric(arrow.length), "mm")),
        size = line.width
      ) +
      scale_size_continuous(limit = c(min(dat_chart$proportion), max(dat_chart$proportion)),
                            range = c(min.node.size, max.node.size)) +
      labs(x = NULL,
           y = NULL,
           size = "Proportion of cells",
           fill = marker)  +
      theme(
        legend.key = element_blank(),
        axis.text.x = element_text(
          colour = "black",
          size = font.size,
          angle = 90,
          vjust = 0.3,
          hjust = 1
        ),
        axis.text.y = element_text(colour = "black", size = font.size),
        legend.text = element_text(size = font.size, colour = "black"),
        legend.title = element_text(size = font.size),
        panel.background = element_blank(),
        panel.border = element_rect(
          colour = "black",
          fill = NA,
          size = 1.2
        ),
        legend.position = legend.position,
        plot.title = element_text(
          color = "Black",
          face = "bold",
          size = font.size,
          hjust = 0
        ),
        plot.margin = grid::unit(c(0, 0, 0, 0), "mm")
      ) +
      scale_fill_gradientn(colours = c(col.scheme(50))) +
      scale_y_discrete(limits = rev(levels(factor(dat_chart[[to_plot_node]]))))

    if (!(is.null(plot.title))) {
      plt <- plt + ggtitle(plot.title)
    }

    if (map.cluster.id) {
      filename <- paste0("Timeseries_heatmap_mapped_by_", marker, ".", file.format)
    } else {
      filename <- paste0("Timeseries_heatmap_by_", marker, ".", file.format)
    }

    ggsave(filename = filename,
           plot = plt,
           width = plot.width,
           height = plot.height,
           limitsize = FALSE)
  }

  if(!is.null(population.col)) {

    num_populations <- length(unique(dat_chart[[population.col]]))

    message(paste("Drawing timeseries heatmap coloured by", population.col))
    plt <-
      ggplot(dat_chart, aes_string(x = timepoint.col, y = cluster.col)) +
      geom_point(
        aes_string(size = "proportion", fill = population.col),
        alpha = 0.75,
        shape = 21
      ) +
      geom_line(aes(x = timepoint, y = cluster, group = group),
                edges_chart,
                arrow = arrow(length = unit(as.numeric(arrow.length), "mm"))) +
      scale_size_continuous(limit = c(min(dat_chart$proportion), max(dat_chart$proportion)),
                            range = c(min.node.size, max.node.size)) +
      labs(x = "",
           y = "",
           size = "Proportion of cells",
           fill = population.col)  +
      theme(
        legend.key = element_blank(),
        axis.text.x = element_text(
          colour = "black",
          size = font.size,
          angle = 90,
          vjust = 0.3,
          hjust = 1
        ),
        axis.text.y = element_text(colour = "black", size = font.size),
        legend.text = element_text(size = font.size, colour = "black"),
        legend.title = element_text(size = font.size),
        panel.background = element_blank(),
        panel.border = element_rect(
          colour = "black",
          fill = NA,
          size = 1.2
        ),
        legend.position = legend.position,
        plot.title = element_text(
          color = "Black",
          face = "bold",
          size = font.size,
          hjust = 0
        ),
        plot.margin = grid::unit(c(0, 0, 0, 0), "mm")
      ) +
      scale_fill_manual(values = c(col.scheme(num_populations))) +
      guides(fill = guide_legend(override.aes = list(size = min.node.size))) +
      scale_y_discrete(limits = rev(levels(dat_chart[[cluster.col]])))

    if (!(is.null(plot.title))) {
      plt <- plt + ggtitle(plot.title)
    }

    if (map.cluster.id) {
      filename <- paste0("Timeseries_heatmap_mapped_by_", population.col, ".", file.format)
    } else {
      filename <- paste0("Timeseries_heatmap_by_", population.col, ".", file.format)
    }

    ggsave(filename = filename,
           plot = plt,
           width = plot.width,
           height = plot.height,
           limitsize = FALSE)
  }

}
