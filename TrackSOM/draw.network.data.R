library(Spectre)
library(gtools)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(RColorBrewer)
library(viridis)

dirname(rstudioapi::getActiveDocumentContext()$path)
PrimaryDirectory <- '/Users/givanna/Documents/phd/covid_data/immune_panel/all_patients'

setwd(PrimaryDirectory)
list.files()
dat.list <- Spectre::read.files()

# try making network from 2 time points
dat <- rbindlist(dat.list, fill=TRUE)
timepoints <- unique(dat$FileName)

get.idx.roundclsbracket <- function(cl) {
  tail(gregexpr("\\)", cl)[[1]], n=1)
}


# get all the transitions first
edge.df <- data.frame(from=character(),
                      to=character())


for (tp.idx in c(2:length(timepoints))) {
  prev.tp.dat <- dat[FileName == timepoints[tp.idx-1], ]
  prev.tp.clust <- mixedsort(unique(prev.tp.dat$TrackSOM_metacluster_lineage_tracking))
  
  complex.clusters <- prev.tp.clust[lapply(prev.tp.clust, get.idx.roundclsbracket) > -1]
  simple.clusters <- setdiff(prev.tp.clust, complex.clusters)
  
  # Order the vector based on the length of each element. 
  ## Need to use rev as whatever is in it order it in order it in increasing manner. STUPID!
  simple.clusters <- rev(simple.clusters[order(sapply(simple.clusters,length))])
  
  curr.tp.dat <- dat[FileName == timepoints[tp.idx], ]
  curr.tp.clust <- mixedsort(unique(curr.tp.dat$TrackSOM_metacluster_lineage_tracking))
  
  for (cl in curr.tp.clust) {
    clean.cl <- cl
    ## Match the complex clusters first
    for (cls in complex.clusters) {
      cls.found <- grepl(cls, clean.cl, fixed = TRUE)
      if (cls.found) {
        df <- data.frame(paste0(tp.idx-1,'_', cls), paste0(tp.idx, '_', cl))
        names(df) <- c("from","to")
        edge.df <- rbind(edge.df, df)
        clean.cl <- gsub(cls, "", clean.cl)
      }
    }
    ## Then simple clusters
    for (cls in simple.clusters) {
      cls.found <- grepl(cls, clean.cl, fixed = TRUE)
      if (cls.found) {
        df <- data.frame(paste0(tp.idx-1,'_',cls), paste0(tp.idx,'_',cl))
        names(df) <- c("from","to")
        edge.df <- rbind(edge.df, df)
        clean.cl <- gsub(cls, "", clean.cl)
      }
    }
  }
}

# then get extra details on the nodes
all.clust <- lapply(c(1:length(timepoints)), function(tp.idx) {
  tp.dat <- dat[FileName == timepoints[tp.idx], ]
  tp.clust <- mixedsort(unique(tp.dat$TrackSOM_metacluster_lineage_tracking))
  tp.clust.name <- sapply(tp.clust, function(c) {
    return(paste0(tp.idx, '_', c))
  })
  return(tp.clust.name)
})
all.clust <- unlist(all.clust)

cluster.ids <- sapply(all.clust, function(c) {
  strsplit(c, '_')[[1]][2]
})
timepoints.col <- sapply(all.clust, function(c) {
  strsplit(c, '_')[[1]][1]
})
node.dat <- data.frame(nodeId = all.clust,
                       clusterId = cluster.ids,
                       timepoints = timepoints.col)
node.dat$id <- seq.int(nrow(node.dat))
# node.dat <- node.dat[, c(4,2, 1, 3)]

# convert edges so it's to and from the numeric id of the node
names(edge.df) <- c('source', 'destination')
edges <- edge.df %>%
  left_join(node.dat, by = c("source" = "nodeId")) %>%
  rename(from = id)
edges <- edges %>%
  left_join(node.dat, by = c("destination" = "nodeId")) %>%
  rename(to = id)
# filter out other node information
keep.cols <- c('from', 'to')
edges <- data.table(edges)
edges <- edges[,..keep.cols]

node.ids <- as.vector(node.dat$nodeId)
# cluster.mapping <- sapply(node.ids, function(nodeId) {
#   nodeId.split <- strsplit(nodeId, '_')[[1]]
#   filename <- timepoints[as.numeric(nodeId.split[1])]  ## this indicate the filename of the cell in dat
#   cluster.id <- nodeId.split[2]
# 
#   # filter the data out
#   subset.dat <- dat[FileName == filename & TrackSOM_metacluster_lineage_tracking == cluster.id, ]
#   pop.names <- subset.dat$Population
#   cnt.pop.names <- table(pop.names) / sum(table(pop.names))
#   most.frequent.pop.cnt <- max(cnt.pop.names)
# 
#   # if (most.frequent.pop.cnt < 0.5) {
#   #   print(paste0("Cluster ", cl, ' proportion: ', most.frequent.pop.cnt))
#   # }
# 
#   most.frequent.pop <- names(which(cnt.pop.names == most.frequent.pop.cnt))
# 
#   if (length(most.frequent.pop) > 1) {
#     print(paste0("Tie for cluster ",
#                  nodeId,
#                  ". Mapping: ",
#                  paste(most.frequent.pop, collapse=', '), ". Picking first one!"))
#   }
# 
#   return(most.frequent.pop[1])
# })
# 
# 
# node.dat$population <- unlist(cluster.mapping, use.names=FALSE)

# marker.index <- c(60:75)  # granzyme panel
marker.index <- c(61:74)  # immune panel

### Compute each cluster's mean expression
for (idx in marker.index) {
  marker <- names(dat)[idx]
  marker.mean <- sapply(as.vector(node.dat$nodeId), function(id) {
    id_split <- strsplit(id, '_')[[1]]
    timepoint <- timepoints[as.numeric(id_split[1])]
    cluster_id <- id_split[2]
    sub.dat <- dat[FileName == timepoint & TrackSOM_metacluster_lineage_tracking == cluster_id,]
    mean(sub.dat[[marker]])
  })
  node.dat[[marker]] <- marker.mean
}

### Compute number of cells in each cluster
cell.count <- sapply(as.vector(node.dat$nodeId), function(nId) {
  id_split <- strsplit(nId, '_')[[1]]
  timepoint <- timepoints[as.numeric(id_split[1])]
  cluster_id <- id_split[2]
  sub.dat <- dat[FileName == timepoint & 
                   TrackSOM_metacluster_lineage_tracking == cluster_id,]
  nrow(sub.dat)
})
node.dat$cellCount <- cell.count

### Add an extra column which will colour just the node using the very first cluster
cluster.tp1 <- node.dat[node.dat$timepoints == 1,]$clusterId
cluster.colours <- sapply(c(1:nrow(node.dat)), function(i) {
  row <- node.dat[i,]
  cluster.remain <- row$clusterId %in% cluster.tp1
  
  if (cluster.remain) {
    return(as.character(row$clusterId))
  } else {
    return("NA")
  }
})
node.dat$origin <- cluster.colours

### Start plotting

# routes_tidy <- igraph::graph_from_data_frame(edges,
#                                              vertices = node.dat,
#                                              directed = TRUE) %>% as_tbl_graph()

routes_tidy <- tbl_graph(nodes = as.data.frame(node.dat),
                         edges = as.data.frame(edges),
                         directed = TRUE)
node.size <- 7

setwd(PrimaryDirectory)
dir.create("inferno")
setwd("inferno")


ggraph(routes_tidy, layout = 'kk', maxiter = 10000) +
  geom_edge_link(arrow = arrow(length = unit(1, 'mm')), end_cap = circle(2, 'mm')) +
  geom_node_point(aes(colour = timepoints), size=node.size) +
  scale_colour_manual(values = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(10)),
                      breaks = as.character(unique(node.dat$timepoints))) +
  theme(text = element_text(size=20),
        panel.background = element_rect(fill = 'white'),
        aspect.ratio = 1) + 
  labs(color='DPO bin')

ggsave(paste0("network_timepoints.pdf"),
       width = 15,
       height = 15,
       dpi = 1000,
       limitsize = FALSE)

## Colour plot based on the origin cluster
## +1 for the colour brewer because we have NA
ggraph(routes_tidy, layout = 'kk', maxiter = 10000) +
  geom_edge_link(arrow = arrow(length = unit(1, 'mm')), end_cap = circle(2, 'mm')) +
  geom_node_point(aes(colour = origin), size=node.size) +
  scale_colour_manual(values = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(length(cluster.tp1) + 1)),
                      breaks = as.character(unique(node.dat$origin))) +
  theme(text = element_text(size=20),
        panel.background = element_rect(fill = 'white'),
        aspect.ratio = 1) + 
  labs(color='Cluster origin')

ggsave(paste0("network_origin.pdf"),
       width = 15,
       height = 15,
       dpi = 1000,
       limitsize = FALSE)

# ggraph(routes_tidy, layout = 'sugiyama', maxiter=1000, hgap=4, vgap=0.5) +
# ggraph(routes_tidy, layout = 'kk', maxiter = 10000) +
#   geom_edge_link(arrow = arrow(length = unit(1.5, 'mm')), end_cap = circle(2, 'mm')) +
#   geom_node_point(aes(size = cellCount, colour = cellCount)) +
#   scale_colour_viridis_c(option='inferno') +
#   theme(text = element_text(size=20),
#         panel.background = element_rect(fill = 'white'))
# ggsave('network_cellCount.pdf',
#        width = 20,
#        height = 10,
#        dpi = 1000,
#        limitsize = FALSE)


# routes_tidy <- igraph::graph_from_data_frame(edges,
#                                              vertices = node.dat,
#                                              directed = TRUE) %>% as_tbl_graph()

for (idx in marker.index) {
  marker <- names(dat)[idx]
  
  # ggraph(routes_tidy, layout = 'sugiyama', maxiter = 1000, hgap=4, vgap=0.5) +
  ggraph(routes_tidy, layout = 'kk', maxiter = 10000) +
    geom_edge_link(arrow = arrow(length = unit(1, 'mm')), end_cap = circle(2, 'mm')) +
    geom_node_point(aes_string(colour = marker), size=node.size) +
    scale_colour_viridis_c(option='inferno') +
    # scale_colour_gradientn(colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(50))) +
    labs(color=stringr::str_split(marker, "_")[[1]][1]) +
    theme(text = element_text(size=20),
          panel.background = element_rect(fill = 'white'))
  
  ggsave(paste0(marker, '.pdf'),
         width = 15,
         height = 15,
         dpi = 1000,
         limitsize = FALSE)
}


##################################################################
# Plot 1 chart per population annotated by Tom
##################################################################

# node.dat <- data.table(node.dat)
# edges <- data.table(edges)
# 
# setwd(PrimaryDirectory)
# dir.create("plots_per_population")
# 
# populations <- unique(dat$Population)
# 
# for (population in populations) {
#   subset.dat <- data.table(dat[Population == population,])
#   exists.nodeId <- apply(subset.dat, 1, function(x) {
#     timepoint <- which(timepoints == x[['FileName']])
#     clusterId <- x[['TrackSOM_metacluster_lineage_tracking']]
#     return(paste0(timepoint, '_', clusterId))
#   })
#   
#   nodes.subset <- node.dat[nodeId %in% exists.nodeId, c(1:4)]
#   edges.subset <- edges[from %in% unique(nodes.subset$id) & to %in% unique(nodes.subset$id),]
#   
#   ## plot colour by time point first
#   routes_tidy <- igraph::graph_from_data_frame(edges.subset,
#                                                vertices = nodes.subset,
#                                                directed = TRUE) %>% as_tbl_graph()
#   
#   setwd(PrimaryDirectory)
#   setwd("plots_per_population")
#   save.dir <- population
#   dir.create(save.dir)
#   setwd(save.dir)
#   
#   ggraph(routes_tidy, layout = 'sugiyama', maxiter = 1000, hgap=4, vgap=0.5) +
#     geom_edge_link(arrow = arrow(length = unit(1.5, 'mm')), end_cap = circle(2, 'mm')) +
#     geom_node_point(aes(colour = timepoints), size=4) +
#     theme_classic() + 
#     theme(text = element_text(size=20))
#   
#   ggsave("network_timepoints.pdf",
#          width = 20,
#          height = 10,
#          dpi = 1000,
#          limitsize = FALSE)
#   
#   for (idx in marker.index) {
#     marker <- names(dat)[idx]
#     marker.mean <- sapply(as.vector(nodes.subset$nodeId), function(id) {
#       id_split <- strsplit(id, '_')[[1]]
#       timepoint <- timepoints[as.numeric(id_split[1])]
#       cluster_id <- id_split[2]
#       sub.dat <- subset.dat[FileName == timepoint & TrackSOM_metacluster_lineage_tracking == cluster_id,]
#       mean(sub.dat[[marker]])
#     })
#     nodes.subset[[marker]] <- marker.mean
#   }
#   
#   
#   routes_tidy <- igraph::graph_from_data_frame(edges.subset,
#                                                vertices = nodes.subset,
#                                                directed = TRUE) %>% as_tbl_graph()
#   
#   for (idx in marker.index) {
#     marker <- names(dat)[idx]
#     
#     ggraph(routes_tidy, layout = 'sugiyama', maxiter = 1000, hgap=4, vgap=0.5) +
#       geom_edge_link(arrow = arrow(length = unit(1.5, 'mm')), end_cap = circle(2, 'mm')) +
#       geom_node_point(aes_string(colour = marker), size=4) +
#       scale_colour_viridis_c(option='inferno') +
#       labs(color="Mean") +
#       theme_classic() + 
#       theme(text = element_text(size=20))
#     
#     ggsave(paste0(marker, '.pdf'),
#            width = 20,
#            height = 10,
#            dpi = 1000,
#            limitsize = FALSE)
#   }
# }


