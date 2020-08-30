library(Spectre)
library(gtools)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(RColorBrewer)
library(viridis)

dirname(rstudioapi::getActiveDocumentContext()$path)
PrimaryDirectory <- '~/Documents/phd/code/TrackSOM/covid19_paper/TrackSOM_result'

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

marker.index <- c(61:74)

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