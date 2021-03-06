% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TrackingFunctions.R
\name{GetNodesPerCluster}
\alias{GetNodesPerCluster}
\title{GetNodesPerCluster}
\usage{
GetNodesPerCluster(mc)
}
\arguments{
\item{mc}{A factor of the meta-cluster per node}
}
\value{
A list of node IDs for each metacluster ID
}
\description{
Function which generates a list of the nodes for each metacluster, used to
generate data matrix for tracking.
}
