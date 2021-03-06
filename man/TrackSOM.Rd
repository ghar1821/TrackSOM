% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TrackSOM.R
\name{TrackSOM}
\alias{TrackSOM}
\title{TrackSOM}
\usage{
TrackSOM (inputFiles, colsToUse, maxMeta, ...)
}
\arguments{
\item{inputFiles}{NO DEFAULT. Vector of filenames, each corresponds to the
file containing data for a time point}

\item{colsToUse}{NO DEFAULT. Vector of column names to be used to train SOM.
One list which must be same for all time periods.}

\item{tracking}{DEFAULT TRUE. Whether to track cluster changes (TRUE), or
just cluster (NULL).}

\item{noMerge}{DEFAULT FALSE. Whether to allow meta-cluster merging (FALSE) 
or not (TRUE).}

\item{maxMeta}{DEFAULT NULL. Numeric. Maximum number of meta clusters for all
time point. TrackSOM use this to limit the number of meta clusters allowed
when it is given the freedom to decide the optimum number of meta clusters.}

\item{nClus}{DEFAULT NULL. Numeric or Vector of numbers. If single number,
TrackSOM will produce same number of meta clusters per time point.
Otherwise, it will use the entries in the vector. Number of meta clusters
must be at least 3 for each time point.}

\item{seed}{DEFAULT 42. Random seed number.}

\item{xdim}{DEFAULT 10. SOM grid size.}

\item{ydim}{DEFAULT 10. SOM grid size.}

\item{dataFileType}{DEFAULT data.frame. What is the type of your inputFiles?}

\item{...}{other parameters for FlowSOM functions ReadInput, BuildSOM and
BuildMST. See FlowSOM vignette for specific parameter information}
}
\value{
A list containing the following elements:
1) FlowSOM: a list which is very similar to output of FlowSOM function.
2) Metaclustering: a list of two lists:
  2.1) Metaclustering: contains the clustering levels for each node in grid
  for each time step.
  2.2) mcPerCell: contains the clusering levels for each datapoint for each
  time step.
3) Tracking: a list of containing the unique IDs as determined by the
tracking function for each time step
}
\description{
Run the TrackSOM algorithm.
Method to run general TrackSOM workflow.
Will run Tracking by default with no merging allowed.
}
\examples{

data_files <- sapply(c(0:4), function(i) {
system.file("extdata", paste0("synthetic_d", i, ".fcs"), package="TrackSOM")
})
use_cols <- c("x", "y", "z")
tracksom_result <- TrackSOM(inputFiles = data_files,
                            colsToUse = use_cols,
                            nClus = 10,
                            dataFileType = ".fcs"
)
tracksom_res <- ExportClusteringDetailsOnly(tracksom_result)

}
