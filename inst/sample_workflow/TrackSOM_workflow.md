TrackSOM workflow
================

In this tutorial, we outline a step by step instruction on how to use
TrackSOM to cluster and track time-series data. For simplicity, we will
use the synthetic dataset outlined in our manuscript downloadable from
[bioArxiv](https://www.biorxiv.org/content/10.1101/2021.06.08.447468v1).
The dataset is stored under `inst` directory:
[link](https://github.com/ghar1821/TrackSOM/tree/master/inst/extdata).

# Importing TrackSOM package

If you have not installed the TrackSOM package, please use `devtools` to
install the package from the following [TrackSOM github
repo](https://github.com/ghar1821/TrackSOM). For devtools, the repo
parameter will be: `ghar1821/TrackSOM`.

The following code shall import the TrackSOM package:

``` r
library(TrackSOM)
```

# Importing dataset

TrackSOM supports dataset stored as either CSV or FCS files or as
`data.table` object (enhanced version of R’s native data.frame. See
[data.table
vignette](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html)
for more details). The following sections shall show you how to read in
those files and pass them to the TrackSOM function.

## Dataset as CSV files

Here, we assume that each file contains data belonging to one
time-point. Hence you should have more than 1 CSV files. If this is not
the case, please reformat your data files.

To import the dataset stored as CSV files, TrackSOM needs to know where
the files are stored. These files’ location must be stored within a
vector which get passed on to the TrackSOM function.

*Important:* the vector must be organised such as the first element is
the data for the very first time-point, the 2nd element for the 2nd
time-point, and so on.

First, we start with specifying the CSV files are:

``` r
InputDirectory <- "~/Documents/GitHub/TrackSOM/inst/extdata/"
data.files <- list.files(InputDirectory, ".csv")
data.files.fullpath <- sapply(data.files, function(fname) {
    full.fname <- paste(InputDirectory, fname, sep="/")
})
```

Let’s inspect the content of `data.files.fullpath`:

``` r
print(data.files.fullpath)
```

    ##                                             synthetic_d0.csv 
    ## "~/Documents/GitHub/TrackSOM/inst/extdata//synthetic_d0.csv" 
    ##                                             synthetic_d1.csv 
    ## "~/Documents/GitHub/TrackSOM/inst/extdata//synthetic_d1.csv" 
    ##                                             synthetic_d2.csv 
    ## "~/Documents/GitHub/TrackSOM/inst/extdata//synthetic_d2.csv" 
    ##                                             synthetic_d3.csv 
    ## "~/Documents/GitHub/TrackSOM/inst/extdata//synthetic_d3.csv" 
    ##                                             synthetic_d4.csv 
    ## "~/Documents/GitHub/TrackSOM/inst/extdata//synthetic_d4.csv"

You can see it contains a list of *absolute path* of multiple CSV files,
each belonging to a time-point.

## Dataset as FCS files

The procedure is very similar to reading CSV files above, except we’re
going to be storing FCS files’ path rather than CSV.

Again, here, we assume that each file contains data belonging to one
time-point. Hence you should have more than 1 FCS files. If this is not
the case, please reformat your data files.

To import the dataset stored as FCS files, TrackSOM needs to know where
the files are stored. These files’ location must be stored within a
vector which get passed on to the TrackSOM function. *Important:* the
vector must be organised such as the first element is the data for the
very first time-point, the 2nd element for the 2nd time-point, and so
on.

First, we start with specifying the FCS files are:

``` r
InputDirectory <- "~/Documents/GitHub/TrackSOM/inst/extdata/"
data.files.fcs <- list.files(InputDirectory, ".fcs")
data.files.fullpath.fcs <- sapply(data.files.fcs, function(fname) {
    full.fname <- paste(InputDirectory, fname, sep="/")
})
```

Let’s inspect the content of `data.files.fullpath`:

``` r
print(data.files.fullpath.fcs)
```

    ##                                             synthetic_d0.fcs 
    ## "~/Documents/GitHub/TrackSOM/inst/extdata//synthetic_d0.fcs" 
    ##                                             synthetic_d1.fcs 
    ## "~/Documents/GitHub/TrackSOM/inst/extdata//synthetic_d1.fcs" 
    ##                                             synthetic_d2.fcs 
    ## "~/Documents/GitHub/TrackSOM/inst/extdata//synthetic_d2.fcs" 
    ##                                             synthetic_d3.fcs 
    ## "~/Documents/GitHub/TrackSOM/inst/extdata//synthetic_d3.fcs" 
    ##                                             synthetic_d4.fcs 
    ## "~/Documents/GitHub/TrackSOM/inst/extdata//synthetic_d4.fcs"

You can see it contains a list of *absolute path* of multiple FCS files,
each belonging to a time-point.

## Dataset as `data.table` object

Sometimes, it is convenient to have the dataset stored as the
`data.table` object files, e.g. when you need to run some code to
preprocess your data using R! As an example, supposed the synthetic
dataset CSV files were already read in as `data.table` object prior to
running TrackSOM (say you did some preliminary clean up or filtering).
What you need to do is organise them in a list such that each element is
a `data.table` object for the dataset in a time-point. *Important:* the
list must be organised such as the first element is the data for the
very first time-point, the 2nd element for the 2nd time-point, and so
on.

``` r
library(data.table)
InputDirectory <- "~/Documents/GitHub/TrackSOM/inst/extdata/"
data.files <- list.files(InputDirectory, ".csv")
dat <- lapply(data.files, function(f) fread(paste0(InputDirectory, f)))
```

Let’s do a quick preview of the data:

``` r
dat
```

    ## [[1]]
    ##               x         y         z timepoint
    ##    1:  9.598217 11.728198 11.382064      Mock
    ##    2:  8.937549  8.653750 11.107314      Mock
    ##    3:  9.632932  9.657278  9.115804      Mock
    ##    4:  9.770393 12.478038 11.451082      Mock
    ##    5: 30.426699 31.290929 28.901204      Mock
    ##   ---                                        
    ## 7096: 10.533174  9.462175  9.867354      Mock
    ## 7097: 10.544833 10.426697  9.571881      Mock
    ## 7098: 10.140938  8.727880 10.247086      Mock
    ## 7099: 30.967131 30.261241 30.469085      Mock
    ## 7100: 10.063984 10.289172 10.382876      Mock
    ## 
    ## [[2]]
    ##               x         y         z timepoint
    ##    1: 33.850857 30.653891 29.257884     SYN-1
    ##    2: 10.429818 13.880273 10.878985     SYN-1
    ##    3: 10.009247 10.327127 10.994237     SYN-1
    ##    4:  9.751578  8.555277  8.914872     SYN-1
    ##    5: 10.456213 10.972776 10.011441     SYN-1
    ##   ---                                        
    ## 7096:  9.237145 10.203747 12.475563     SYN-1
    ## 7097: 12.957000  9.147637  8.807634     SYN-1
    ## 7098: 12.297580 11.034141  9.828642     SYN-1
    ## 7099:  9.393480 11.766797 11.137156     SYN-1
    ## 7100: 11.684874  8.864774 11.771160     SYN-1
    ## 
    ## [[3]]
    ##               x         y         z timepoint
    ##    1: 10.107081 10.371958  8.862552     SYN-2
    ##    2: 12.271113 28.710774 13.945666     SYN-2
    ##    3: 10.248482 10.249500  7.918076     SYN-2
    ##    4:  9.341533 10.228882 12.209994     SYN-2
    ##    5: 10.116334  9.805188  9.742291     SYN-2
    ##   ---                                        
    ## 7196: 10.682697 13.643521 10.637965     SYN-2
    ## 7197: 27.598520 30.275747 29.691414     SYN-2
    ## 7198: 10.728449 10.020652  9.517445     SYN-2
    ## 7199: 10.167862  9.754739  9.416342     SYN-2
    ## 7200:  9.554067 11.296156 12.535182     SYN-2
    ## 
    ## [[4]]
    ##               x         y         z timepoint
    ##    1: 13.864150  8.884197  8.742413     SYN-3
    ##    2: 25.420014 40.961726 26.453685     SYN-3
    ##    3: 23.968729 31.313815 31.314857     SYN-3
    ##    4: 10.538730  9.855319  9.061674     SYN-3
    ##    5:  8.814686 12.068781 10.729973     SYN-3
    ##   ---                                        
    ## 7196: 10.358157  8.856995 10.440222     SYN-3
    ## 7197: 10.093926 11.400187 10.429873     SYN-3
    ## 7198: 14.720125 10.599644  9.599185     SYN-3
    ## 7199: 24.545254 29.734070 29.035897     SYN-3
    ## 7200:  9.967504 14.351193 11.619231     SYN-3
    ## 
    ## [[5]]
    ##               x         y         z timepoint
    ##    1: 10.088854 10.171541 10.122280     SYN-4
    ##    2: 14.776449 10.189092 10.621575     SYN-4
    ##    3:  9.652116 12.567544  9.282672     SYN-4
    ##    4: 10.436809 10.725664 10.016877     SYN-4
    ##    5: 10.235657 10.966951  9.371363     SYN-4
    ##   ---                                        
    ## 7096: 21.234807 30.522073 30.590687     SYN-4
    ## 7097:  9.703690  8.319128  9.579323     SYN-4
    ## 7098: 17.654003 30.891041 29.390515     SYN-4
    ## 7099:  9.089960  9.744563 12.483508     SYN-4
    ## 7100: 36.910728 29.551236 29.926378     SYN-4

As you can see, there are 5 elements in the list, each containing a
dataset belonging to a time-point.

# Running TrackSOM

Depending on how your dataset is stored, the parameter `inputFiles` is
either a vector of absolute path of your CSV or FCS files or a list of
`data.table` object. Additionally, you need to specify the type as the
parameter `dataFileType`. It can be either `.csv`, `.fcs`, or
`data.frame` depending on how your dataset is stored.

## Other parameters

***Note:*** examples here assume datasets are stored as CSV files.

The TrackSOM function have various parameters. Of them, only
`inputFiles` and `colsToUse` have no default values. The `inputFiles`
parameter has been explained in the previous section.

The `colsToUse` parameter specify the columns in your dataset to be used
for clustering and tracking. This must be a vector. For the synthetic
dataset, the columns are denoted as `x`, `y`, and `z`. For cytometry
data, this should be a vector of markers.

The TrackSOM functions have the following parameters which are
pre-filled with default values:

-   tracking: TRUE (default) or FALSE, whether to perform tracking.
-   noMerge: TRUE or FALSE (default), whether to allow merging of
    meta-clusters (FALSE) or not (TRUE). If noMerge is set to TRUE,
    meta-clusters are only allowed to split, but not merge.
-   maxMeta: NULL (default) or round number. Maximum number of meta
    clusters for all time point. Only used for Autonomous Adaptive
    operation.
-   nClus: NULL (default) or numeric or a vector of numbers. If single
    number, TrackSOM will run Prescribed Invariant operation. Otherwise,
    it will run Prescribed Variant operation.
-   seed: 42 (default). Random seed number.
-   xdim: 10 (default). SOM grid size.
-   ydim: 10 (default). SOM grid size.

The function also accept parameters that are built into
[FlowSOM](https://bioconductor.org/packages/release/bioc/html/FlowSOM.html)
`ReadInput`, `BuildSOM` and `BuildMST` functions. See FlowSOM’s vignette
for specific parameter information.

Let’s run TrackSOM with the following settings:

-   No merging of meta-clusters allowed (`noMerge = TRUE`)
-   Enable tracking (`tracking = TRUE`)
-   Prescribed Variant producing 3, 3, 9, 7, 15 meta-clusters for each
    time-point.

The remaining parameters will be set to the default values.

``` r
tracksom.result <- TrackSOM(inputFiles = data.files.fullpath,
                            colsToUse = c('x', 'y', 'z'),
                            tracking = TRUE,
                            noMerge = TRUE,
                            nClus = c(3,3,9,7,15),
                            dataFileType = ".csv"
)
```

    ## Building SOM

    ## Mapping data to SOM

    ## Building MST

    ## Extracting SOM nodes for each time point

    ## Running meta clustering

    ## Meta clustering time point 1 with 82 SOM nodes

    ## Meta clustering time point 2 with 93 SOM nodes

    ## Meta clustering time point 3 with 90 SOM nodes

    ## Meta clustering time point 4 with 95 SOM nodes

    ## Meta clustering time point 5 with 95 SOM nodes

# Extracting meta-clusters ID and SOM nodes

TrackSOM result is stored in an object. To facilitate the extraction of
meta-clusters ID and SOM nodes for each cell, we provide the following
functions:

-   `ConcatenateClusteringDetails`: assuming your dataset is stored as a
    `data.table` object, the function attaches the meta-cluster ID and
    SOM nodes as separate columns.
-   `ExportClusteringDetailsOnly`: this function simply extract the
    meta-cluster ID and SOM nodes for each cell as a `data.table`
    object. The cells are ordered based on the ordering of your dataset.

To use `ConcatenateClusteringDetails`, you need to first read in all
your datasets as one giant `data.table`.

``` r
library(data.table)
InputDirectory <- "~/Documents/GitHub/TrackSOM/inst/extdata/"
data.files <- list.files(InputDirectory, ".csv")
dat <- lapply(data.files, function(f) fread(paste0(InputDirectory, f)))
dat <- rbindlist(dat)
```

Double check the data is read in properly:

``` r
head(dat)
```

    ##            x         y         z timepoint
    ## 1:  9.598217 11.728198 11.382064      Mock
    ## 2:  8.937549  8.653750 11.107314      Mock
    ## 3:  9.632932  9.657278  9.115804      Mock
    ## 4:  9.770393 12.478038 11.451082      Mock
    ## 5: 30.426699 31.290929 28.901204      Mock
    ## 6: 30.789497 30.069909 30.767145      Mock

``` r
tail(dat)
```

    ##            x         y         z timepoint
    ## 1:  7.954781  9.241713  9.429919     SYN-4
    ## 2: 21.234807 30.522073 30.590687     SYN-4
    ## 3:  9.703690  8.319128  9.579323     SYN-4
    ## 4: 17.654003 30.891041 29.390515     SYN-4
    ## 5:  9.089960  9.744563 12.483508     SYN-4
    ## 6: 36.910728 29.551236 29.926378     SYN-4

To run `ConcatenateClusteringDetails`, you need to pass the following
parameters:

-   `timepoint.col`: which column in your data define the time-point.
-   `timepoints`: what are the time-points (in order). This must be a
    vector.

Now let’s attach the meta-cluster ID and the SOM nodes:

``` r
dat.clust <- ConcatenateClusteringDetails(
    tracksom.result = tracksom.result,
    dat = dat,
    timepoint.col = "timepoint",
    timepoints = c('Mock', 'SYN-1', 'SYN-2', 'SYN-3', 'SYN-4')
    )
```

Inspect the content:

``` r
head(dat.clust)
```

    ##            x         y         z timepoint TrackSOM_cluster
    ## 1:  9.598217 11.728198 11.382064      Mock               62
    ## 2:  8.937549  8.653750 11.107314      Mock               53
    ## 3:  9.632932  9.657278  9.115804      Mock               86
    ## 4:  9.770393 12.478038 11.451082      Mock               61
    ## 5: 30.426699 31.290929 28.901204      Mock               16
    ## 6: 30.789497 30.069909 30.767145      Mock                6
    ##    TrackSOM_metacluster TrackSOM_metacluster_lineage_tracking
    ## 1:                    2                                     B
    ## 2:                    2                                     B
    ## 3:                    2                                     B
    ## 4:                    2                                     B
    ## 5:                    1                                     A
    ## 6:                    1                                     A

The function attaches extra 3 columns:

-   `TrackSOM_cluster`: The SOM node of each cell.
-   `TrackSOM_metacluster`: The meta-cluster assignment produced by
    FlowSOM’s meta-clustering.
-   `TrackSOM_metacluster_lineage_tracking`: The tracking of
    meta-clusters’ evolution produced by TrackSOM. **This gives you the
    changes undergone by the meta-clusters over time**.

The function gives back a `data.table` object which you can export as
CSV file using `data.table`’s `fwrite` function.

# Visualisations

TrackSOM offers 2 visualisation mediums:

1.  Network plot.
2.  Timeseries heatmap.

## Network plot

To draw network plots, we need to call the `DrawNetworkPlot` function.
Using the clustered and tracked data from previous sections, the
following is an example on how to use the function:

``` r
DrawNetworkPlot(dat = dat.clust,
                timepoint.col = "timepoint",
                timepoints = c('Mock', 'SYN-1', 'SYN-2', 'SYN-3', 'SYN-4'),
                cluster.col = 'TrackSOM_metacluster_lineage_tracking',
                marker.cols = c('x', 'y', 'z'))
```

    ## Calculating edges

    ## Computing node details

    ## Calculating marker's average per node

    ## Saving node and edge details

    ## Start drawing plots

    ## Warning: Existing variables `x`, `y` overwritten by layout variables

    ## Drawing plots coloured by time point

    ## Drawing plots coloured by origin

    ## Drawing plots coloured by x

    ## Drawing plots coloured by y

    ## Drawing plots coloured by z

The `marker.cols` will determine the *markers* which mean/median
expression will be drawn on the network plots. There is no need to
specify all the markers in the dataset, just the ones that you want the
network plots to be coloured on.

The function won’t preview any plots, but it will instead store the
plots as image files and some extra information (median/mean expression
of markers) as CSV files:

``` r
list.files()
```

    ##  [1] "network_colBy_origin.pdf"      "network_colBy_timepoints.pdf" 
    ##  [3] "network_colBy_x.pdf"           "network_colBy_y.pdf"          
    ##  [5] "network_colBy_z.pdf"           "network_plot_edge_details.csv"
    ##  [7] "network_plot_node_details.csv" "TrackSOM_workflow.md"         
    ##  [9] "TrackSOM_workflow.nb.html"     "TrackSOM_workflow.R"          
    ## [11] "TrackSOM_workflow.Rmd"

In this example, the plots are saved as PDF files. This can be changed,
e.g. to save the plots as PNG files, by specifying the desired file
format as the `file.format` parameter.

## Timeseries heatmap

To draw a timeseries heatmap, we need to call the
`DrawTimeseriesHeatmap` function. Using the clustered and tracked data
from previous sections, the following is an example on how to use the
function:

``` r
DrawTimeseriesHeatmap(dat = dat.clust,
                timepoint.col = "timepoint",
                timepoints = c('Mock', 'SYN-1', 'SYN-2', 'SYN-3', 'SYN-4'),
                cluster.col = 'TrackSOM_metacluster_lineage_tracking',
                marker.cols = c('x', 'y', 'z'))
```

    ## Computing node details

    ## Computing edge details

    ## Saving node and edge details

    ## Drawing timeseries heatmap coloured by x

    ## Drawing timeseries heatmap coloured by y

    ## Drawing timeseries heatmap coloured by z

The `marker.cols` will determine the *markers* which mean/median
expression will be drawn on the network plots. There is no need to
specify all the markers in the dataset, just the ones that you want the
network plots to be coloured on.

The function won’t preview any plots, but it will instead store the
plots as image files and some extra information (median/mean expression
of markers) as CSV files:

``` r
list.files()
```

    ##  [1] "network_colBy_origin.pdf"            
    ##  [2] "network_colBy_timepoints.pdf"        
    ##  [3] "network_colBy_x.pdf"                 
    ##  [4] "network_colBy_y.pdf"                 
    ##  [5] "network_colBy_z.pdf"                 
    ##  [6] "network_plot_edge_details.csv"       
    ##  [7] "network_plot_node_details.csv"       
    ##  [8] "Timeseries_heatmap_by_x.pdf"         
    ##  [9] "Timeseries_heatmap_by_y.pdf"         
    ## [10] "Timeseries_heatmap_by_z.pdf"         
    ## [11] "timeseries_heatmap_edges_details.csv"
    ## [12] "timeseries_heatmap_node_details.csv" 
    ## [13] "TrackSOM_workflow.md"                
    ## [14] "TrackSOM_workflow.nb.html"           
    ## [15] "TrackSOM_workflow.R"                 
    ## [16] "TrackSOM_workflow.Rmd"

In this example, the plots are saved as PDF files. This can be changed,
e.g. to save the plots as PNG files, by specifying the desired file
format as the `file.format` parameter.

# Concluding remarks

That is pretty much it folks! We’re actively updating TrackSOM and
welcome feedbacks!

Thank you for your interest!
