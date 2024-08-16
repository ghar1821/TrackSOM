#' S3 ParseInput function
#' 
#' Internal Function.
#' Convert the input data to FlowSOM object ready to be processed.
#' Many of the parameters are used by FlowSOM's ReadInput function.
#' The description for these parameters are taken from FlowSOM.
#' 
#' @param inputFiles The list containing the input data
#' @param compensate Whether to run data compensation
#' @param spillover The spillover matrix to compensate
#' @param transform Whether to run data transformation
#' @param toTransform column names or indices to transform
#' @param transformFunction The transform function to run
#' @param scale Whether to scale the data
#' @param scale.center See \link{base::scale}
#' @param scale.scale See \link{base::scale}
#' @param silent Whether to show progress updates
#' @param colToUse Which column to be used by TrackSOM 
#' @param ... Additional arguments used if inputFiles is a data.table object
#' 
#' @return FlowSOM object
#' 
ParseInput <- function(
        inputFiles,
        compensate,
        spillover,
        transform,
        toTransform,
        transformFunction,
        scale,
        scaled.center,
        scaled.scale,
        silent,
        colsToUse,
        ...
    ) {
    UseMethod("ParseInput")
}

#' @describeIn ParseInput Process vector of characters input
#' 
ParseInput.character <- function(
        inputFiles,
        compensate,
        spillover,
        transform,
        toTransform,
        transformFunction,
        scale,
        scaled.center,
        scaled.scale,
        silent,
        colsToUse,
        ...
    ) {
    
    # if it is a vector of filenames, it'll go to here.
    
    inputFilesType <- unique(tools::file_ext(inputFiles))
    
    if (length(inputFilesType) > 1) {
        stop(paste(
            "Multiple file types are detected (",
            paste(inputFilesType, collapse = ", "),
            "). Please pass only one file types in inputFiles!"
        ))
    }
    
    inputFilesType <- inputFilesType[1]
    
    if (! inputFilesType %in% c("fcs", "csv")) {
        stop(paste(
            "inputFiles are .", inputFilesType,
            "file type. Only fcs or csv files are supported."
        ))
    }
    
    if (inputFilesType == 'csv') {
        data_ff <- lapply(inputFiles, function(input_file) {
            input_dat <- fread(input_file)
            input_dat_ff <-create_flowframe(input_dat, colsToUse)
            return(input_dat_ff)
        })
        data_flowset <- as(data_ff, "flowSet")
        
        fsom <- ReadInput(
            data_flowset,
            compensate  =  compensate,
            spillover = spillover,
            transform = transform,
            toTransform = toTransform,
            transformFunction = transformFunction,
            scale = scale,
            scaled.center = scaled.center,
            scaled.scale = scaled.scale,
            silent = silent
        )
    } else {
        # will only get here if data files are fcs file
        fsom <- ReadInput(
            inputFiles,
            pattern = ".fcs",
            compensate  =  compensate,
            spillover = spillover,
            transform = transform,
            toTransform = toTransform,
            transformFunction = transformFunction,
            scale = scale,
            scaled.center = scaled.center,
            scaled.scale = scaled.scale,
            silent = silent
        )
    }
    
    return(fsom)
    
}

#' @describeIn ParseInput Process list of data.frames input
#' @import data.table
#' 
ParseInput.list <- function(
        inputFiles,
        compensate,
        spillover,
        transform,
        toTransform,
        transformFunction,
        scale,
        scaled.center,
        scaled.scale,
        silent,
        colsToUse,
        ...
    ) {
    
    # preliminary error checking
    for (i in length(inputFiles)) {
        if (! is.data.table(inputFiles[[i]])) {
            stop(paste(
                "element at position", i, 
                "is not a data.table Please check."
            ))
        }
    }
    
    data_ff <- lapply(inputFiles, function(input_file) {
        input_dat_ff <- create_flowframe(data.table(input_file), colsToUse)
        return(input_dat_ff)
    })
    data_flowset <- as(data_ff, "flowSet")
    
    fsom <- ReadInput(
        data_flowset,
        compensate  =  compensate,
        spillover = spillover,
        transform = transform,
        toTransform = toTransform,
        transformFunction = transformFunction,
        scale = scale,
        scaled.center = scaled.center,
        scaled.scale = scaled.scale,
        silent = silent
    )
    
    return(fsom)
}

#' @describeIn ParseInput Process data.table input
#' @param timepoints A vector of time-points.
#' @param timepointCol A character denoting the column in data.frame that
#' corresponds to the time-point each cell belongs to
#' 
#' @import data.table
#' 
ParseInput.data.table <- function(
        inputFiles,
        compensate,
        spillover,
        transform,
        toTransform,
        transformFunction,
        scale,
        scaled.center,
        scaled.scale,
        silent,
        colsToUse,
        timepointCol,
        timepoints
) {
    if (! timepointCol %in% names(inputFiles)) {
        stop(paste(timepointCol, "column doesn't exist in inputFiles!"))
    }
    
    inputFiles <- data.table(inputFiles, keep.rownames = TRUE)
    
    data_ff <- lapply(timepoints, function(t) {
        input_dat <- inputFiles[get(timepointCol) == eval(t)]
        
        input_dat_ff <- create_flowframe(input_dat, colsToUse)
        return(input_dat_ff)
    })
    data_flowset <- as(data_ff, "flowSet")
    
    fsom <- ReadInput(
        data_flowset,
        compensate  =  compensate,
        spillover = spillover,
        transform = transform,
        toTransform = toTransform,
        transformFunction = transformFunction,
        scale = scale,
        scaled.center = scaled.center,
        scaled.scale = scaled.scale,
        silent = silent
    )
    
    return(fsom)
}

#' Create flowframe
#'
#' @param input_dat a data.table
#' @param colsToUse a vector of columns to isolate from input_dat
#' 
#' @import data.table
#' 
#' @return flowFrame object
#' 
create_flowframe <- function(input_dat, colsToUse) {
    input_dat <- input_dat[, colsToUse, with = FALSE]
    
    # convert to flowFrame. Pain in the butt
    metadata <-
        data.frame(
            name = dimnames(input_dat)[[2]],
            desc = paste('column', dimnames(input_dat)[[2]], 'from dataset')
        )
    # in order to create a flow frame, data (exprs param) needs to be read as
    # matrix
    input_dat_ff <- new(
        "flowFrame",
        exprs = as.matrix(input_dat),
        parameters = Biobase::AnnotatedDataFrame(metadata)
    )
    
    return(input_dat_ff)
}