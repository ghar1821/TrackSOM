context("Test reading input files")
library(TrackSOM)
library(data.table)

test_that("CSV files are parsed", {
    data_files <- sapply(c(0:4), function(i) {
        system.file("extdata", paste0("synthetic_d", i, ".csv"), package = "TrackSOM")
    })
    fsom <- ParseInput(
        inputFiles = data_files,
        colsToUse = c("x", "y", "z"),
        dataFileType = '.csv',
        compensate = FALSE,
        spillover = NULL,
        transform = FALSE,
        toTransform = NULL,
        transformFunction = flowCore::logicleTransform(),
        scale = FALSE,
        scaled.center = TRUE,
        scaled.scale = TRUE,
        silent = TRUE
    )
    expect_equal(nrow(fsom$data), 35700)
    expect_equal(ncol(fsom$data), 3)
})

test_that("FCS files are parsed", {
    data_files <- sapply(c(0:4), function(i) {
        system.file("extdata", paste0("synthetic_d", i, ".fcs"), package = "TrackSOM")
    })
    fsom <- ParseInput(
        inputFiles = data_files,
        colsToUse = c("x", "y", "z"),
        dataFileType = '.fcs',
        compensate = FALSE,
        spillover = NULL,
        transform = FALSE,
        toTransform = NULL,
        transformFunction = flowCore::logicleTransform(),
        scale = FALSE,
        scaled.center = TRUE,
        scaled.scale = TRUE,
        silent = TRUE
    )
    expect_equal(nrow(fsom$data), 35700)
    expect_equal(ncol(fsom$data), 3)
})

test_that("Data table is parsed", {
    data <- lapply(c(0:4), function(i) {
        dat_file <-
            system.file("extdata", paste0("synthetic_d", i, ".csv"), package = "TrackSOM")
        return(fread(dat_file))
    })
    fsom <- ParseInput(
        inputFiles = data,
        colsToUse = c("x", "y", "z"),
        dataFileType = 'data.frame',
        compensate = FALSE,
        spillover = NULL,
        transform = FALSE,
        toTransform = NULL,
        transformFunction = flowCore::logicleTransform(),
        scale = FALSE,
        scaled.center = TRUE,
        scaled.scale = TRUE,
        silent = TRUE
    )
    expect_equal(nrow(fsom$data), 35700)
    expect_equal(ncol(fsom$data), 3)
})
