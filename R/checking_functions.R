#' @importFrom methods is

## ideogram data
.checkIdeogramData <- function(ideogram_data)
{
    ## type check and convert to data frame
    if (methods::is(ideogram_data, "GRanges")) {
        ideogram_data <- .GRangesToBed(ideogram_data)

    } else if (!methods::is(ideogram_data, "data.frame")) {
        stop("ideogram_data should be either a data frame or GRanges")
    }

    ## check that the columns and the values are ok:
    stopifnot(exprs={
        ## have all of the columns we need
        all(c("chr", "start", "end") %in% colnames(ideogram_data))

        ## columns are the right types
        !any(is.na(as.character(ideogram_data$chr)))
        !any(is.na(as.numeric(ideogram_data$start)))
        !any(is.na(as.numeric(ideogram_data$end)))

        ## start and end values make sense
        all(ideogram_data$end - ideogram_data$start >= 0)

        ## we don't have any duplicates in sector names
        length(ideogram_data$chr) == length(unique(ideogram_data$chr))
    })
    return(ideogram_data)
}

## label data
.checkLabelData <- function(label_data)
{
    ## type check and convert to data frame
    if (methods::is(label_data, "GRanges")) {
        label_data <- .GRangesToBed(label_data)
    } else if (!methods::is(label_data, "data.frame")) {
        stop("label_data should be either a data frame or GRanges")
    }

    ## check that the columns and the values are ok:
    stopifnot(exprs={
        ## have all of the columns we need
        all(c("chr", "start", "end", "label") %in% colnames(label_data))

        ## columns are the right type
        !any(is.na(as.character(label_data$chr)))
        !any(is.na(as.numeric(label_data$start)))
        !any(is.na(as.numeric(label_data$end)))
        !any(is.na(as.character(label_data$label)))

        ## start and end values make sense
        all(label_data$end - label_data$start >= 0)

    })

    ## warn people if the colour column is factor: the colours get converted to
    ## numbers and they'll come out wrong
    if ("colour" %in% colnames(label_data) &
        methods::is(label_data$colour, "factor")) {
        warning("The 'colour' column is a factor. Colours may not appear as you
                intended - if this is the case please change it to character")
    }
    return(label_data)
}

## feature data
.checkFeatureData <- function(feature_data)
{
    ## type check and covert to data frame
    if (methods::is(feature_data, "GRanges")) {
        feature_data <- .GRangesToBed(feature_data)

    } else if (!methods::is(feature_data, "data.frame")) {
        stop("feature_data should be either a data frame or GRanges")
    }

    ## check for necessary columns
    stopifnot(all(c("chr", "start", "end") %in% colnames(feature_data)))

    ## check for the optional columns (label, colour, shape & track) and set
    ## default values if they haven't been supplied

    ## default label: none
    if (!"label" %in% colnames(feature_data)) {
        feature_data$label <- rep("", nrow(feature_data))
    }

    ## default shape: rectangle
    if (!"shape" %in% colnames(feature_data)) {
        feature_data$shape <- rep("rectangle", nrow(feature_data))
    }

    ## default track: 1
    if (!"track" %in% colnames(feature_data)) {
        feature_data$track <- rep(1, nrow(feature_data))
    }

    ## default colour: rich_colours, one colour per feature
    if (!"colour" %in% colnames(feature_data)) {
        feature_data$colour <- rich_colours[seq_len(nrow(feature_data))]
    }

    ## check that the columns and the values are ok
    stopifnot(exprs={
        ## columns are the right types
        !any(is.na(as.character(feature_data$chr)))
        !any(is.na(as.numeric(feature_data$start)))
        !any(is.na(as.numeric(feature_data$end)))
        !any(is.na(as.character(feature_data$label)))
        !any(is.na(as.character(feature_data$shape)))
        !any(is.na(as.numeric(feature_data$track)))
        !any(is.na(as.character(feature_data$colour)))

        ## values make sense
        all(feature_data$end - feature_data$start >= 0)
        all(feature_data$track >= 0)
        all(feature_data$shape %in% c(
            "rectangle", "forward_arrow", "reverse_arrow", "upwards_triangle",
            "downwards_triangle"))

    })

    ## warn people about the factor colour problem
    if (methods::is(feature_data$colour, "factor")) {
        warning("The 'colour' column is a factor. Colours may not appear as
                you intended - if this is the case please change it to
                character")
    }

    ## currently only rectangles are supported on the outermost track
    ## give a warning if people try to use another shape here
    for (i in seq_along(feature_data$track)) {
        if (feature_data$track[[i]] == 0 &
            feature_data$shape[[i]] != "rectangle") {
            warning("Currently only rectangles ('rectangle') are supported on
                    the outermost track (track 0). Your feature will still be
                    plotted, but as a rectangle.")
            break()
        }
    }
    return(feature_data)
}

## insertion data
.checkInsertionData <- function(insertion_data, multiple = FALSE) {

    ## type check and convert to data frame
    if (methods::is(insertion_data, "GRanges")) {
        insertion_data <- .GRangesToBed(insertion_data)

    } else if (!methods::is(insertion_data, "data.frame")) {
        stop("insertion_data should be either a data frame or GRanges")
    }

    ## check for necessary columns
    stopifnot(all(c("chr", "start", "end", "name", "length") %in%
                    colnames(insertion_data)))

    ## check for optional columns (shape, in_tandem, colour) and set defaults
    ## default shape: forward_arrow
    if (!"shape" %in% colnames(insertion_data)) {
        insertion_data$shape <- rep("forward_arrow", nrow(insertion_data))
    }

    ## default in_tandem: 1
    if (!"in_tandem" %in% colnames(insertion_data)) {
        insertion_data$in_tandem <- rep(1, nrow(insertion_data))
    }

    ## default colour: rich_colours, one colour per feature
    if (!"colour" %in% colnames(insertion_data)) {
        insertion_data$colour <- rich_colours[seq_len(nrow(insertion_data))]
    }

    ## check that the columns and the values are ok
    stopifnot(exprs={
        ## columns are the right types
        !any(is.na(as.character(insertion_data$chr)))
        !any(is.na(as.numeric(insertion_data$start)))
        !any(is.na(as.numeric(insertion_data$end)))
        !any(is.na(as.character(insertion_data$name)))
        !any(is.na(as.character(insertion_data$shape)))
        !any(is.na(as.numeric(insertion_data$in_tandem)))
        !any(is.na(as.character(insertion_data$colour)))

        ## values make sense
        all(insertion_data$end - insertion_data$start >= 0)
        all(insertion_data$length > 0)
        all(insertion_data$in_tandem > 0)
        all(insertion_data$shape %in% c("rectangle", "forward_arrow",
                                        "reverse_arrow"))

        ## we don't have any duplicates in insertion names
        length(insertion_data$name) == length(unique(insertion_data$name))
    })

    ## warn people about the factor colour column problem
    if (methods::is(insertion_data$colour, "factor")) {
        warning("The 'colour' column is a factor. Colours may not appear as
                you intended - if this is the case please change it to
                character")
    }

    ## can't have multiple chromosomes in one insertion diagram (use the MID
    ## instead)
    if (length(unique(insertion_data$chr)) > 1 & multiple == FALSE) {
        stop(
            "Sorry, having multiple different target sequences in one
            insertion diagram is not supported. Please instead use the
            multipleInsertionDiagram function to plot this.")
    }
    return(insertion_data)
}

## coverage data
.checkCoverageData <- function(coverage_data)
{
    ## type check and convert to data frame
    if (methods::is(coverage_data, "GRanges")) {
        coverage_data <- .GRangesToBed(coverage_data)

    } else if (!methods::is(coverage_data, "data.frame")) {
        stop("coverage_data should be either a data frame or GRanges")
    }

    ## check that the columns and the values are ok
    stopifnot(exprs={
        ## have all of the columns we need
        all(c("chr", "start", "end", "coverage") %in% colnames(coverage_data))

        ## columns are the right types
        !any(is.na(as.character(coverage_data$chr)))
        !any(is.na(as.numeric(coverage_data$start)))
        !any(is.na(as.numeric(coverage_data$end)))
        !any(is.na(as.numeric(coverage_data$coverage)))

        ## values make sense
        all(coverage_data$end - coverage_data$start >= 0)
        all(coverage_data$coverage >= 0)
    })
    return(coverage_data)
}

## check that the sector you're trying to plot into exists & warn if not
.checkSectorsMatch <- function(to_plot)
{
    if (any(!to_plot$chr %in% get.all.sector.index())) {
        warning("Not all of the rows in the provided data match a sector in the
                plot. Please double check that the spelling is exactly the same
                or these rows will not be plotted.")
    } else {
        TRUE
    }
}
