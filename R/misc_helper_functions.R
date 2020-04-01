#' @importFrom GenomicRanges reduce
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges strand

## turn 1 into chr1, X into chrX etc
.addChrToNumber <- function(add_to)
{
    add_to <- as.character(add_to)
    for (i in seq_along(add_to)) {
        test1 <- grepl("^.$", add_to[[i]])
        test2 <- grepl("^..$", add_to[[i]])

        if (test1 | test2) {
            add_to[[i]] <- paste0("chr", add_to[[i]]) #turn 1 into chr1 etc
        }
    }
    add_to <- unname(as.vector(add_to), force=TRUE)
    return(as.factor(add_to))
}

.makeRegex <- function(chrs)
{
    ## we should add $ so we don't match chr14 and chr22 for chr1 and chr2
    chrs <- paste0("^", chrs, "$")

    ## then collapse to single string with w/ |
    chrs <- paste(chrs, collapse = "|")
    return(chrs)
}

## convert GRanges to data frame
.GRangesToBed <- function(GRange, ideogram = FALSE)
{
    ## if it's the ideogram, then we only want chr, start and end
    if (ideogram == TRUE) {
        bed <- as.data.frame(GenomicRanges::reduce(GRange))
        bed <- subset(bed, select = c(seqnames, start, end))
        colnames(bed) <- c("chr", "start", "end")

    ## otherwise, take all columns except width and strand
    ## this assumes that the metadata columns provided are relevant
    } else {
        bed <- as.data.frame(GRange, row.names=NULL)
        bed <- subset(bed, select=-c(width, strand))
        new_col_names <- colnames(bed)
        new_col_names[1] <- "chr"
        colnames(bed) <- new_col_names

    }
    return(bed)
}

## from https://stackoverflow.com/questions/43635846/calculating-mean-for-
## every-n-values-from-a-vector
## honestly i don't understand this at all but it works lol
.makeWindows <- function(vector, every, na.rm = FALSE)
{
    vector_length <- length(vector)
    windowed_vector <- .colMeans(vector, every, vector_length %/% every, na.rm)
    remainder <- vector_length %% every
    if (remainder) {
        windowed_vector <- c(
            windowed_vector, mean.default(vector[(
                vector_length - remainder + 1):vector_length], na.rm=na.rm))
    }
    return(windowed_vector)
}

## vectorised rep and seq for use in coverage data import
.vectorisedSeq <- Vectorize(
    seq.default, SIMPLIFY=FALSE, vectorize.args=c("from", "length.out"))
.vectorisedRep <- Vectorize(
    rep_len, SIMPLIFY=FALSE, vectorize.args=c("x", "length.out"))

## subset out this the coverage for this particular sector
.getThisSectorCoverage <- function(coverage_data, currentSector)
{
    this_sector_coverage <- subset(
        coverage_data, as.character(coverage_data$chr) %in%
            as.character(currentSector))
    return(this_sector_coverage)
}

.insertionsToFeatures <- function(insertion_data, insertion_label="default")
{
    feature_data <- data.frame(stringsAsFactors = FALSE)
    for (i in seq_along(insertion_data$name)) {
        ## choose the label: name of insert or number of the insert
        if (insertion_data$in_tandem[i] == 1 && insertion_label ==
            "default") {
            labels <- insertion_data$name[i]
        } else if (insertion_label == "default") {
            labels <- seq(
                from = 1, by = 1, length.out = insertion_data$in_tandem[i])
        } else {
            labels <- insertion_label[i]
        }

        ## construct the df
        this_insertion <- data.frame(
            chr = rep(insertion_data$name[i], insertion_data$in_tandem[i]),
            start = seq(
                from = 0, by = insertion_data$length[i],
                length.out = insertion_data$in_tandem[i]),
            end = seq(
                from = insertion_data$length[i], by = insertion_data$length[i],
                length.out = insertion_data$in_tandem[i]),
            colour = rep(
                insertion_data$colour[i], insertion_data$in_tandem[i]),
            label = labels,
            shape = rep(insertion_data$shape[i], insertion_data$in_tandem[i]),
            track = rep(1, insertion_data$in_tandem[i]),
            stringsAsFactors = FALSE)
        feature_data <- rbind(feature_data, this_insertion)
    }
    return(feature_data)
}

.insertionsToIdeogram <- function(insertion_data, either_side)
{
    ## how long the original sequence sector should be is based on the start &
    ## end positions of the insert as well as either_side
    if (methods::is(either_side, "GRanges")) { # custom start/end as GRanges
        original_start <- as.numeric(min(BiocGenerics::start(either_side)))
        original_end <- as.numeric(max(BiocGenerics::end(either_side)))

    } else if (methods::is(either_side[1], "character") &
                either_side[1] == "default") { # default value
        either_side <- (max(insertion_data$end) - min(
            insertion_data$start)) * 0.15
        original_start <- min(insertion_data$start) - either_side
        original_end <- max(insertion_data$end) + either_side

    } else if (methods::is(either_side, "numeric") &
                length(either_side) == 1) { # a set number either side
        original_start <- min(insertion_data$start) - either_side
        original_end <- max(insertion_data$end) + either_side

    } else if (methods::is(either_side, "numeric") &
                length(either_side) == 2) { # custom start and end as vector
        original_start <- either_side[1]
        original_end <- either_side[2]

    } else {
        stop(
            "either_side should be one of: (1) 'default', (2) a single number,
            (3) a vector of length 2 or (4) a GRanges")
    }

    ## how long the inserted sequence sector should be depends on how many
    ## copies we have of the insertion
    insertion_start <- 0 - 0.4 * insertion_data$length
    insertion_end <- (insertion_data$length * insertion_data$in_tandem) +
        0.4 * insertion_data$length
    ideogram <- data.frame(
        chr = as.factor(c(
            as.character(insertion_data$chr[1]),
            as.character(insertion_data$name))),
        start = c(original_start, insertion_start),
        end = c(original_end, insertion_end), stringsAsFactors = FALSE)
    levels(ideogram$chr) <- c(
        as.character(insertion_data$chr[1]), as.character(insertion_data$name))
    return(ideogram)
}

## small little function to convert angle to one between 0 and 360 degrees
.simplifyAngle <- function(angle) {
    while(angle > 360){
        angle <- angle - 360
    }

    while(angle < 0){
        angle <- angle + 360
    }
    return(angle)
}
