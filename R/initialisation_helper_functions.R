## small helper function to add the zoomed row to the ideo data
## adapted from circlize book chapter 8.4 https://jokergoo.github.io/circlize_
## book/book/initialize-genomic-plot.html#zooming-chromosomes
.makeZoomedIdeogramData <- function(
    bed, chromosome, prefix = "zoomed_", remove_unzoomed=TRUE)
{
    zoom_bed <- bed[bed[[1]] %in% chromosome, , drop=FALSE]
    zoom_bed[[1]] <- paste0(prefix, zoom_bed[[1]])

    if (remove_unzoomed == TRUE) {
        bed <- subset(bed, !(bed[[1]] %in% chromosome), drop=TRUE)
    }

    zoomed <- rbind(bed, zoom_bed)
    zoomed$chr <- as.factor(zoomed$chr)
    return(zoomed)
}

## a function that will handle the rest of the zooming
.setupZoomedIdeogramData <- function(
    ideogram_data, zoom_sectors, prefix="zoomed_", zoom_size=0.055,
    remove_unzoomed=TRUE)
{
    ## make the zoomed ideogram data
    zoomedIdeo <- .makeZoomedIdeogramData(
        ideogram_data, zoom_sectors, prefix, remove_unzoomed)

    ## sort and reorder factor levels
    zoomedIdeo <- .sortIdeogramData(zoomedIdeo)

    ## setup a few things
    xrange <- zoomedIdeo$end - zoomedIdeo$start # chromosome length
    normal_chr_index <- grep(
        paste0("^", prefix), zoomedIdeo$chr, invert=TRUE)
    zoomed_chr_index <- grep(paste0("^", prefix), zoomedIdeo$chr)

    ## normalise lengths separately for zoomed and normal chromosomes
    normalSectors <- xrange[normal_chr_index] / sum(xrange[normal_chr_index]) *
        (1 - zoom_size * length(zoom_sectors))
    zoomedSectors <- xrange[zoomed_chr_index] / sum(xrange[zoomed_chr_index]) *
        (zoom_size * length(zoom_sectors))
    sector_width = c(normalSectors, zoomedSectors)
    return(list(zoomedIdeo, sector_width))
}

## helper function to sort the ideogram data
.sortIdeogramData <- function(ideogram_data)
{
    ## make a sort column from the number of the chromosome
    ideogram_data$sort <- suppressWarnings( # letters become NA -> sort last
        as.numeric(substr(ideogram_data$chr, 4, 50)))
    sorted_ideogram_data <- ideogram_data[order(ideogram_data$sort) ,]

    sorted_ideogram_data <- sorted_ideogram_data[, -4] # remove sort column
    if (methods::is(sorted_ideogram_data$chr, "factor")) {
        sorted_ideogram_data$chr <- droplevels(sorted_ideogram_data$chr)
    }
    return(sorted_ideogram_data)
}
