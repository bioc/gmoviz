#' @title Import coverage data from .bam file
#'
#' @description Uses RSamtools to import coverage data from .bam file and
#' format it appropriately for plotting with gmoviz.
#'
#' @param bam_file Location of the bam file from which to read coverage data.
#' @param regions_of_interest either a \linkS4class{GRanges} of regions OR a
#' character vector of sequences/chromosomes to find the coverage for (please
#' be careful that the names here match the spelling/format of those in the
#' bam file).
#' @param window_size The size of the window to for calculating coverage
#' (default is 1; per base coverage). Use \code{smoothCoverage} to smooth the
#' data, this is more for reducing time taken to read in and plot coverage
#' over a large number of bases.
#' @param smoothing_window_size If supplied, then moving average smoothing will
#' be applied using the \code{\link[pracma]{movavg}} function from the package
#' \pkg{pracma} (please make sure it's installed). Note: smoothing may take
#' some time when there are many points involved. Please be patient, or
#' proceed without smoothing.
#'
#' @export
#' @importFrom Rsamtools BamFile
#' @importFrom Rsamtools seqinfo
#' @importFrom Rsamtools ScanBamParam
#' @importFrom IRanges IRanges
#' @importFrom IRanges shift
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors elementNROWS
#' @importFrom BiocGenerics unlist
#' @importFrom GenomicAlignments coverage
#' @importFrom pracma movavg
#'
#' @return A \linkS4class{GRanges} containing the coverage data in the metadata
#' column 'coverage'.
#'
#' @seealso The \code{\link{gmovizInitialise}},
#' \code{\link{drawLinegraphTrack}}, \code{\link{insertionDiagram}} and
#' \code{\link{featureDiagram}} functions which can plot the coverage data.
#'
#' @examples
#' ## the example .bam file
#' path <- system.file('extdata', 'ex1.bam', package='Rsamtools')
#'
#' ## example without smoothing or windowing
#' getCoverage(regions_of_interest='seq1', bam_file=path)
#'
#' ## using windowing
#' getCoverage(regions_of_interest='seq1', bam_file=path, window_size=5)
#'
#' ## using smoothing
#' getCoverage(regions_of_interest='seq1', bam_file=path,
#' smoothing_window_size=3)
#'
#' ## specifying only a particular region to read in using GRanges
#' small_range <- GRanges('seq1', IRanges(0, 500))
#' getCoverage(regions_of_interest=small_range, bam_file=path)
#' ## please be very careful that the sequence names are spelled exactly like
#' ## in the bam file or you'll get an error! The following WON'T WORK.
#' \dontrun{
#' getCoverage(regions_of_interest='chr1', bam_file=path)}
#'
getCoverage <- function(regions_of_interest, bam_file,
    window_size = 1, smoothing_window_size = NULL) {
    ## check the inputs
    stopifnot(exprs = {
        methods::is(regions_of_interest, "GRanges") |
            methods::is(regions_of_interest, "character")
        methods::is(bam_file, "character")
        window_size >= 1
        is.null(smoothing_window_size) | smoothing_window_size >=
            1
    })

    ## read in the bam file
    bam <- Rsamtools::BamFile(bam_file)
    sequence_info <- Rsamtools::seqinfo(bam)

    ## if regions_of_interest is a character vector,
    ## make a GRanges for it
    if (methods::is(regions_of_interest, "character")) {
        if (!all(regions_of_interest %in% names(sequence_info))) {
            stop("Make sure all of the chromsomes in regions_of_interest are in
                the bam file and spelled exactly the same as in the bam")
        }
        coverage_range <- GRanges(regions_of_interest,
            IRanges::IRanges(start = rep(0, length(regions_of_interest)),
                end = GenomeInfoDb::seqlengths(sequence_info)[
                    regions_of_interest] - 1))

    } else {
        coverage_range <- regions_of_interest
    }

    ## check the sequence names given are in the bam
    ## file
    if (!all(GenomeInfoDb::seqlevels(coverage_range) %in%
        names(sequence_info))) {
        stop("Make sure all of the chromsomes in regions_of_interest are in the
            bam file and spelled exactly the same as in the bam")
    }

    ## get the coverage information
    coverage <- GenomicAlignments::coverage(bam,
        param = Rsamtools::ScanBamParam(which = coverage_range))
    if (start(coverage_range) == 0) {
        # idk why but this is important
        coverage_vector_list <- coverage[shift(coverage_range,
            1)]
    } else {
        coverage_vector_list <- coverage[coverage_range]
    }

    coverage_vector_list <- methods::as(coverage_vector_list,
        "IntegerList")

    ## apply windowing if needed
    if (window_size > 1) {
        coverage_vector_list <- lapply(coverage_vector_list,
            .makeWindows, every = window_size)
    }

    ## make the correctly formatted GRanges object
    starts <- unlist(.vectorisedSeq(from = start(coverage_range),
        by = window_size,
        length.out = S4Vectors::elementNROWS(coverage_vector_list)))
    ends <- unlist(.vectorisedSeq(from = (start(coverage_range) +
        window_size), by = window_size,
        length.out = S4Vectors::elementNROWS(coverage_vector_list)))
    coverage_data <- GRanges(
        seqnames = BiocGenerics::unlist(.vectorisedRep(
            names(coverage_vector_list),
        S4Vectors::elementNROWS(coverage_vector_list))),
        ranges = IRanges::IRanges(start = starts,
            end = ends), coverage = BiocGenerics::unlist(coverage_vector_list))

    ## smooth the coverage data if needed
    if (!is.null(smoothing_window_size)) {
        if (!requireNamespace("pracma", quietly = TRUE)) {
            stop("The package 'pracma' is needed for smoothing functionality.
                Please install it.",
                call. = FALSE)
        } else {
            # TODO: try to find something faster
            coverage_data$coverage <- pracma::movavg(
                x = coverage_data$coverage,
                n = smoothing_window_size, type = "s")
        }
    }

    ## change 1 to chr1 etc
    GenomeInfoDb::seqlevels(coverage_data) <- as.character(
        .addChrToNumber(GenomeInfoDb::seqlevels(coverage_data)))

    ## warn if there's more than 10-15k points
    if (length(coverage_data) > 10000) {
        warning("Your coverage_data has more than 10,000 points; this might
                take quite a while to plot. Consider increasing window_size to
                reduce the number of points.")
    }
    return(coverage_data)
}

#' @title Import transgenic genome data from .bam or .fasta file
#'
#' @description Read in the seqname, start & end from .bam or .fasta file and
#' format correctly for plotting with gmoviz.
#'
#' @param bam_file,fasta_file,fasta_folder Location of either a .bam file,
#' .fasta file or folder of .fasta files to read in. You only need to supply
#' one of these file types; .bam files are recommended because it is much
#' faster than using .fasta files. Also note that the filters
#' \code{unwatedChr}, \code{wanted_chr} and \code{just_pattern} won't work with
#' single .fasta files (only with .bam or .fasta folders).
#' @param unwanted_chr If supplied, these sequences won't be read in
#' @param wanted_chr If supplied, only these sequences will be read in
#' @param just_pattern If supplied, this pattern (regex) will be used to select
#' the sequences to read in
#' @param add_chr If \code{TRUE}, 'chr' will be added to the start of sequence
#' names with one or two characters (e.g. X will become chrX and 10 will become
#' chr10 but plasmid will remain as-is)
#'
#' @export
#' @importFrom Biostrings readDNAStringSet
#' @importFrom BiocGenerics width
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom IRanges IRanges
#' @importFrom Rsamtools BamFile
#' @importFrom GenomicRanges GRanges
#'
#' @return A \linkS4class{GRanges} containing the ideogram data (sequence
#' names, starts & ends).
#'
#' @seealso The \code{\link{gmovizInitialise}} and
#' \code{\link{featureDiagram}} functions which can be used to plot
#' this data.
#'
#' @examples
#' ## the example .bam file
#' path <- system.file('extdata', 'ex1.bam', package='Rsamtools')
#'
#' ## just starting with 'seq'
#' getIdeogramData(bam_file=path, just_pattern='^seq')
#'
#' ## only seq1
#' getIdeogramData(bam_file=path, wanted_chr='seq1')
#'
#' ## not seq2 (same as above)
#' getIdeogramData(bam_file=path, unwanted_chr='seq2')
#'
#' ## you can mix and match any of the filters
#' getIdeogramData(bam_file=path, unwanted_chr='seq2', just_pattern='^seq')
#'
#' ## the function also works to read in individual .fasta files, but please
#' ## note that for now the filters won't work (so if you have multiple
#' ## sequences in one .fasta file then they will all be read in)
#' path <- system.file('extdata', 'someORF.fa', package='Biostrings')
#' getIdeogramData(fasta_file=path)
#'
#' ## we can also read in selected .fasta files from a folder of .fasta files,
#' ## based on the filters shown above for the .bam file
#' path <- system.file('extdata', 'fastaFolder', package='gmoviz')
#' getIdeogramData(fasta_folder=path)

getIdeogramData <- function(bam_file = NULL, fasta_file = NULL,
    fasta_folder = NULL, just_pattern = NULL, unwanted_chr = NULL,
    wanted_chr = NULL, add_chr = TRUE) {
    ## check the inputs
    stopifnot(exprs = {
        methods::is(bam_file, "character") | is.null(bam_file)
        methods::is(fasta_file, "character") |
            is.null(fasta_file)
        methods::is(fasta_folder, "character") |
            is.null(fasta_folder)
        methods::is(just_pattern, "character") |
            is.null(just_pattern)
        methods::is(unwanted_chr, "character") |
            is.null(unwanted_chr)
        methods::is(wanted_chr, "character") |
            is.null(wanted_chr)
    })

    ## if we have only a fasta file, read it in
    if (!is.null(fasta_file)) {
        if (!requireNamespace("Biostrings", quietly = TRUE)) {
            stop("The package 'Biostrings' is needed to read .fasta files.
                Please install it.",
                call. = FALSE)
        } else {
            ## warn people in case they tried to use the
            ## filters
            if (!is.null(just_pattern) | !is.null(unwanted_chr) |
                !is.null(wanted_chr)) {
                warning("The filters 'just_pattern', 'unwanted_chr' and
                        'wanted_chr' only work for the 'fasta_folder' and
                        'bam_file' inputs. Reading in without any filters.")
            }
            fasta <- Biostrings::readDNAStringSet(fasta_file)
            return(GRanges(names(fasta), IRanges::IRanges(start = rep(0,
                length(names(fasta))), end = BiocGenerics::width(fasta))))
        }
    }

    ## otherwise, read in all of the sequences names
    ## and then filter
    if (!is.null(bam_file)) {
        sequence_info <- GenomeInfoDb::seqinfo(Rsamtools::BamFile(bam_file))
        sequence_names <- GenomeInfoDb::seqnames(sequence_info)

    } else if (!is.null(fasta_folder)) {
        sequence_names <- list.files(fasta_folder,
            pattern = "*.fasta", full.names = TRUE)
        if (length(sequence_names) == 0) {
            stop("Please make sure there are .fasta files inside the folder you
                have specified")
        } else {
            sequence_names <- as.vector(sequence_names)
        }

    } else {
        stop("Please supply either a bam file, a fasta file or a folder
            containing fasta files")
    }

    ## apply filters
    if (!is.null(just_pattern)) {
        sequence_names <- subset(sequence_names,
            grepl(just_pattern, sequence_names,
                ignore.case = TRUE))
    }

    if (!is.null(wanted_chr)) {
        sequence_names <- subset(sequence_names,
            grepl(.makeRegex(wanted_chr), sequence_names,
                ignore.case = TRUE))
    }

    if (!is.null(unwanted_chr)) {
        sequence_names <- subset(sequence_names,
            !grepl(.makeRegex(unwanted_chr), sequence_names,
                ignore.case = TRUE))
    }

    if (length(sequence_names) == 0) {
        stop("Sorry, there are no sequences to read in matching your filters.
            Please ensure that when using the just_pattern and wanted_chr
            argument that the spelling matches exactly the spelling in the .bam
            file or the name of the .fasta file, depending on what input you
            are using")
    }

    ## get the start/end of each sequence and format
    ## 'ideogram_data' as GRanges
    if (!is.null(bam_file)) {
        GenomeInfoDb::seqlevels(sequence_info) <- sequence_names
        ideogram_data <- GRanges(seqnames = sequence_names,
            ranges = IRanges::IRanges(start = rep(0,
                length(sequence_names)),
                end = GenomeInfoDb::seqlengths(sequence_info)))

    } else if (!is.null(fasta_folder)) {
        fasta <- Biostrings::readDNAStringSet(sequence_names)
        ideogram_data <- GRanges(seqnames = names(fasta),
            ranges = IRanges::IRanges(start = rep(0,
                length(names(fasta))), end = BiocGenerics::width(fasta)))
    }

    if (add_chr == TRUE) {
        GenomeInfoDb::seqlevels(ideogram_data) <- as.character(
            .addChrToNumber(GenomeInfoDb::seqlevels(ideogram_data)))
    }

    return(ideogram_data)
}

#' @title Generate a GRanges containing 'features' from .gff files
#'
#' @description Uses a .gff file to create a \linkS4class{GRanges} of
#' 'features' (e.g. genes or other regions of interest within the genome)
#' which can then be plotted with the \code{\link{featureDiagram}} or
#' \code{\link{drawFeatureTrack}} functions.
#'
#' @param gff_file Location of the gff file to read in.
#' @param colour_by_type If \code{TRUE}, the features will be coloured
#' according to the 'type' field of the gff file. If \code{FALSE}, colours will
#' be assigned based on the name of the feature (each uniquely named feature
#' gets its own colour).
#' @param colours A character vector of colours to be used to colour code
#' the features.
#'
#' @export
#' @importFrom rtracklayer readGFF
#' @importFrom GenomicRanges GRanges
#'
#' @return A \linkS4class{GRanges} containing the 'features'. See
#' \code{\link{drawFeatureTrack}} for a detailed description of the format.
#'
#' @seealso \code{\link{getLabels}} for a function which reads the entries
#' of a .gff file into labels rather than 'features'. Also
#' \code{\link{featureDiagram}} or \code{\link{drawFeatureTrack}} for
#' functions which can plot this data.
#'
#' @examples
#' ## the example .gff
#' path <- system.file('extdata', 'example.gff3', package='gmoviz')
#'
#' ## coloured by type
#' getFeatures(path)
#'
#' ## not coloured by type (each uniquely named feature gets its own colour)
#' getFeatures(path, colour_by_type=FALSE)

getFeatures <- function(gff_file, colours = nice_colours,
    colour_by_type = TRUE) {
    ## check the inputs
    stopifnot(exprs = {
        methods::is(gff_file, "character")
        methods::is(colours, "character")
        methods::is(colour_by_type, "logical")
    })

    ## read in the gff file
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
        stop("The package 'rtracklayer' is needed to read .gff files. Please
            install it.",
            call. = FALSE)
    } else {
        gff <- rtracklayer::readGFF(gff_file)
    }

    ## set up the correctly formatted GRanges
    feature_data <- GRanges(seqnames = gff$seqid,
        ranges = IRanges::IRanges(gff$start, gff$end),
        label = gff$ID, track = rep(1, length(gff$seqid)))

    ## only want the type column if we are
    ## colour-coding by type:
    if (colour_by_type == TRUE) {
        feature_data$type <- gff$type
    }

    ## colour by type or by name
    shape_column <- vector(length = length(feature_data))
    colour_column <- vector(length = length(feature_data))
    if (colour_by_type == TRUE) {
        colour_by_these <- gff$type
        types <- unique(gff$type)
    } else {
        colour_by_these <- feature_data$label
        types <- unique(feature_data$label)
    }

    ## check that we have enough colours to assign
    stopifnot(length(colours) >= length(types))

    for (i in seq_along(feature_data$label)) {
        ## assign shapes: arrows for genes, otherwise
        ## rectangles
        if (grepl(gff$type[i], "gene", ignore.case = TRUE)) {
            if (is.na(gff$strand[i])) {
                shape_column[i] <- "forward_arrow"
            } else if (gff$strand[i] == "-") {
                shape_column[i] <- "reverse_arrow"
            } else {
                shape_column[i] <- "forward_arrow"
            }
        } else {
            shape_column[i] <- "rectangle"
        }

        ## assign colours: either based on the 'type'
        ## specified in the gff OR give every uniquely
        ## named feature its own colour
        for (j in seq_along(types)) {
            if (colour_by_these[i] == types[j]) {
                colour_column[i] <- colours[j]
            }
        }
    }
    feature_data$shape <- as.character(shape_column)
    feature_data$colour <- as.character(colour_column)
    return(feature_data)
}

#' @title Generate a GRanges of labels from .gff files
#'
#' @description Uses a .gff file to create a GRanges of labels which can then
#' be plotted with the \code{label_data} argument of many functions in this
#' package such as \code{\link{gmovizInitialise}},
#' \code{\link{insertionDiagram}} or \code{\link{featureDiagram}}.
#'
#' @param gff_file Location of the gff file to read in.
#' @param colour_code If \code{TRUE}, the labels will be assigned colours
#' according to the 'type' field of the gff file. If \code{FALSE}, colours will
#' not be assigned.
#' @param colours A character vector of colours to be used to colour code
#' the labels (if \code{colour_code} is \code{TRUE}).
#'
#' @export
#' @importFrom rtracklayer readGFF
#' @importFrom GenomicRanges GRanges
#'
#' @return A GRanges containing the gene label data. See
#' \code{\link{gmovizInitialise}} for a detailed description of the format.
#'
#' @seealso \code{\link{getFeatures}} for a function which reads the entries
#' of a .gff file into 'features' rather than labels. Also
#' \code{\link{gmovizInitialise}}, \code{\link{insertionDiagram}} and
#' \code{\link{featureDiagram}} for functions which can plot this data.
#'
#' @examples
#' ## example .gff
#' path <- system.file('extdata', 'example.gff3', package='gmoviz')
#'
#' ## colour coded
#' getLabels(path)
#'
#' ## not colour coded (all black)
#' getLabels(path, colour_code=FALSE)

getLabels <- function(gff_file, colour_code = TRUE,
    colours = bright_colours_opaque) {
    ## check inputs:
    stopifnot(exprs = {
        methods::is(colour_code, "logical")
        methods::is(colours, "character")
    })

    ## read in gff file
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
        stop("The package 'rtracklayer' is needed to read .gff files. Please
            install it.",
            call. = FALSE)
    } else {
        gff <- rtracklayer::readGFF(gff_file)
    }

    ## set up the correctly formatted GRanges
    label_data <- GRanges(seqnames = gff$seqid,
        ranges = IRanges::IRanges(gff$start, gff$end),
        label = gff$ID, track = rep(1, length(gff$seqid)))

    ## assign colours
    if (colour_code == TRUE) {
        label_colour <- vector(length = length(label_data))
        colour_by_these <- unique(gff$type)
        for (i in seq_along(label_data$label)) {
            for (j in seq_along(colour_by_these)) {
                if (gff$type[i] == colour_by_these[j]) {
                    label_colour[i] <- colours[j]
                }
            }
        }
        label_data$colour <- as.character(label_colour)
    }
    return(label_data)
}

#' @title Initialise the layout of the circular plot
#'
#' @description Draws the ideogram (sectors around a circle representing
#' sequences of interest, like chromosomes), labels and genomic axis in
#' preparation for the addition of other tracks like
#' \code{\link{drawFeatureTrack}} or \code{\link{drawLinegraphTrack}}.
#'
#' @param ideogram_data Either a \linkS4class{GRanges} representing regions of
#' interest or a data frame in bed format (containing the \code{chr},
#' \code{start} and \code{end} columns). If you want to read in data from file,
#' please see the \code{\link{getIdeogramData}} function.
#' @param start_degree Where on the circle the first sector will start being
#' drawn from (90 = 12 o'clock).
#' @param space_between_sectors Space between each sector.
#' @param zoom_sectors A character vector of sectors that should be 'zoomed'
#' (made bigger than usual, useful to show shorter sequences like plasmids
#' alongside longer sequences like chromosomes).
#' @param zoom_size The size of the zoomed chromosome, as a proportion of the
#' entire circle (0 = invisible, 1 = entire circle filled). The default value
#' of \code{0.055} is good for displaying something small (e.g. plasmid)
#' alongside something large (e.g. chromosomes).
#' @param remove_unzoomed If \code{TRUE}, the sectors in \code{zoom_sectors}
#' will only appear in their zoomed form. If \code{FALSE}, both the zoomed
#' and unzoomed versions will be plotted.
#' @param zoom_prefix A character prefix that will be applied to zoomed
#' sequences to distinguish them from non-zoomed ones.
#' @param custom_sector_width Normally, the size of each sector is proportional
#' to its relative length, but \code{custom_sector_width} can change this. It
#' is a vector of sector sizes (as proportions of the entire circle), given in
#' the same order in which sectors are plotted: firstly 'chr1', 'chr2' ...
#' through to 'chrX' and 'chrY' followed by any differently named sectors e.g.
#' 'gene 1', plasmid' in alphabetical order.
#' @param track_height The height (vertical distance around the circle) that
#' will be taken up by this track. Should be a number between 0 (none) and 1
#' (entire circle).
#' @param sector_colours Either a single colour (which will be applied to all
#' sectors) or a vector with the same length as the number of sectors/regions.
#' This package includes 5 colour sets: \code{nice_colours},
#' \code{pastel_colours}, \code{bright_colours_transparent},
#' \code{bright_colours_opaque} and \code{rich_colours}. See
#' \code{\link{colourSets}} for more information about these.
#' @param sector_border_colours Same as \code{sector_colours}, only for the
#' border of each sector.
#' @param sector_labels If \code{TRUE}, labels ('chr1', 'chr2' etc.) will be
#' drawn for each sector (recommended).
#' @param sector_label_size Size of the sector labels.
#' @param sector_label_colour Colour of the sector labels.
#' @param xaxis If \code{TRUE}, an x (genomic position) xaxis  will be plotted.
#' @param xaxis_spacing Space between the x axis labels, in degrees.
#' Alternatively, the string 'start_end' will place a label at the start and
#' end of each sector only.
#' @param xaxis_spacing_unit Either \code{"deg"} to draw ticks every certain
#' number of degrees around the circle or \code{"bp"} to draw ticks every
#' certain bp around the circle (be warned that when the scales for each sector
#' are very different, it's best to use \code{"deg"})
#' @param xaxis_orientation Either \code{'top'} to put the x axis on the
#' outside of the circle or \code{'bottom'} to put it on the inside.
#' @param xaxis_label_size Size of the x axis labels.
#' @param xaxis_colour Colour of the x axis labels.
#' @param label_data Data frame or \linkS4class{GRanges} containing the
#' labels. If a GRanges, \code{label} should be a metadata column containing
#' the character strings of the labels. \code{type} and \code{colour} can also
#' be used to store additional information about the type (e.g. 'gene' or
#' 'promoter') and colour of the label. This information can be used to
#' \bold{colour code} the labels by supplying the colour column as the
#' \code{label_colour} parameter. Data frames should additionally include the
#' \code{chr}, \code{start}, \code{end} which dictate the position of the
#' label.
#' @param label_colour Colour of the labels, can be either a single value
#' (applied to all labels) or a vector with the same length as the number of
#' labels (for colour-coding).
#' @param label_size Size of the labels.
#' @param space_between_labels Space between the labels
#' @param label_orientation \code{'outside'} to put the labels on the
#' outside of the circle, \code{'inside'} to put them on the inside.
#' @param coverage_rectangle A vector containing the name(s) of any sector(s)
#' that you would like to depict as 'coverage rectangles': filled in shapes
#' that are a plot of the coverage data over that sector. See the example below
#' or the vignette for an example of this.
#' @param coverage_data A GRanges (or data frame) containing the coverage data
#' to plot for those sectors in \code{coverage_rectangle}. To read this data in
#' from a BAM file, please see the \code{\link{getCoverage}} function.
#' @param custom_ylim A vector of length two containing the y (coverage) axis
#' limits. No need to set if not using coverage rectangles or if you're happy
#' with the default: c(0, maximum coverage).
#' @param sort_sectors If \code{TRUE}, the sectors will be plotted around the
#' circle in alphabetical order. Otherwise, they will be in the order in which
#' they appear in \code{ideogram_data}
#' @export
#' @import circlize
#'
#' @return Generates an image of the initial ideogram track which can then be
#' added to with various other functions.
#' @seealso The \code{\link{drawScatterplotTrack}},
#' \code{\link{drawFeatureTrack}} and \code{\link{drawLinegraphTrack}}, which
#' can be used to add information to this plot. Also
#' \code{\link{getIdeogramData}} which can be used to read in the needed
#' ideogram data for this function.
#'
#' @examples
#' ## normal/standard usage
#' ideogram <- data.frame(chr=paste0('chr', c(1:19, 'X', 'Y')),
#' start=rep(0, 21),
#' end=c(195471971, 182113224, 160039680, 156508116, 151834684, 149736546,
#' 145441459, 129401213, 124595110, 130694993, 122082543, 120129022,
#' 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566,
#' 171031299, 91744698))
#' gmovizInitialise(ideogram)
#'
#' ## zooming a sector
#' gmovizInitialise(ideogram, zoom_sectors='chr19', zoom_size=0.2)
#'
#' ## custom sector width
#' small_ideogram <- data.frame(chr=c('region 1', 'region 2', 'region 3'),
#' start=c(0, 0, 0), end=c(10000, 12000, 10000))
#' gmovizInitialise(small_ideogram, custom_sector_width=c(0.3, 0.3, 0.3))
#'
#' ## coverage rectangle
#' path <- system.file('extdata', 'ex1.bam', package='Rsamtools')
#' ideo <- getIdeogramData(path, wanted_chr='seq1')
#' coverage <- getCoverage(bam_file=path, regions_of_interest='seq1',
#' window_size=30)
#' gmovizInitialise(ideo, coverage_rectangle='seq1', coverage_data=coverage)

gmovizInitialise <- function(ideogram_data, start_degree = 90,
    space_between_sectors = 1, zoom_sectors = NULL,
    zoom_size = 0.055, remove_unzoomed = TRUE,
    zoom_prefix = "zoomed_", custom_sector_width = NULL,
    track_height = 0.1, sector_colours = nice_colours,
    sector_border_colours = nice_colours, coverage_rectangle = NULL,
    coverage_data = NULL, custom_ylim = NULL, sector_labels = TRUE,
    sector_label_size = 0.9, sector_label_colour = "black",
    xaxis = TRUE, xaxis_orientation = "top", xaxis_label_size = 0.75,
    xaxis_colour = "#747577", xaxis_spacing = 10, xaxis_spacing_unit = "deg",
    label_data = NULL, label_colour = "black",
    label_size = 0.85, space_between_labels = 0.4,
    label_orientation = "outside", sort_sectors = TRUE) {
    ## check the the data
    ideogram_data <- .checkIdeogramData(ideogram_data)
    if (!is.null(coverage_rectangle) & !is.null(coverage_data)) {
        coverage_data <- .checkCoverageData(coverage_data)
    } else if (!is.null(coverage_rectangle) & is.null(coverage_data)) {
        stop("If you want to represent the coverage of the sectors given in
            coverage_rectangle, you need to supply a GRanges or data frame
            containing the coverage to coverage_data. See ?getCoverage for
            help on how to do this")
    }

    ## check other inputs
    stopifnot(exprs = {
        start_degree >= 0 & start_degree <= 360
        space_between_sectors >= 0
        is.null(zoom_sectors) | all(zoom_sectors %in%
            ideogram_data$chr)
        zoom_size > 0 & zoom_size < 1
        methods::is(remove_unzoomed, "logical")
        methods::is(zoom_prefix, "character")
        is.null(custom_sector_width) | length(custom_sector_width) ==
            nrow(ideogram_data)
        track_height > 0 & track_height < 1
        methods::is(sector_colours, "character")
        methods::is(sector_border_colours, "character")
        is.null(coverage_rectangle) | methods::is(coverage_rectangle,
            "character")
        is.null(custom_ylim) | length(custom_ylim ==
            2)
        methods::is(sector_labels, "logical")
        sector_label_size >= 0
        xaxis_orientation %in% c("top", "bottom")
        xaxis_label_size >= 0
        (methods::is(xaxis_spacing, "numeric") &
            xaxis_spacing > 0 | xaxis_spacing == "start_end")
        xaxis_spacing_unit %in% c("deg", "bp")
        methods::is(label_colour, "character")
        label_size >= 0
        space_between_labels > 0 & space_between_labels <
            1
        label_orientation %in% c("inside", "outside")
    })

    ## initialisation/setup
    circos.clear()
    par(xpd = NA)
    circos.par(start.degree = start_degree, unit.circle.segments = 1000,
        gap.after = space_between_sectors)

    ## if zooming is required, do it now:
    if (!is.null(zoom_sectors)) {
        zoomData <- .setupZoomedIdeogramData(ideogram_data,
            zoom_sectors, prefix = zoom_prefix,
            zoom_size, remove_unzoomed = remove_unzoomed)
        ideogram_data <- zoomData[[1]]
        circos.initializeWithIdeogram(ideogram_data,
            plotType = NULL, sort.chr = FALSE,
            sector.width = unlist(zoomData[2]))
        ## custom sector widths
    } else if (!is.null(custom_sector_width)) {
        ## need to manually order ideogram data before
        ## plotting in this case
        ideogram_data <- suppressWarnings(.sortIdeogramData(ideogram_data))
        circos.initializeWithIdeogram(ideogram_data,
            plotType = NULL, sector.width = custom_sector_width,
            sort.chr = FALSE)
        ## regular plotting
    } else {
        if (sort_sectors == TRUE){
            circos.initializeWithIdeogram(ideogram_data, plotType = NULL)
        } else {
            circos.initializeWithIdeogram(ideogram_data, plotType = NULL,
                                          sort.chr = FALSE)
        }

    }

    ## plot labels, if needed
    if (!is.null(label_data)) {
        label_data <- .checkLabelData(label_data)
        circos.genomicLabels(bed = label_data,
            labels.column = 4, side = label_orientation,
            col = label_colour, cex = label_size,
            line_col = label_colour)
    }

    ## draw the ideogram
    .plotNewTrack(current_sectors = ideogram_data$chr,
        custom_ylim, coverage_rectangle, coverage_data,
        track_height = track_height)

    .plotIdeogram(current_sectors = get.all.sector.index(),
        coverage_data, coverage_rectangle,
        current_track = CELL_META$track.index,
        sector_colours, sector_border_colours,
        sector_labels, sector_label_size, sector_label_colour)

    if (!is.null(xaxis_colour)) {
        .plotXaxis(xaxis, xaxis_label_size, xaxis_colour,
            xaxis_spacing, sectors = as.character(ideogram_data$chr),
            xaxis_spacing_unit = xaxis_spacing_unit)
    }
}

#' @title Add a scatterplot track to an existing plot
#'
#' @description Adds a scatterplot track to the existing plot. Must have
#' initialised the circular plot (by \code{\link{gmovizInitialise}} first).
#'
#' @param plot_data Either: (1) a \linkS4class{GRanges} object with a metadata
#' column of y values to plot OR (2) a data frame with four columns; \code{chr}
#' (should match those supplied when initialising the plot); \code{start} and
#' \code{end} (x values of the point: can both be the same if you only have
#' a single x value for position) and then a fourth column of y values.
#' @param track_border_colour Colour of the border of the plotting region.
#' @param track_height The proportion (between 0 and 1) of the circle
#' taken up by this track.
#' @param point_bicolour_cutoff A numeric threshold for the colour of the
#' points (points above/below this number will be different colours).
#' @param point_colour The fill colour of the points. If \code{
#' point_bicolour_cutoff != NULL} then this should be a vector with two
#' elements.
#' @param point_outline_colour The colour of the outline of the points. If
#' using point_bicolour_cutoff then this should be a vector with two elements.
#' @param point_size Size of the points.
#' @param point_type Type (shape) of the points, same as base R.
#' @param ylim Vector of length 2; upper and lower limits for y axis.
#' @param yaxis_increment The increment the y axis and gridlines will use.
#' @param show_yaxis If \code{TRUE}, a y axis will be drawn.
#' @param yaxis_label_size Size of the labels on the y axis.
#' @param yaxis_tick_size Size of the ticks on the y axis.
#' @param yaxis_location Sector the y axis is drawn on.
#' @param yaxis_side Side of the sector the y axis is on; either \code{'left'}
#' or \code{'right'}.
#' @param yaxis_colour Colour of the y axis.
#' @param show_gridlines If \code{TRUE} then gridlines will be drawn.
#' @param gridline_colour Colour of the gridlines.
#'
#' @return Adds a scatterplot track to existing visualisation.
#' @seealso \code{\link{gmovizInitialise}}, which must be used to initialise
#' the graph before this function. Also \code{\link{drawLinegraphTrack}} for a
#' similar function which displays data as a line graph instead.
#'
#' @export
#' @import circlize
#'
#' @examples
#' ## you must initialise first!
#' small_ideo <- data.frame(chr=c('region 1', 'region 2', 'region 3'),
#'                          start=c(0, 0, 0), end=c(10000, 12000, 10000))
#' gmovizInitialise(small_ideo, custom_sector_width=c(0.3, 0.3, 0.3))
#'
#' ## make the data
#' smallplot_data <- data.frame(
#' chr = sample(c('region 1', 'region 2','region 3'), size=40, replace=TRUE),
#' start = seq(0, 10000, length.out=40), end = seq(0, 10000, length.out=40),
#' val = rnorm(40, 2, 0.5))
#'
#' ## scatterplot where all points are the same colour
#' drawScatterplotTrack(smallplot_data)
#'
#' ## scatterplot with bi-colour cutoff of 2
#' drawScatterplotTrack(smallplot_data, point_bicolour_cutoff=2,
#'                      point_colour=c('red', 'blue'),
#'                      point_outline_colour=c('black', 'black'))

drawScatterplotTrack <- function(plot_data, track_border_colour = "black",
    track_height = 0.3, point_bicolour_cutoff = NULL,
    point_colour = "black", point_outline_colour = "black",
    point_size = 0.55, point_type = 21, ylim = NULL,
    yaxis_increment = NULL, show_yaxis = TRUE,
    yaxis_label_size = 0.6, yaxis_tick_size = 0.5,
    yaxis_location = CELL_META$sector.index, yaxis_side = "left",
    yaxis_colour = "black", show_gridlines = TRUE,
    gridline_colour = "#aaaaaa") {
    ## check inputs
    stopifnot(exprs = {
        methods::is(plot_data, "GRanges") | methods::is(plot_data,
            "data.frame")
        track_height > 0 & track_height < 1
        is.null(point_bicolour_cutoff) | methods::is(point_bicolour_cutoff,
            "numeric")
        length(point_colour) == 1 | length(point_colour) ==
            2
        length(point_outline_colour) == 1 | length(point_outline_colour) ==
            2
        point_size > 0
        point_type %in% 0:25
        (length(ylim) == 2 & methods::is(ylim,
            "numeric")) | is.null(ylim)
        methods::is(show_yaxis, "logical")
        yaxis_label_size >= 0
        yaxis_tick_size >= 0
        yaxis_location %in% get.all.sector.index()
        yaxis_side %in% c("left", "right")
        methods::is(show_gridlines, "logical")
    })

    ## convert into a data frame
    if (methods::is(plot_data, "GRanges")) {
        plot_data <- .GRangesToBed(plot_data)
    }

    ## let people know if they're plotting into a
    ## sector that doesn't exist
    .checkSectorsMatch(plot_data)

    ## set default ylim and yaxis_increments if
    ## needed
    if (is.null(ylim)) {
        ylim <- c(min(plot_data[, 4]), max(plot_data[,
            4]))
    }

    if (is.null(yaxis_increment)) {
        yaxis_increment <- (max(plot_data[, 4]) -
            min(plot_data[, 4]))/4
        yaxis_increment <- round(yaxis_increment,
            digits = 1)
    }
    stopifnot(yaxis_increment <= ylim[2] | yaxis_increment >=
        0)

    ## make the plotting track
    circos.genomicTrackPlotRegion(plot_data, ylim = c(ylim[1],
        ylim[2]), bg.border = track_border_colour,
        panel.fun = function(region, value, ...) {
            .plotGridlines(ylim, yaxis_increment,
                show_gridlines, gridline_colour)
            ## apply the colour cutoff for points, if
            ## necessary
            if (!is.null(point_bicolour_cutoff)) {
                point_colour <- ifelse(value[1] >
                    point_bicolour_cutoff, point_colour[1],
                    point_colour[2])
                point_outline_colour <- ifelse(value[1] >
                    point_bicolour_cutoff, point_outline_colour[1],
                    point_outline_colour[2])
            }

            circos.genomicPoints(region, value,
                bg = point_colour, pch = point_type,
                col = point_outline_colour, cex = point_size)

        }, track.height = track_height)

    .plotYaxis(show_yaxis, ylim, yaxis_increment,
        yaxis_label_size, yaxis_tick_size, yaxis_side,
        yaxis_location, yaxis_colour)
}

#' @title Add a line graph track to an existing plot
#'
#' @description Adds a line graph track to the existing plot. Must have
#' initialised the circular plot (by \code{\link{gmovizInitialise}} first).
#'
#' @inheritParams drawScatterplotTrack
#' @param line_shade_colour The colour the will be used to fill in under the
#' line. Set this to NULL if you just want the line rather than the area.
#' @param line_colour The colour of the line itself.
#'
#' @export
#' @import circlize
#'
#' @return Adds a line graph track to existing visualisation.
#' @seealso \code{\link{gmovizInitialise}}, which must be used to initialise
#' the graph before this function. Also \code{\link{drawScatterplotTrack}} for
#' a similar function which displays data as a scatterplot rather than as a
#' line graph.
#'
#' @examples
#' ## you must initialise first!
#' small_ideo <- data.frame(chr=c('region 1', 'region 2', 'region 3'),
#'                          start=c(0, 0, 0), end=c(10000, 12000, 10000))
#' gmovizInitialise(small_ideo, custom_sector_width=c(0.3, 0.3, 0.3))
#'
#' ## make the data
#' smallplot_data <- data.frame(
#' chr = sample(c('region 1', 'region 2','region 3'), size=300, replace=TRUE),
#' start = seq(0, 10000, length.out=300), end = seq(0, 10000, length.out=300),
#' val = rnorm(300, 2, 0.5))
#' ## line graph with no shading (just the line)
#' drawLinegraphTrack(smallplot_data, line_shade_colour=NULL)
#'
#' ## line graph with shading (a filled in shape)
#' drawLinegraphTrack(smallplot_data, line_shade_colour='#db009599')

drawLinegraphTrack <- function(plot_data, track_border_colour = "black",
    track_height = 0.3, yaxis_increment = NULL,
    ylim = NULL, line_shade_colour = "#5ab4ac",
    line_colour = "black", yaxis_label_size = 0.5,
    show_yaxis = TRUE, yaxis_tick_size = 0.4, yaxis_side = "left",
    yaxis_colour = "black", yaxis_location = CELL_META$sector.index,
    show_gridlines = TRUE, gridline_colour = "#aaaaaa") {

    ## check inputs
    stopifnot(exprs = {
        methods::is(plot_data, "GRanges") | methods::is(plot_data,
            "data.frame")
        track_height > 0 & track_height < 1
        (length(ylim) == 2 & methods::is(ylim,
            "numeric")) | is.null(ylim)
        methods::is(show_yaxis, "logical")
        yaxis_label_size >= 0
        yaxis_tick_size >= 0
        yaxis_location %in% get.all.sector.index()
        yaxis_side %in% c("left", "right")
        methods::is(show_gridlines, "logical")
    })

    ## convert granges into a data frame
    if (methods::is(plot_data, "GRanges")) {
        plot_data <- .GRangesToBed(plot_data)
    }

    ## let people know if they're plotting into a
    ## sector that doesn't exist
    .checkSectorsMatch(plot_data)

    ## set defaults for ylim and yaxis_increment, if
    ## needed
    if (is.null(ylim)) {
        ylim <- c(min(plot_data[, 4]), max(plot_data[,
            4]))
    }

    if (is.null(yaxis_increment)) {
        yaxis_increment <- (max(plot_data[, 4]) -
            min(plot_data[, 4]))/4
        yaxis_increment <- round(yaxis_increment,
            digits = 1)
    }
    stopifnot(yaxis_increment <= ylim[2] | yaxis_increment >=
        0)

    ## the plotting track
    circos.genomicTrackPlotRegion(data = plot_data,
        ylim = c(ylim[1], ylim[2]), bg.border = track_border_colour,
        panel.fun = function(region, value, ...) {
            .plotGridlines(ylim, yaxis_increment,
                show_gridlines, gridline_colour)

            ## are we shading the area under the curve?
            if (is.null(line_shade_colour)) {
                shade_under <- FALSE
                line_shade_colour <- line_colour

            } else {
                shade_under <- TRUE
            }

            ## plot the line graph
            circos.genomicLines(region, value,
                area = shade_under, col = line_shade_colour,
                border = line_colour)

        }, track.height = track_height)

    .plotYaxis(show_yaxis, ylim, yaxis_increment,
        yaxis_label_size, yaxis_tick_size, yaxis_side,
        yaxis_location, yaxis_colour)
}

#' @title Add a 'feature' track to an existing plot
#'
#' @description Adds to an existing plot a track which displays 'features'
#' (e.g. genes, indels, primer sequences etc) using coloured shapes. Note that
#' you must have initialised the circular plot (by
#' \code{\link{gmovizInitialise}} first).
#'
#' @param feature_data A data frame or \linkS4class{GRanges} containing the
#' 'features' to plot. \itemize{
#' \item GRanges input should contain \code{label}, \code{colour},
#' \code{shape} and \code{track} as metadata columns.
#' \item Data frame should contain \code{label}, \code{colour}, \code{shape}
#' and \code{track}, as well as the additional columns \code{chr}, \code{start}
#' and \code{end}} Please see below for a detailed description of these
#' columns, and \code{\link{getFeatures}} for a function which can read this
#' information in from a .gff file.
#' @param flipped_sector A vector of sectors that will have their genomic
#' position (x values) reversed (ascending in anti-clockwise direction, as
#' opposed to the usual ascending in a clock-wise direction).
#' @param feature_label_cutoff To enhance readability when the shapes are
#' small, those labels belonging to features smaller than
#' \code{feature_label_cutoff} will instead be plotted on a new track closer
#' to the centre of the circle, rather than inside the shapes themselves.
#' @param track_height The height (proportion of the circle) taken up by
#' \bold{each track} of features. The default value of \code{0.1} is
#' appropriate for up to 2 feature tracks; if you get an error due to running
#' out of space please reduce this.
#' @param feature_label_size Size of the feature labels.
#' @param label_track_height Size of the track on which to plot the labels.
#' @param coverage_rectangle,coverage_data If, when initialising the graph you
#' have used coverage_rectangleangle AND you want to plot features on the
#' outermost track (track 0), please fill these in the same as in your
#' \code{gmovizInitialise} function call. \strong{Otherwise, there is no need
#' to supply these}.
#' @param feature_outline Should a black outline be drawn around the feature
#' shape? (It is recommended to set this to \code{FALSE} when dealing with
#' very small features)
#' @param internal For internal use only.
#'
#' @section Feature data format:
#' The feature data \linkS4class{GRanges} contains four metadata columns:
#' \describe{
#' \item{label}{A character string which will be used to label the feature. It
#' is suggested to keep this label relatively short, if possible.}
#' \item{colour}{A character string of a colour to use. Supports hex colours
#' (\emph{e.g. #000000}) and named R colours (\emph{e.g. red}).}
#' \item{shape}{The shape that will be used to represent the feature: \itemize{
#' \item \code{'rectangle'}
#' \item \code{'forward_arrow'}
#' \item \code{'reverse_arrow'}
#' \item \code{'upwards_triangle'} (out of the circle).
#' \item \code{'downwards_triangle'} (into the circle).} It is suggested to use
#' \code{'forward_arrow'} for genes on the forward strand and
#' \code{'reverse_arrow'} for genes on the reverse strand.}
#' \item{track}{The index of the track on which to plot the feature: \itemize{
#' \item 0 represents the outermost track, where the ideogram rectangles that
#' represent sequences/chromosomes are plotted.
#' \item 1 is the conventional (default) track on which to plot a feature.
#' \item 2, 3 and so on are further into the centre of the circle.} It is
#' strongly recommended to keep the tracks below 3, otherwise there may not be
#' enough space in the circle to fit them all.}}
#' These columns are all \strong{optional}. If you don't supply them, then
#' default values will be added as follows: \describe{
#' \item{label}{\code{''}}
#' \item{colour}{a colour allocated from \code{\link{rich_colours}}}
#' \item{shape}{\code{'rectangle'}}
#' \item{track}{\code{1}}
#' }
#'
#' @export
#' @import circlize
#'
#' @return Adds a 'feature' track to an existing plot.
#' @seealso \code{\link{featureDiagram}} for a function that, while
#' slightly less flexible, generates an entire visualisation in one go. Also
#' \code{\link{getFeatures}} for a function that can read the feature data in
#' from a .gff file.
#'
#' @examples
#' ## plasmid map
#' plasmid_ideogram <- data.frame(chr='plasmid', start=0, end=2500)
#'
#' plasmid_features <- GRanges(seqnames=rep('plasmid', 4),
#' ranges=IRanges(start=c(0, 451, 901, 1700), end=c(450, 900, 1400, 2200)),
#' colour = c('#d44a9f', '#4a91d4', '#7ad44a', '#d49d4a'),
#' label = c('promoter', 'gene', 'GFP', 'ampR'),
#' shape = c('rectangle', 'forward_arrow', 'forward_arrow', 'reverse_arrow'),
#' track = rep(1, 4))
#'
#' ## for a simple case like this you might as well use the featureDiagram
#' ## function because it's only 1 function call, whereas here we need two:
#' gmovizInitialise(plasmid_ideogram)
#' drawFeatureTrack(plasmid_features)
#'
#' ## however the drawFeatureTrack function allows more flexibility e.g. if you
#' ## want to add features to a plot containing numerical data for example:
#' ## data
#' scatter_data <- GRanges(rep('plasmid', 50),
#' IRanges(start=sample(1:3000, 50), width=2),
#' scatter=rnorm(50, mean=4, sd=1))
#'
#' ## plotting
#' gmovizInitialise(plasmid_ideogram)
#' drawScatterplotTrack(plot_data=scatter_data)
#' drawFeatureTrack(plasmid_features, track_height = 0.15)

drawFeatureTrack <- function(feature_data, flipped_sector = NULL,
    feature_label_cutoff = 50, track_height = 0.1,
    feature_label_size = 0.9, label_track_height = 0.1 *
        feature_label_size, coverage_rectangle = NULL,
    coverage_data = NULL, internal = FALSE, feature_outline = TRUE) {

    ## check inputs
    feature_data <- .checkFeatureData(feature_data)
    stopifnot(exprs = {
        is.null(flipped_sector) | methods::is(flipped_sector,
            "character")
        #feature_label_cutoff >= 0
        feature_label_size > 0
        track_height > 0 & track_height < 1
        label_track_height > 0 & label_track_height <
            1
        is.null(coverage_rectangle) | methods::is(coverage_rectangle,
            "character")
    })

    # let people know if they're plotting into a
    # sector that doesn't exist
    .checkSectorsMatch(feature_data)

    ## plot the outermost track features
    outer_track_features <- subset(feature_data,
        feature_data$track == 0)
    if (nrow(outer_track_features) != 0) {
        for (i in seq_along(outer_track_features$label)) {
            ## if it's in a 'coverage rectangle'
            if (as.character(outer_track_features$chr[i]) %in%
                as.character(coverage_rectangle)) {
                ## basically we're re-drawing the coverage again
                ## here to be the colour of the feature
                this_sector_coverage <- coverage_data[coverage_data$chr ==
                    outer_track_features$chr[i], ]
                this_bit_of_coverage <- this_sector_coverage[
                    this_sector_coverage$start %in%
                        seq(from = outer_track_features$start[i],
                            to = outer_track_features$end[i]), ]
                circos.lines(x = this_bit_of_coverage$start,
                    y = this_bit_of_coverage$coverage,
                    sector.index = outer_track_features$chr[i],
                    area = TRUE, track.index = CELL_META$track.index,
                    col = outer_track_features$colour[i],
                    border = outer_track_features$colour[i])
                ## if it's in a normal rectangle
            } else {
                ## feature should end at the top of the
                ## ideogram: 1/2 of the way up the track
                if (internal) {
                    y_top <- (CELL_META$ylim[2])/2
                } else {
                    y_top <- CELL_META$ylim[2]
                }
                circos.rect(xleft = outer_track_features$start[[i]],
                    xright = outer_track_features$end[[i]],
                    col = outer_track_features$colour[[i]],
                    border = outer_track_features$colour[[i]],
                    ybottom = CELL_META$ylim[1],
                    ytop = y_top, sector.index = outer_track_features$chr[[i]],
                    track.index = CELL_META$track.index)
            }
        }
    }

    ## now for the rest of the features we need to
    ## know what the track outside this one is so we
    ## can orientate ourselves on the plot
    outer_index <- CELL_META$track.index

    ## for insertionDiagram, the track ylim is
    ## dependent on coverage because of the layout
    ## (ideogram and features on the same track).
    if (internal == "insertionDiagram" & !is.null(coverage_data)) {
        track_ylim <- c(0, max(coverage_data$coverage))
    } else {
        track_ylim <- c(0, 1)
    }

    ## plot track by track
    for (i in seq_along(unique(feature_data$track))) {
        this_track_features <- subset(feature_data,
            feature_data$track == i)
        circos.track(ylim = track_ylim, bg.border = NA,
            bg.col = NA, track.height = track_height)

        ## plot each insert on this track
        for (j in seq_along(this_track_features$chr)) {
            if (this_track_features$chr[j] %in%
                get.all.sector.index()) {
                ## adjust the coordinates if the sector is
                ## flipped
                if (this_track_features$chr[j] %in% flipped_sector) {
                    x_start <- .reverseXaxis(this_track_features$end[j])
                    x_end <- .reverseXaxis(this_track_features$start[j])
                    arrow_position <- ifelse(this_track_features$shape[j] ==
                        "forward_arrow", "start", "end")

                } else {
                    x_start <- this_track_features$start[j]
                    x_end <- this_track_features$end[j]
                    arrow_position <- ifelse(this_track_features$shape[j] ==
                        "forward_arrow", "end", "start")
                }

                ## arrows
                if (this_track_features$shape[j] ==
                    "forward_arrow" | this_track_features$shape[j] ==
                    "reverse_arrow") {
                    circos.arrow(x_start, x_end,
                        arrow.head.width = CELL_META$ylim[2] *
                            0.8, col = this_track_features$colour[j],
                        border = ifelse(feature_outline, "black",
                                        this_track_features$colour[j]),
                        arrow.position = arrow_position,
                        sector.index = this_track_features$chr[j])

                  ## rectangles
                } else if (this_track_features$shape[j] == "rectangle") {
                  circos.rect(
                    ybottom = (CELL_META$ylim[1] + (0.27 * CELL_META$ylim[2])),
                    xleft = x_start,
                    border = ifelse(feature_outline, "black",
                                        this_track_features$colour[j]),
                    ytop = (CELL_META$ylim[1] + (0.75 * CELL_META$ylim[2])),
                    xright = x_end, col = this_track_features$colour[j],
                    sector.index = this_track_features$chr[j])

                  ## triangles
                } else if (this_track_features$shape[j] ==
                    "upwards_triangle" | this_track_features$shape[j] ==
                    "downwards_triangle") {
                  x_midpoint <- (this_track_features$start[j] +
                    this_track_features$end[j])/2
                  y_midpoint <- mean(CELL_META$ylim)

                  if (this_track_features$shape[j] == "upwards_triangle") {
                      triangle_apex_y <- CELL_META$ylim[2] * 0.75

                  } else if (this_track_features$shape[j] ==
                      "downwards_triangle") {
                      triangle_apex_y <- CELL_META$ylim[1] * 0.75
                  }

                  circos.polygon(x = c(x_start, x_midpoint, x_end, x_start),
                      y = c(y_midpoint, triangle_apex_y,
                          y_midpoint, y_midpoint),
                      col = this_track_features$colour[[j]],
                      border = ifelse(feature_outline, "black",
                                        this_track_features$colour[j]),
                      sector.index = this_track_features$chr[[j]])
                } else {
                  stop("the 'shape' column should be one of 'rectangle',
                        'forward_arrow', 'reverse_arrow', 'upwards_triangle'
                        or 'downwards_triangle'. Please see ?featureDiagram
                        for more information")
                }
            } else {
                warning("The feature ", this_track_features$label[j],
                  " could not
                    be plotted because there is no sector matching ",
                  this_track_features$chr[j], ". Please check for typos &
                    that the names are an exact match between the feature and
                    ideogram data")
            }
        }
    }

    ## add the labels (do this last so we don't
    ## create unecessary tracks)
    for (i in seq_along(unique(feature_data$track))) {
        this_track_features <- subset(feature_data,
            feature_data$track == i)
        for (j in seq_along(this_track_features$label)) {
            if (this_track_features$label[j] != "") {
                if (this_track_features$chr[j] %in% get.all.sector.index()) {
                    ## adjust position if the sector has been flipped
                    if (as.character(this_track_features$chr[j]) %in%
                        flipped_sector) {
                        x_start <- .reverseXaxis(this_track_features$end[j])
                        x_end <- .reverseXaxis(this_track_features$start[j])

                    } else {
                        x_start <- this_track_features$start[j]
                        x_end <- this_track_features$end[j]
                    }

                    ## plot the text
                    this_track <- outer_index + this_track_features$track[1]
                    if (this_track_features$end[j] -
                        this_track_features$start[j] <= feature_label_cutoff) {
                        ## if the shape is small, plot the label on the
                        ## next track
                        .plotNewTrackText(
                            label = this_track_features$label[[j]],
                            x_start = x_start, x_end = x_end,
                            this_track = this_track,
                            sector_index = this_track_features$chr[[j]],
                            feature_label_size = feature_label_size,
                            label_track_height = label_track_height,
                            max_track = max(feature_data$track) + outer_index)
                    } else {
                        ## otherwise, the label should go inside the shape
                        circos.text(
                            x = ((x_start + x_end)/2),
                            y = mean(CELL_META$ylim),
                            labels = as.character(
                                this_track_features$label[j]),
                            facing = "bending.outside",
                            cex = feature_label_size,
                            track.index = this_track, niceFacing = TRUE,
                            sector.index = this_track_features$chr[[j]])
                    }
                }
            }
        }
    }
}

#' @title Add a legend
#'
#' @description Makes a legend object using ComplexHeatmap package which can
#' then be plotted using the \code{\link{gmovizPlot}} function.
#'
#' @param label_legend Whether to make a legend for labels (good for
#' colour-coded labels).
#' @param label_data The label data.
#' @param label_legend_title Title for the label legend.
#' @param feature_legend Whether to make a legend for features.
#' @param feature_data The feature data to use for the feature legend.
#' @param feature_legend_title Title for the features legend.
#' @param scatterplot_legend Whether to make a legend for the scatterplot
#' track.
#' @param scatterplot_legend_title Title for scatterplot track legend.
#' @param scatterplot_legend_labels A vector of the name/description of each
#' point e.g. if a point represents methylation, use 'methylation'. If we have
#' red/blue points for copy number gain/loss use c('gain', 'loss').
#' @param point_type,point_colour,point_outline_colour The type and colour of
#' points, as supplied to the \code{\link{drawScatterplotTrack}} function.
#' @param linegraph_legend Whether to plot a legend for a line graph track.
#' @param linegraph_legend_labels A vector of label(s) for what the line graph
#' means (e.g. \code{'Per Base Coverage'} for a line graph track showing
#' coverage).
#' @param linegraph_legend_colours The colour of to the line graph track.
#' @param linegraph_legend_title A title for the line graph legend.
#' @param background_colour The colour of the background (either 'white' or
#' 'black').
#'
#' @export
#' @importFrom grid gpar
#' @importFrom ComplexHeatmap packLegend
#' @importFrom ComplexHeatmap Legend
#' @importFrom GenomicFeatures features
#'
#' @return An object of the Legends class.
#' @seealso If you want more customisation over your legends, please see
#' \url{https://jokergoo.github.io/circlize_book/book/legends.html} for a
#' detailed guide as to how to implement legends alongside the circlize plots.
#' To plot these legends, see \code{\link{gmovizPlot}}
#'
#' @examples
#' ## a gene label legend
#' ## the data
#' labels <- data.frame(chr=c('chr1', 'chr1'), start=c(100, 300),
#' end=c(150, 350), label=c('a', 'b'), type=c('gene', 'lncRNA'),
#' colour=c('red', 'blue'))
#'
#' ## making the legend
#' makeLegends(label_legend=TRUE, label_data=labels)

makeLegends <- function(label_legend = FALSE, label_data = NULL,
    label_legend_title = "Gene Labels", feature_legend = FALSE,
    feature_data = NULL, feature_legend_title = "Features",
    scatterplot_legend = FALSE, scatterplot_legend_labels = c("Gains",
        "Losses"), point_colour = "black", point_outline_colour = "black",
    point_type = 21, scatterplot_legend_title = "Copy Number Variants",
    linegraph_legend = FALSE, linegraph_legend_labels = "Per Base Coverage",
    linegraph_legend_colours = "black", linegraph_legend_title = "Line Graph",
    background_colour = "white") {
    ## check inputs and packages required
    stopifnot(exprs = {
        methods::is(label_legend, "logical")
        methods::is(label_legend_title, "character")
        methods::is(feature_legend, "logical")
        methods::is(feature_legend_title, "character")
        methods::is(scatterplot_legend, "logical")
        methods::is(scatterplot_legend_labels,
            "character")
        methods::is(point_colour, "character")
        methods::is(scatterplot_legend_title, "character")
        methods::is(linegraph_legend, "logical")
        methods::is(linegraph_legend_labels, "character")
        methods::is(linegraph_legend_colours, "character")
        methods::is(linegraph_legend_title, "character")
        background_colour %in% c("black", "white")
    })

    if (!requireNamespace(c("ComplexHeatmap", "grid"),
        quietly = TRUE)) {
        stop("The packages 'ComplexHeatmap' and 'grid' are needed for legend
            functionality. Please install them.",
            call. = FALSE)
    }

    legend_to_plot <- list()

    ## if bg is white, text should be black & vice
    ## versa
    text_colour <- ifelse(background_colour ==
        "white", "black", "white")

    if (label_legend == TRUE) {
        label_data <- .checkLabelData(label_data)

        if (any(!c("type", "colour") %in% colnames(label_data))) {
            ## non colour coded
            label_types <- "Labels"
            label_colours <- text_colour

        } else {
            ## colour coded
            label_types <- unique(label_data$type)
            label_colours <- unique(label_data$colour)
        }
        label_legend <- ComplexHeatmap::Legend(at = label_types,
            border = background_colour,
            labels_gp = grid::gpar(col = text_colour),
            title = label_legend_title,
            legend_gp = grid::gpar(fill = label_colours),
            title_position = "topleft",
            title_gp = grid::gpar(col = text_colour, font = 2))
        legend_to_plot <- append(legend_to_plot,
            list(label_legend))
    }

    if (feature_legend == TRUE) {
        feature_data <- .checkFeatureData(feature_data)
        ## if there is a `type` column then use that for
        ## the legend
        if (!is.null(feature_data$type)) {
            categories <- unique(feature_data$type)
            colours <- unique(feature_data$colour)
        } else {
            categories <- unique(feature_data$label)
            colours <- vector(length = length(categories))
            for (i in seq_along(categories)) {
                colours[i] <- (feature_data[feature_data$label ==
                  categories[i], ])$colour
            }
        }
        feature_legend <- ComplexHeatmap::Legend(at = categories,
            border = background_colour,
            labels_gp = grid::gpar(col = text_colour),
            title_position = "topleft",
            legend_gp = grid::gpar(fill = colours),
            title = feature_legend_title,
            title_gp = grid::gpar(col = text_colour, font = 2))
        legend_to_plot <- append(legend_to_plot,
            list(feature_legend))
    }

    if (scatterplot_legend == TRUE) {
        scatterplot_legend <- ComplexHeatmap::Legend(
            at = scatterplot_legend_labels, background = background_colour,
            labels_gp = grid::gpar(col = text_colour),
            border = background_colour, legend_gp = grid::gpar(
                col = point_outline_colour, fill = point_colour),
            type = "points", title = scatterplot_legend_title,
            title_position = "topleft", pch = point_type,
            title_gp = grid::gpar(col = text_colour, font = 2))
        legend_to_plot <- append(legend_to_plot, list(scatterplot_legend))
    }

    if (linegraph_legend == TRUE) {
        linegraph_legend <- ComplexHeatmap::Legend(
            at = linegraph_legend_labels,
            labels_gp = grid::gpar(col = text_colour),
            title_position = "topleft",
            legend_gp = grid::gpar(fill = linegraph_legend_colours),
            border = background_colour, title = linegraph_legend_title,
            title_gp = grid::gpar(col = text_colour, font = 2))
        legend_to_plot <- append(legend_to_plot,
            list(linegraph_legend))
    }

    if (length(legend_to_plot) == 0) {
        return(FALSE)

    } else {
        legend_to_plot <- ComplexHeatmap::packLegend(list = legend_to_plot)
        return(legend_to_plot)
    }
}
#' @title Display number of copies of an insertion
#'
#' @description Generates a diagram which displays insertions, showing their
#' position, size and copy number. See \code{\link{featureDiagram}} for
#' a more general function which can display other features of interest.
#'
#' @inheritParams gmovizInitialise
#' @param insertion_data A \code{\link{GRanges}} or data frame describing the
#' insertion. See below for the detailed format.
#' @param style How the original sequence and insertions will be positioned
#' around the circle. Choose from options 1, 2, 3 or 4. Please see the examples
#' below or the vignette for what these options represent.
#' @param either_side How much extra of the genome should be shown around the
#' insertion site. Can be either a single number (e.g. \code{1000}, then 1000bp
#' will be shown either side of the insertion site), a vector of length 2
#' (e.g. \code{c(2000, 13000)} in which case from 2000 to 13000 will be shown)
#' or a GRanges (in which case all ranges in the GRanges object will be used
#' to determine the start/end points of the sector)
#' @param insertion_label The label(s) that will be applied to the insertions.
#' If \code{'default'} then the name of the insertion will be used to label
#' single copy insertions and a number will be used for multiple copy number
#' insertions. Otherwise, \code{insertion_label} should be a vector with one
#' element for each row of the insertion data, indicating the label that should
#' be used for that insertion.
#' @param link_colour The colour of the link: this should be a 6 digit hex
#' code, the transparency is automatically added.
#' @param link_ends How far the link extends in either direction. \emph{This is
#' set automatically} but if you want to edit it, provide a vector of length 2
#' with each element being between 0 (centre of circle) and 1 (right at the
#' edge of the circle).
#' @param internal For internal use only.
#'
#' @export
#' @import circlize
#' @importFrom graphics par
#' @importFrom graphics strwidth
#' @importFrom methods is
#'
#' @section Insertion data format:
#' The start, end and seqnames of insertion_data \linkS4class{GRanges} should
#' describe the insertion site. Additionally, there are five metadata columns:
#' \describe{
#' \item{name}{A character string which will be used to label insertion. It
#' is suggested to keep this label relatively short, if possible.}
#' \item{colour}{A character string of a colour to use. Supports hex colours
#' (\emph{e.g. #000000}) and named R colours (\emph{e.g. red}).}
#' \item{shape}{The shape that will be used to represent the feature: \itemize{
#' \item{\code{'rectangle'}} is a rectangle.
#' \item{\code{'forward_arrow'}} for a forwards facing arrow.
#' \item{\code{'reverse_arrow'}} for a backwards (reverse) facing arrow.} It is
#' suggested to use \code{'forward_arrow'}}
#' \item{length}{The length of the insertion}
#' \item{in_tandem}{The number of copies of the insert in tandem}}
#' The columns \strong{in_tandem, colour and shape are all optional}. If you
#' don't supply them, then default values will be added as follows: \describe{
#' \item{in_tandem}{1 (only one copy inserted)}
#' \item{colour}{a colour allocated from \code{\link{rich_colours}}}
#' \item{shape}{\code{'forward_arrow'}}}
#'
#' @section Warning:
#' If you choose to use a data frame to supply the insertion_data, please be
#' careful to add the \code{stringsAsFactors=FALSE} argument. Otherwise, the
#' colours may not be correct.
#'
#' @return Generates an image displaying the copy number of the insertion(s)
#' provided
#' @seealso \code{\link{featureDiagram}} for a more flexible function
#' that takes a similar approach to representing features of interest.
#'
#' @examples
#' ## one insertion with 4 tandem copies
#' ## the data as a data.frame
#' exampleins <- data.frame(
#' chr='chr12', start=70905597, end=70917885, name='plasmid',
#' colour='#7270ea', length=12000, in_tandem=11, shape='forward_arrow',
#' stringsAsFactors=FALSE)
#'
#' ## or we can supply it as GRanges (same thing)
#' exampleins <- GRanges(
#' seqnames='chr12', ranges=IRanges(start=70905597, end=70917885),
#' name='plasmid', colour='#7270ea', length=12000, in_tandem=11,
#' shape='forward_arrow')
#'
#' ## plot it
#' insertionDiagram(exampleins, either_side=c(70855503, 71398284))
#'
#' ## that was the default 'style'. The other 3 styles are:
#' ## style 2
#' insertionDiagram(exampleins, either_side=c(70855503, 71398284), style=2)
#'
#' ## style 3
#' insertionDiagram(exampleins, either_side=c(70855503, 71398284), style=3)
#'
#' ## style 4
#' insertionDiagram(exampleins, either_side=c(70855503, 71398284), style=4)
#'
#' ## 2 different insertions
#' ## the data
#' example2ins <- data.frame(
#' chr=c('chr12', 'chr12'), start=c(70905597, 70705597),
#' end=c(70917885, 70717885), name=c('plasmid1', 'plasmid2'),
#' colour=c('#7270ea', '#ea7082'), length=c(12000, 10000),
#' in_tandem=c(4, 8), shape=c('reverse_arrow', 'forward_arrow'),
#' stringsAsFactors=FALSE)
#'
#' ## plot it
#' insertionDiagram(example2ins, link_colour='#ffe677', start_degree=45)

insertionDiagram <- function(insertion_data, style = 1,
    either_side = "default", insertion_label = "default",
    sector_colours = nice_colours, sector_border_colours = nice_colours,
    start_degree = 180, custom_sector_width = NULL,
    coverage_rectangle = NULL, coverage_data = NULL,
    custom_ylim = NULL, space_between_sectors = 15,
    sector_labels = TRUE, sector_label_size = 1.3,
    sector_label_colour = "black", label_data = NULL,
    label_colour = "black", link_colour = "default",
    label_size = 1.1, xaxis = TRUE, xaxis_label_size = 0.9,
    xaxis_colour = "#747577", xaxis_spacing = 10, xaxis_spacing_unit = "deg",
    link_ends = "default", track_height = 0.15,
    internal = FALSE) {
    ## check inputs
    insertion_data <- .checkInsertionData(insertion_data)
    if (!is.null(coverage_data)) {
        coverage_data <- .checkCoverageData(coverage_data)
    }
    stopifnot(exprs = {
        insertion_label == "default" | length(insertion_label) ==
            nrow(insertion_data)
        style %in% c(1, 2, 3, 4)
        methods::is(sector_colours, "character")
        methods::is(sector_border_colours, "character")
        start_degree >= 0 & start_degree <= 360
        is.null(custom_sector_width) | length(custom_sector_width) ==
            length(unique(insertion_data$chr)) +
                1
        is.null(coverage_rectangle) | methods::is(coverage_rectangle,
            "character")
        is.null(custom_ylim) | length(custom_ylim ==
            2)
        space_between_sectors >= 0
        methods::is(sector_labels, "logical")
        sector_label_size >= 0
        methods::is(label_colour, "character")
        label_size >= 0
        (methods::is(xaxis_spacing, "numeric") &
            xaxis_spacing > 0 | xaxis_spacing == "start_end")
        xaxis_spacing_unit %in% c("deg", "bp")
        link_ends == "default" | (length(link_ends) ==
            2 & methods::is(link_ends, "numeric"))
        track_height > 0 & track_height < 1
    })

    ## convert the insertions data to feature data
    ## to reuse drawFeatureTrack
    feature_data <- .insertionsToFeatures(insertion_data,
        insertion_label)

    ## make the ideogram_data data from the
    ## insertion data
    ideogram_data <- .insertionsToIdeogram(insertion_data,
        either_side)
    original_sequence <- as.character(insertion_data$chr[1])
    inserted_sequence <- insertion_data$name
    inserted_ideogram <- subset(ideogram_data,
        ideogram_data$chr %in% inserted_sequence)

    ## the link data: this is just based off the
    ## start and end of the integration event
    link_side_1 <- data.frame(chr = insertion_data$chr,
        start = insertion_data$start, end = insertion_data$end,
        stringsAsFactors = FALSE)
    link_side_2 <- data.frame(chr = inserted_sequence,
        start = inserted_ideogram$start, end = inserted_ideogram$end,
        stringsAsFactors = FALSE)

    ## sector width is based on number of sectors to
    ## be displayed
    number_of_sectors <- nrow(ideogram_data)
    if (is.null(custom_sector_width)) {
        if (number_of_sectors == 2) {
            sector_width <- c(0.35, 0.65)
        } else if (number_of_sectors %in% 3:5) {
            site_width <- (1 - 0.3)/number_of_sectors
            sector_width <- c(0.3, rep(site_width,
                number_of_sectors - 1))
        } else {
            site_width <- (1 - 0.2)/number_of_sectors
            sector_width <- c(0.2, rep(site_width,
                number_of_sectors - 1))
        }
    } else {
        ## or can use custom values
        sector_width <- custom_sector_width
    }

    ## initialise
    circos.clear()
    par(xpd = NA)
    circos.par(start.degree = start_degree, gap.after = space_between_sectors,
        unit.circle.segments = 1000, track.margin = c(0,
            0), cell.padding = c(0, 0, 0, 0))

    circos.initializeWithIdeogram(cytoband = ideogram_data,
        plotType = NULL, sector.width = sector_width)

    ## gene labels, if needed
    if (!is.null(label_data)) {
        label_data <- .checkLabelData(label_data)
        if (internal) {
            # if using multipleInsertionDiagram, make v.
            # small
            label_height <- convert_height(0.01,
                "cm")
            connection_height <- convert_height(1,
                "mm")
        } else {
            # these are the default values
            label_height <- min(c(convert_height(1.5,
                "cm"), max(strwidth(label_data[,
                4], cex = label_size, font = par("font")))))
            connection_height <- convert_height(5,
                "mm")
        }

        circos.genomicLabels(bed = label_data,
            labels.column = 4, col = label_colour,
            cex = label_size, line_col = label_colour,
            side = "outside", connection_height = connection_height,
            labels_height = label_height, track.margin = c(0,
                0))

    }

    ## plot the links first so features can be drawn
    ## on top. We need to know where they start and
    ## end (determined by the style of the diagram)
    ## give option for manual override make
    ## allowance for label track
    if (link_ends == "default") {
        if (!is.null(label_data)) {
            label_track_height <- label_height +
                connection_height
        } else {
            label_track_height <- 0
        }

        ## now depending on the style of plot: (0.01 is
        ## the margin)
        if (style == 1) {
            insertion_bit <- 1 - label_track_height -
                track_height - 0.01
            original_bit <- 1 - label_track_height -
                (2 * track_height) - 0.01
        } else if (style == 2) {
            insertion_bit <- 1 - label_track_height -
                2 * (track_height + 0.01)
            original_bit <- 1 - label_track_height -
                track_height - 0.01
        } else if (style == 3) {
            insertion_bit <- 1 - label_track_height -
                track_height - 0.02
            original_bit <- 1 - label_track_height -
                (2 * track_height) - 0.01
            link_side_2$start <- link_side_2$start +
                (0.4 * insertion_data$length)
            link_side_2$end <- link_side_2$end -
                (0.4 * insertion_data$length)
        } else if (style == 4) {
            insertion_bit <- 1 - label_track_height -
                track_height - 0.01
            original_bit <- 1 - label_track_height -
                track_height - 0.01
            link_side_2$start <- link_side_2$start +
                (0.4 * insertion_data$length)
            link_side_2$end <- link_side_2$end -
                (0.4 * insertion_data$length)
        }
        link_ends <- c(original_bit, insertion_bit)
    }
    if (link_colour == "default") {
        set_link_colour <- sector_colours[seq(1,
            nrow(insertion_data))]
    } else {
        set_link_colour <- link_colour
    }
    circos.genomicLink(link_side_1, link_side_2,
        border = set_link_colour, rou1 = link_ends[1],
        col = paste0(set_link_colour, "33"), rou2 = link_ends[2])

    ## the order we plot things in depends on the
    ## style (plotting first = on the outside of the
    ## circle)
    if (style == 1) {
        ## outer box for insertion(s)
        .plotNewTrack(inserted_sequence, custom_ylim,
            coverage_rectangle, coverage_data,
            track_height)
        .plotIdeogram(inserted_sequence, coverage_data,
            coverage_rectangle, get.cell.meta.data("track.index"),
            sector_colours[-1], sector_border_colours[-1],
            sector_labels, sector_label_size, sector_label_colour)
        .plotXaxis(xaxis, xaxis_label_size, xaxis_colour,
            xaxis_spacing, sectors = inserted_sequence,
            xaxis_spacing_unit = xaxis_spacing_unit)

        ## insertions (features)
        drawFeatureTrack(feature_data, feature_label_cutoff = 5,
            track_height = track_height, internal = "insertionDiagram",
            coverage_data = coverage_data)

        # inner box for original sequence
        .plotIdeogram(original_sequence, coverage_data,
            coverage_rectangle, get.cell.meta.data("track.index"),
            sector_colours[1], sector_border_colours[1],
            sector_labels, sector_label_size, sector_label_colour)
        .plotXaxis(xaxis, xaxis_label_size, xaxis_colour,
            xaxis_spacing, sectors = inserted_sequence,
            xaxis_spacing_unit = xaxis_spacing_unit)

    } else if (style == 2) {
        ## outer box for original sequence
        .plotNewTrack(original_sequence, custom_ylim,
            coverage_rectangle, coverage_data,
            track_height)
        .plotIdeogram(original_sequence, coverage_data,
            coverage_rectangle, get.cell.meta.data("track.index"),
            sector_colours[1], sector_border_colours[1],
            sector_labels, sector_label_size, sector_label_colour)
        .plotXaxis(xaxis, xaxis_label_size, xaxis_colour,
            xaxis_spacing, sectors = inserted_sequence,
            xaxis_spacing_unit = xaxis_spacing_unit)

        ## inner box for insertion(s)
        .plotNewTrack(inserted_sequence, custom_ylim,
            coverage_rectangle, coverage_data,
            track_height)
        .plotIdeogram(inserted_sequence, coverage_data,
            coverage_rectangle, get.cell.meta.data("track.index"),
            sector_colours[-1], sector_border_colours[-1],
            sector_labels, sector_label_size, sector_label_colour)
        .plotXaxis(xaxis, xaxis_label_size, xaxis_colour,
            xaxis_spacing, sectors = inserted_sequence,
            xaxis_spacing_unit = xaxis_spacing_unit)

        ## insertions (features)
        drawFeatureTrack(feature_data, feature_label_cutoff = 5,
            track_height = track_height, internal = "insertionDiagram",
            coverage_data = coverage_data)

    } else if (style == 3) {
        ## insertions (features)
        .plotNewTrack(current_sectors = inserted_sequence,
            custom_ylim, coverage_rectangle, coverage_data,
            track_height)
        drawFeatureTrack(feature_data, feature_label_cutoff = 5,
            track_height = track_height, internal = "insertionDiagram",
            coverage_data = coverage_data)

        ## inner box for original sequence
        .plotNewTrack(current_sectors = original_sequence,
            custom_ylim, coverage_rectangle, coverage_data,
            track_height)
        .plotIdeogram(original_sequence, coverage_data,
            coverage_rectangle, get.cell.meta.data("track.index"),
            sector_colours[1], sector_border_colours[1],
            sector_labels, sector_label_size, sector_label_colour)
        .plotXaxis(xaxis, xaxis_label_size, xaxis_colour,
            xaxis_spacing, sectors = inserted_sequence,
            xaxis_spacing_unit = xaxis_spacing_unit)

    } else if (style == 4) {
        ## outer box for original sequence
        .plotNewTrack(original_sequence, custom_ylim,
            coverage_rectangle, coverage_data,
            track_height)
        .plotIdeogram(original_sequence, coverage_data,
            coverage_rectangle, get.cell.meta.data("track.index"),
            sector_colours[1], sector_border_colours[1],
            sector_labels, sector_label_size, sector_label_colour)
        .plotXaxis(xaxis, xaxis_label_size, xaxis_colour,
            xaxis_spacing, sectors = inserted_sequence,
            xaxis_spacing_unit = xaxis_spacing_unit)

        ## insertions (features)
        drawFeatureTrack(feature_data, feature_label_cutoff = 5,
            track_height = track_height, internal = "insertionDiagram",
            coverage_data = coverage_data)
    }
}
#' @title Display 'features' of interest in a diagram
#'
#' @description Generates a diagram which displays 'features' (e.g. genes,
#' indels, primer sequences etc) using coloured shapes. See
#' \code{\link{insertionDiagram}} for a similar function which specialises in
#' plotting insertions or \code{\link{drawFeatureTrack}} to add a feature track
#' to an existing graph.
#'
#' @inheritParams gmovizInitialise
#' @inheritParams drawFeatureTrack
#' @param link_data If you would like to draw a link between two sectors of the
#' circle, \code{link_data} should be a data frame with two rows: one for each
#' end of the link. There should be 3 columns: \code{chr}, \code{start} &
#' \code{end} which describe the position of each end of the link.
#' @param link_colour The colour of the link: this should be a 6 digit hex
#' code, the transparency is automatically added.
#' @param link_ends How far the link extends in either direction. \emph{This is
#' set automatically} but if you want to edit it, provide a vector of length 2
#' with each element being between 0 (centre of circle) and 1 (right at the
#' edge of the circle).
#' @section Warning:
#' If you choose to use a data frame to supply the feature data, please be
#' careful to add the \code{stringsAsFactors = FALSE} argument. Otherwise, the
#' colours may not be correct.
#'
#' @export
#' @import circlize
#'
#' @return Generates an image of the feature data supplied.
#' @seealso \code{\link{insertionDiagram}} for a more specialised function
#' which shows the copy number of insertions. Also
#' \code{\link{drawFeatureTrack}} to add the exact same feature information to
#' an existing plot and \code{\link{getFeatures}} for a function that can
#' read in the feature information from a .gff file.
#'
#' @examples
#' ## plasmid map
#' plasmid_ideogram <- data.frame(chr='plasmid', start=0, end=2500)
#'
#' plasmid_features <- GRanges(seqnames=rep('plasmid', 4),
#' ranges=IRanges(start=c(0, 451, 901, 1700), end=c(450, 900, 1400, 2200)),
#' colour=c('#d44a9f', '#4a91d4', '#7ad44a', '#d49d4a'),
#' label=c('promoter', 'gene', 'GFP', 'ampR'),
#' shape=c('rectangle', 'forward_arrow', 'forward_arrow', 'reverse_arrow'),
#' track=rep(1, 4))
#'
#' featureDiagram(plasmid_ideogram, plasmid_features)

featureDiagram <- function(ideogram_data, feature_data,
    start_degree = 180, coverage_rectangle = NULL,
    coverage_data = NULL, custom_sector_width = NULL,
    space_between_sectors = 4, flipped_sector = NULL,
    sector_colours = nice_colours, 
    sector_border_colours = nice_colours,
    sector_labels = TRUE, sector_label_size = 1.3,
    sector_label_colour = "black", label_data = NULL,
    label_size = 1.1, label_colour = "black", xaxis = TRUE,
    xaxis_label_size = 0.9, xaxis_colour = "#747577",
    xaxis_spacing = 10, feature_label_cutoff = 50, xaxis_spacing_unit = "deg",
    track_height = 0.1, feature_label_size = 0.9,
    link_data = NULL, link_colour = "#84c6d6",
    link_ends = "default", custom_ylim = NULL,
    label_track_height = 0.1 * feature_label_size,
    feature_outline = TRUE) {
    ## check data
    feature_data <- .checkFeatureData(feature_data)
    ideogram_data <- .checkIdeogramData(ideogram_data)

    if (!is.null(coverage_data)) {
        coverage_data <- .checkCoverageData(coverage_data)
    }

    ## check other inputs:
    stopifnot(exprs = {
        start_degree >= 0 & start_degree <= 360
        is.null(coverage_rectangle) | methods::is(coverage_rectangle,
            "character")
        is.null(custom_sector_width) | length(custom_sector_width) ==
            nrow(ideogram_data)
        space_between_sectors >= 0
        is.null(flipped_sector) | all(flipped_sector %in%
            ideogram_data$chr)
        methods::is(sector_colours, "character")
        methods::is(sector_border_colours, "character")
        methods::is(sector_labels, "logical")
        sector_label_size >= 0
        label_size >= 0
        methods::is(label_colour, "character")
        xaxis_label_size >= 0
        (methods::is(xaxis_spacing, "numeric") &
            xaxis_spacing > 0 | xaxis_spacing == "start_end")
        xaxis_spacing_unit %in% c("deg", "bp")
        feature_label_cutoff > 0
        track_height > 0 & track_height < 1
        feature_label_size >= 0
        is.null(link_data) | methods::is(link_data,
            "data.frame")
        link_ends == "default" | (length(link_ends) ==
            2 & methods::is(link_ends, "numeric"))
        is.null(custom_ylim) | length(custom_ylim ==
            2)
        label_track_height > 0 & label_track_height <
            1
    })

    ## initialise
    circos.clear()
    par(xpd = NA)
    circos.par(start.degree = start_degree, gap.after = space_between_sectors,
        unit.circle.segments = 1000, track.margin = c(0,
            0), cell.padding = c(0, 0, 0, 0))
    circos.initializeWithIdeogram(cytoband = ideogram_data,
        plotType = NULL, sector.width = custom_sector_width)

    ## let people know if they're plotting into a
    ## sector that doesn't exist
    .checkSectorsMatch(feature_data)

    ## labels, if needed
    if (!is.null(label_data)) {
        label_data <- .checkLabelData(label_data)
        circos.genomicLabels(bed = label_data,
            labels.column = 4, col = label_colour,
            cex = label_size, line_col = label_colour,
            side = "outside")
    }

    ## plot the ideogram
    .plotNewTrack(current_sectors = ideogram_data$chr,
        custom_ylim, coverage_rectangle, coverage_data,
        track_height = 0.1)
    .plotIdeogram(current_sectors = get.all.sector.index(),
        coverage_data = coverage_data, sector_labels = sector_labels,
        coverage_rectangle = coverage_rectangle,
        border_colour = sector_border_colours,
        sector_label_colour = sector_label_colour,
        current_track = get.cell.meta.data("track.index"),
        colour = sector_colours, sector_label_size = sector_label_size)

    if (!is.null(xaxis_colour)) {
        .plotXaxis(xaxis, xaxis_label_size = xaxis_label_size,
            xaxis_colour = xaxis_colour, xaxis_spacing = xaxis_spacing,
            flipped_sector = flipped_sector,
            sectors = as.character(ideogram_data$chr),
            xaxis_spacing_unit = xaxis_spacing_unit)
    }

    ## add link if needed (needs to be done now so
    ## the features can go on top)
    if (!is.null(link_data)) {
        if (link_ends == "default") {
            # option for manual override allowance for gene
            # label track
            if (!is.null(label_data)) {
                label_track_height <- min(c(convert_height(1.5,
                  "cm"), max(strwidth(label_data[,
                  4], cex = label_size, font = par("font"))))) +
                  convert_height(5, "mm")
            } else {
                label_track_height <- 0
            }
            inside_of_ideogram <- 1 - label_track_height -
                0.1 - 0.02
            link_ends <- c(inside_of_ideogram,
                inside_of_ideogram)
        }

        circos.genomicLink(link_data[1, ], link_data[2,
            ], border = link_colour, rou1 = link_ends[1],
            rou2 = link_ends[2], col = paste0(link_colour,
                "33"))
    }

    ## draw the actual features
    drawFeatureTrack(feature_data = feature_data,
        flipped_sector = flipped_sector,
        feature_label_cutoff = feature_label_cutoff,
        track_height = track_height, coverage_data = coverage_data,
        feature_label_size = feature_label_size,
        label_track_height = label_track_height,
        coverage_rectangle = coverage_rectangle,
        internal = TRUE, feature_outline = feature_outline)
}
#' @title Generate an entire circular plot
#'
#' @description Saves code supplied to \code{plotting_functions} a plot (with
#' optional title and legends) as either .png, .svg or .ps.
#'
#' @param file_name The name of the file to be saved.
#' @param file_type The type of image file to produce: either \code{'png'},
#' \code{'svg'} or \code{'ps'}.
#' @param plotting_functions The functions you want to plot (e.g.
#' \code{\link{insertionDiagram}} or \code{\link{gmovizInitialise}}).
#' @param legends A legend object to plot, generated by
#' \code{\link{makeLegends}}.
#' @param title Text for the title, leave as \code{NULL} for no title.
#' @param width Width of the image.
#' @param height Height of the image.
#' @param units Units for the width and height of the image. One of
#' \code{'mm'}, \code{'cm'} or \code{'in'} (inches).
#' @param res Resolution of the image (only needed for .png files).
#' @param background_colour Colour of the image background.
#' @param title_x_position,title_y_position X and Y positions of the title on
#' the image.
#' @param title_font_face Font face of the title: bold, italic or bold-italic.
#' @param title_size Size of the title.
#' @param title_colour Colour of the title.
#' @param point_size Pointsize (for postscript output only).
#'
#' @export
#' @importFrom grid unit
#' @importFrom grid pushViewport
#' @importFrom grid viewport
#' @importFrom gridBase gridOMI
#' @importFrom grid upViewport
#' @importFrom ComplexHeatmap draw
#' @importFrom graphics grconvertY
#' @importFrom graphics grconvertX
#'
#' @return Saves a plot to disk in the specified format.
#' @seealso \code{\link{makeLegends}} for a function that generates the legend
#' objects.
#'
#' @examples
#' ## make some example data
#' small_ideogram <- data.frame(chr=c('region 1', 'region 2', 'region 3'),
#' start=c(0, 0, 0), end=c(10000, 12000, 10000))
#' small_plot_data <- data.frame(
#' chr=sample(c('region 1', 'region 2', 'region 3'), size=40, replace=TRUE),
#' start=sample(0:10000, 40), end=sample(0:10000, 40),
#' val=rnorm(40, 2, 0.5))
#'
#' ## plot it
#' \dontrun{
#' gmovizPlot('test.png', {
#' gmovizInitialise(small_ideogram, custom_sector_width=c(0.3, 0.3, 0.3))
#' drawScatterplotTrack(small_plot_data)}, title='scatterplot')}

# test of the new function to handle plotting
gmovizPlot <- function(file_name, file_type = "png",
    plotting_functions, legends = NULL, title = NULL,
    width = 338.7, height = 238.7, units = "mm",
    res = 300, background_colour = "transparent",
    title_x_position = 0.5, title_y_position = 0.9,
    title_font_face = "bold", title_size = 1.1,
    title_colour = "black", point_size = 11) {
    ## check we have the right packages and inputs
    if (!requireNamespace(c("grid", "gridBase"),
        quietly = TRUE)) {
        stop("The packages 'grid' and 'gridBase are needed for this function.
            Please install them.",
            call. = FALSE)
    }

    stopifnot(exprs = {
        file_type %in% c("png", "svg", "ps")
        is.null(legends) | methods::is(legends,
            "Legends")
        is.null(title) | !is.na(as.character(title))
        width > 0
        height > 0
        units %in% c("mm", "cm", "in")
        res > 0
        title_x_position >= 0 & title_x_position <=
            1
        title_y_position >= 0 & title_y_position <=
            1
        title_font_face %in% c("", "bold", "italic",
            "bold-italic")
        !is.na(as.numeric(title_size))
        !is.na(as.numeric(point_size))
    })

    ## choose either png or svg or ps devices
    if (file_type == "png") {
        grDevices::png(filename = file_name, width = width,
            height = height, units = units, res = res,
            bg = background_colour)

    } else {
        ## svg and ps only support inches so convert
        if (units == "mm") {
            w <- width/25.4
            h <- height/25.4
        } else if (units == "cm") {
            w <- width/2.54
            h <- height/2.54
        } else {
            w <- width
            h <- height
        }

        if (file_type == "svg") {
            grDevices::svg(filename = file_name,
                width = w, height = h, bg = background_colour)

        } else if (file_type == "ps") {
            grDevices::postscript(file = file_name,
                bg = background_colour, paper = "special",
                width = w, height = h, pointsize = point_size)
        }
    }

    ## plot the legend and the graph:
    graphics::plot.new()
    circle_size <- grid::unit(1, "snpc")
    grid::pushViewport(grid::viewport(x = 0, y = 0.5,
        width = circle_size, height = circle_size,
        just = c("left", "center")))
    graphics::par(omi = gridBase::gridOMI(), new = TRUE)

    plotting_functions

    grid::upViewport()
    if (!is.null(legends)) {
        ComplexHeatmap::draw(legends, x = circle_size,
            just = "left")
    }

    if (!is.null(title)) {
        .plotTitle(title, title_x_position, title_y_position,
            title_font_face, title_size, title_colour)
    }
    grDevices::dev.off()
}

#' @title Display multiple insertion events around a genome
#'
#' @description Generates a diagram which displays multiple insertion events
#' (as displayed using the \code{\link{insertionDiagram}} function) around a
#' central genome
#'
#' @inheritParams gmovizInitialise
#' @param insertion_data A \code{\link{GRanges}} or data frame describing each
#' of the insertion events. Please see \code{\link{insertionDiagram}} for a
#' detailed description of the format.
#' @param genome_ideogram_data Either a \code{\link{GRanges}} representing
#' regions of interest or a data frame in bed format (containing the chr,
#' start and end columns). If you want to read in data from file, please see
#' the \code{\link{getIdeogramData}} function.
#' @param either_side How much extra of the genome should be shown around the
#' insertion site. See \code{\link{insertionDiagram}} for a description of the
#' different ways you can specify \code{either_side}, but note that for this
#' function you need to supply either one value (which will apply to all of
#' the events) or a named list of values (with one element per event. The names
#' should be the names of the insertion _NOT_ the names of the chromosomes).
#' @param track_height The height (vertical distance around the circle) that
#' will be taken up by this track. Should be a number between 0 (none) and 1
#' (entire circle) that will apply to all of the events.
#' @param style How the original sequence and insertions are positioned around
#' the circle (style 1, 2, 3 or 4). Please see the examples of the
#' \code{\link{insertionDiagram}} function or the vignette for what each of
#' these options look like. This should be either a single value (which will
#' apply to all of the events) or a named vector of values (with one element
#' per event).
#' @param colour_set The set of colours that will be used to create the
#' diagram. For simplicity, it isn't possible to specify precisely the colour
#' of each sector and link in the diagram (but you can easily edit them by
#' saving the diagram in a vectorised format and opening it in any vector
#' graphics editing program). See \code{\link{colourSets}} for the built-in
#' gmoviz colour sets or make your own (should be a vector of hex colours; must
#' have a length greater than or equal to the number of rows
#' of genome_ideogram_data)
#' @param xaxis_spacing Space between the x axis labels, in degrees.
#' Alternatively, the string 'start_end' will place a label at the start and
#' end of each sector only. Accepts only a single value which will be applied
#' to all events.
#'
#' @section Warning:
#' Due to space limitations, it isn't possible to display more than 8 events
#' or more than 3 events in the same quarter of the circle. If you have more
#' events than this, please consider splitting them across two or more figures.
#'
#'
#' @export
#' @import circlize
#' @importFrom grid grid.newpage
#' @importFrom grid pushViewport
#' @importFrom grid viewport
#' @importFrom grid grid.polygon
#' @importFrom grid gpar
#' @importFrom gridBase gridOMI
#' @importFrom colorspace lighten
#' @importFrom pracma deg2rad
#' @importFrom graphics grconvertY
#' @importFrom graphics grconvertX
#'
#' @return Generates an image of the multiple insertion events provided.
#' @seealso \code{\link{insertionDiagram}} for a function which generates the
#' figures for each of the individual events and \code{\link{gmovizInitialise}}
#' for the function which draws the central genome
#'
#' @examples
#' ## the data
#' ideogram_data <- GRanges(
#' seqnames=paste0('chr', 1:6), ranges=IRanges(start=rep(0, 6),
#'  end=rep(12000, 6)))
#' insertion_data <- GRanges(
#' seqnames = c('chr1', 'chr5'),
#' ranges = IRanges(start = c(4000, 2000), end = c(4100, 2200)),
#' name = c('ins1', 'ins5'), length = c(100, 200))
#'
#' ## the plot
#' multipleInsertionDiagram(insertion_data=insertion_data,
#'                          genome_ideogram_data=ideogram_data)
#' ## with coverage and labels
#' example_labels <- GRanges(seqnames=c('chr1', 'chr5'),
#'                           ranges=IRanges(start=c(4000, 2000),
#'                           end=c(4120, 2200)),
#'                           label=c('Gene A', 'Gene B'),
#'                           colour=c('red', 'blue'))
#'
#' example_coverage <- GRanges(
#' seqnames = c(rep('chr1', 100), rep('chr5', 100)),
#' ranges = IRanges(start=c(seq(4000, 4099, length.out=100),
#'                          seq(2000, 2199, length.out=100)),
#'                  end=c(seq(4001, 4100, length.out=100),
#'                        seq(2001, 2200, length.out=100))),
#'                  coverage=c(runif(100, 0, 25), runif(100, 0, 15)))
#' multipleInsertionDiagram(insertion_data=insertion_data,
#'                          genome_ideogram_data=ideogram_data,
#'                          label_data=example_labels,
#'                          label_colour=example_labels$colour,
#'                          coverage_rectangle=c('chr1', 'chr5'),
#'                          coverage_data=example_coverage)
#' ## changing either_side and style
#' either_side_GRange <- GRanges('chr5', IRanges(1000, 3200))
#' multipleInsertionDiagram(insertion_data=insertion_data,
#'                          genome_ideogram_data=ideogram_data,
#'                          style=c('ins1'=1, 'ins5'=4),
#'                          either_side=list('ins1'=500,
#'                                             'ins5'=either_side_GRange))
#'
multipleInsertionDiagram <- function(insertion_data,
    genome_ideogram_data, either_side = "default",
    track_height = 0.15, style = 1, colour_set = nice_colours,
    coverage_rectangle = NULL, coverage_data = NULL,
    label_data = NULL, label_colour = "black",
    label_size = 1, xaxis_spacing = "start_end") {
    ## check our key inputs:
    insertion_data <- .checkInsertionData(insertion_data,
        multiple = TRUE)
    genome_ideogram_data <- .checkIdeogramData(genome_ideogram_data)

    ## check other inputs the style of the diagrams
    ## can be set per event or overall
    if (length(style) == 1) {
        style_vector <- rep(style, nrow(insertion_data))
    } else if (length(style) == nrow(insertion_data)) {
        style_vector <- style
    } else {
        stop("style should be either a single number or a named vector with one
            element per event")
    }
    ## use the names of the vector to allocate
    ## specific styles to specific events
    if (is.null(names(style_vector))) {
        names(style_vector) <- insertion_data$name
    }

    ## same thing for either_side
    if (length(either_side) == 1) {
        either_side_list <- list()
        for (i in seq_along(insertion_data$name)) {
            either_side_list[insertion_data$name[i]] <- either_side
        }
    } else if (length(either_side) == nrow(insertion_data)) {
        either_side_list <- either_side
    } else {
        stop("either_side should be either a single number or a named list with
            one element per event")
    }

    ## firstly, we need to figure out where each
    ## plot/event lies around the genome circle in
    ## the middle so we know which slot around the
    ## page to plot it in (there are 9 slots total,
    ## laid out in a 3 x 3 grid). We can use the
    ## circlize function to do this:

    circos.clear()  # we need to initialise the whole genome to use circlize
    circos.par(start.degree = 90, unit.circle.segments = 1000,
        gap.after = 1)
    circos.initializeWithIdeogram(genome_ideogram_data,
        plotType = NULL)

    layout_points <- list()  # point in polar coordinates (theta + radius)
    layout_angles <- vector()  # angles and event names for easy access
    for (i in seq_along(insertion_data$name)) {
        this_event <- insertion_data[i, ]
        point <- circlize(x = mean(this_event$start,
            this_event$end), y = 1, sector.index = this_event$chr)

        ## we need to manually adjust by 90 degrees
        ## because our start degree for the genome
        ## circle is 90 degrees. Also bring it into the
        ## 0-360 range if need be with .simplifyAngle
        point[1] <- .simplifyAngle(point[1] - 90)
        angle <- point[1]

        ## add to list for use later
        suppressWarnings(layout_angles[this_event$name] <- angle)
        suppressWarnings(layout_points[this_event$name] <- list(point))
    }

    ## now we know where each event is, we can begin
    ## to allocate them to slots first, allocate
    ## them to areas (quarters) of the circle. (sort
    ## them because slot_assigner relies on the
    ## events being in ascending order)
    area_one <- sort(subset(layout_angles, layout_angles >=
        0 & layout_angles < 90))
    area_two <- sort(subset(layout_angles, layout_angles >=
        90 & layout_angles < 180))
    area_three <- sort(subset(layout_angles, layout_angles >=
        180 & layout_angles < 270))
    area_four <- sort(subset(layout_angles, layout_angles >=
        270 & layout_angles <= 360))

    ## assign slots
    layout_slots <- .slotAssigner(area_one, area_two,
        area_three, area_four)

    ## next, set up the rotation of each plot so
    ## that the chromosomes face the unzoomed
    ## counterpart on the genome plot (just depends
    ## on the slot)
    layout_rotation <- vector()
    rotation_angles <- c(`2` = 330, `1` = 0, `4` = 45,
        `7` = 90, `8` = 150, `9` = 180, `6` = 230,
        `3` = 270)
    for (i in seq_along(layout_slots)) {
        for (j in seq_along(names(rotation_angles))) {
            if (layout_slots[i] == names(rotation_angles)[j]) {
                layout_rotation[i] <- rotation_angles[j]
            }
        }
    }

    ## now we can move on to setting up the
    ## parameters for plotting:
    op = par(no.readonly = TRUE)  # store the original par() for later
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(name = "circles"))
    par(omi = gridBase::gridOMI())  # set up gridBase
    par(mar = c(0.8, 0, 0.8, 0), mfrow = c(3, 3))

    ## choose our colours & name them so we know
    ## each chromosome's colour
    all_colours <- colour_set[seq(1, nrow(genome_ideogram_data))]
    names(all_colours) <- as.character(genome_ideogram_data$chr)

    ## now plot, going slot by slot
    layout_link_points <- list()
    for (i in seq(1, 9)) {
        if (i %in% layout_slots) {
            # if we've allocated a plot here, plot it we
            # need to get a lot of information from our
            # lists so do that here and give them useful
            # names
            layout_slot_index <- grep(i, layout_slots)
            event_name <- names(layout_slots)[layout_slot_index]
            event_data <- insertion_data[insertion_data$name ==
                event_name, ]
            event_chr <- as.character(event_data$chr[1])  # all same chr
            event_either_side <- either_side_list[[event_name]]

            ## because of the way the labels are
            ## implemented, we need to subset out only this
            ## plot's labels. If we're doing this the we
            ## need to do the same subsetting for the colour
            ## since colours are just assigned by the order
            if (!is.null(label_data)) {
                event_labels <- label_data[seqnames(label_data) ==
                  event_chr]
                if (length(label_colour) > 1) {
                  event_label_colour <- subset(label_colour,
                    as.vector(grepl(event_chr,
                      seqnames(label_data))))
                } else {
                  event_label_colour <- label_colour
                }

            } else {
                event_labels <- NULL
                event_label_colour <- "black"
            }

            ## choose colours for this diagram (match chr
            ## colour to the centre circle; insertion colour
            ## is a tint of the insert's colour):
            event_colours <- c(all_colours[event_chr],
                colorspace::lighten(event_data$colour,
                  0.5))

            insertionDiagram(event_data,
                style = as.numeric(style_vector[event_name]),
                start_degree = layout_rotation[layout_slot_index],
                track_height = track_height, xaxis_spacing = xaxis_spacing,
                coverage_rectangle = coverage_rectangle,
                coverage_data = coverage_data,
                label_colour = event_label_colour,
                label_size = label_size, label_data = event_labels,
                internal = TRUE, either_side = event_either_side,
                sector_colours = event_colours,
                sector_border_colours = event_colours)

            ## calculate some points around the circle so
            ## the link follows along it (we'll use these
            ## later) first find where exactly the chr
            ## sector (cs) that we are drawing our link to
            ## is located (as in which track it's on):
            if (style_vector[event_name] %in% c(2,
                4)) {
                cs_track <- 1  # cs_track is the track where we draw a link to
            } else if (style_vector[event_name] ==
                3) {
                cs_track <- 3
            } else {
                cs_track <- 2
            }

            ## labels add two tracks: one for connector, one
            ## for the text
            if (!is.null(label_data)) {
                cs_track <- cs_track + 2
            }
            cs_ylim <- get.cell.meta.data("ylim",
                sector.index = event_chr, track.index = cs_track)
            cs_xlim <- get.cell.meta.data("xlim",
                sector.index = event_chr, track.index = cs_track)

            ## the ideogram box only goes 1/2 way up the
            ## track
            chr_box_top <- cs_ylim[1] + 0.5 * (cs_ylim[2] -
                cs_ylim[1])

            ## we need 3 types of points: the top corners of
            ## the ideogram box where the links will join in
            ## with this circle and the points along the
            ## bottom of the box where the line stops.  for
            ## the bottom of the box
            circle_points <- circlize(x = seq(from = cs_xlim[1],
                to = cs_xlim[2], length.out = 100),
                y = rep(cs_ylim[1], 100), sector.index = event_chr,
                track.index = cs_track)

            ## join the top left corner before the bottom of
            ## box
            circle_points <- rbind(circlize(x = cs_xlim[1],
                y = chr_box_top, sector.index = event_chr,
                track.index = cs_track), circle_points)

            ## join the top right corner after the bottom of
            ## box
            circle_points <- rbind(circle_points,
                circlize(x = cs_xlim[2], y = chr_box_top,
                  sector.index = event_chr, track.index = cs_track))

            ## convert the polar coordinates to cartesian
            circle_y_base <- sin(pracma::deg2rad(circle_points[,
                1])) * (circle_points[, 2] + 0.008)
            circle_x_base <- cos(pracma::deg2rad(circle_points[,
                1])) * (circle_points[, 2] + 0.008)
            circos.clear()

            ## convert the native coordinates to device
            ## coordinates (which we can later use in the
            ## grid graphics system)
            y_grid <- graphics::grconvertY(circle_y_base,
                "user", "ndc")
            x_grid <- graphics::grconvertX(circle_x_base,
                "user", "ndc")

            ## then add them all to a list
            suppressWarnings(layout_link_points[event_name] <- list(list(
                x_grid, y_grid)))

        } else {
            ## the 5th plot will always be the central
            ## genome circle:
            if (i == 5) {
                gmovizInitialise(genome_ideogram_data,
                  xaxis_spacing = "start_end",
                  sector_colours = all_colours,
                  sector_border_colours = all_colours)
                ## while we're plotting it, we need to find the
                ## coordinates where the links will connect with
                ## this circle.
                all_mc_points <- list()
                for (j in seq_along(insertion_data$name)) {
                  event_name <- insertion_data$name[j]
                  ## use the points we calculated earlier for slot
                  ## allocation
                  mc_point <- unlist(layout_points[event_name])

                  ## convert polar to cartesian, then user to
                  ## device coords (not sure why but when we do
                  ## this we need to reverse the -90 adjustment we
                  ## made when allocating the slots)
                  mc_y_base <- sin(pracma::deg2rad((mc_point[1] +
                    90))) * mc_point[2]
                  mc_x_base <- cos(pracma::deg2rad((mc_point[1] +
                    90))) * mc_point[2]
                  mc_y <- grconvertY(mc_y_base,
                    "user", "ndc")
                  mc_x <- grconvertX(mc_x_base,
                    "user", "ndc")
                  suppressWarnings(all_mc_points[event_name] <- list(list(mc_x,
                    mc_y)))
                }
                circos.clear()
            } else {
                ## otherwise, plot something blank so we can
                ## skip to the next slot (doesn't matter what it
                ## is)
                circos.initialize(factors = letters[seq(1,
                  5)], xlim = c(0, 1))
                circos.clear()
            }
        }
    }

    ## to plot the links we need to switch to grid
    ## graphics
    grid::pushViewport(grid::viewport(x = 0.5,
        y = 0.5, width = 1, height = 1, just = "centre",
        clip = "off"))
    for (i in seq_along(insertion_data$name)) {
        ## again just name some important info so we can
        ## use it easily
        event_name <- insertion_data$name[i]
        event_data <- insertion_data[insertion_data$name ==
            event_name, ]
        event_chr <- as.character(event_data$chr[1])  # all same chr

        ## generate the points for the link
        link_points <- .makeLinkPoints(slot = layout_slots[event_name],
            x_grid = layout_link_points[[i]][[1]],
            y_grid = layout_link_points[[i]][[2]],
            mc_x = all_mc_points[[i]][[1]], mc_y = all_mc_points[[i]][[2]])

        ## finally, plot the link polygon
        grid::grid.polygon(x = link_points[[1]],
            y = link_points[[2]], gp = grid::gpar(col = all_colours[event_chr],
                fill = paste0(all_colours[event_chr],
                  "33")))
    }

    ## return par to original values
    par(op)
}
