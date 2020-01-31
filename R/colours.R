utils::globalVariables(c(
    "nice_colours", "pastel_colours", "rich_colours",
    "bright_colours_transparent", "bright_colours_opaque"))
#' gmoviz colour sets
#'
#' @name colourSets
#' @details Due to the often high number of sectors being plotted with gmoviz
#' (e.g. 20+ when plotting each chromosome), a number of 'colour sets' have
#' been included for convenience.
#' @format Character vectors of 34 hex colours.
#' @source Many of the colours are from, or inspired by ColorBrewer
#' \url{http://colorbrewer2.org/}.
NULL

#' \code{nice_colours}: The default colour set. Medium brightness, light
#' colours. Designed for use on a white/light coloured background.
#' @rdname colourSets
"nice_colours"

#' \code{pastel_colours}:A set of pale/pastel colours, modelled on the
#' \code{nice_colours} set but less saturated. Designed for use on a
#' white/light background.
#' @rdname colourSets
"pastel_colours"

#' \code{rich_colours}: A set of bright, vibrant colours (but not neon, like
#' the bright_colours_transparent). Designed for use on any sort of background.
#' @rdname colourSets
"rich_colours"

#' \code{bright_colours_transparent}: A set of very bright (neon) colours with
#' slight transparency. Designed for use on a black background.
#' @rdname colourSets
"bright_colours_transparent"

#' \code{bright_colours_opaque}: A set of very bright (neon) colours without
#' transparency. Designed for use on a black background.
#' @rdname colourSets
"bright_colours_opaque"
