context("Initialisation functions")
library(gmoviz)
#### sortIdeogramData ####
sorted_ideogram_data <- data.frame(
  chr = c("chr1", "chr3", "chrX"), start = c(0, 0, 0),
  end = c(1000, 3000, 2000))
unsorted_ideogram_data <- data.frame(
  chr = c("chr1", "chrX", "chr3"), start = c(0, 0, 0),
  end = c(1000, 2000, 3000))

test_that("Sorting of the regions to plot w/ sortIdeogramData", {
  expect_equal(
    gmoviz:::.sortIdeogramData(sorted_ideogram_data), sorted_ideogram_data)
  expect_equivalent(
    gmoviz:::.sortIdeogramData(unsorted_ideogram_data), sorted_ideogram_data)
})

#### makeZoomedIdeogramData ####
single_zoom_not_removed <- data.frame(
  chr = c("chr1", "chr3", "chrX", "zoomed_chr3"), start = c(0, 0, 0, 0),
  end = c(1000, 3000, 2000, 3000))
multi_zoom_not_removed <- data.frame(
  chr = c("chr1", "chr3", "chrX", "zoomed_chr3", "zoomed_chrX"),
  start = c(0, 0, 0, 0, 0), end = c(1000, 3000, 2000, 3000, 2000))
single_zoom_removed <- data.frame(
  chr = c("chr1", "chrX", "zoomed_chr3"), start = c(0, 0, 0),
  end = c(1000, 2000, 3000))
multi_zoom_removed <- data.frame(
  chr = c("chr1", "zoomed_chr3", "zoomed_chrX"), start = c(0, 0, 0),
  end = c(1000, 3000, 2000))

test_that("Zoomed data frame sets up correctly w/ makeZoomedIdeogramData", {
  expect_equivalent(
    gmoviz:::.makeZoomedIdeogramData(sorted_ideogram_data, chromosome = "chr3",
                                    remove_unzoomed = FALSE),
    single_zoom_not_removed)

  expect_equivalent(
    gmoviz:::.makeZoomedIdeogramData(sorted_ideogram_data,
                                    chromosome = c("chrX", "chr3"),
                                    remove_unzoomed = FALSE),
    multi_zoom_not_removed)

  expect_equivalent(
    gmoviz:::.makeZoomedIdeogramData(sorted_ideogram_data, chromosome = "chr3",
                                    remove_unzoomed = TRUE),
    single_zoom_removed)
  expect_equivalent(
    gmoviz:::.makeZoomedIdeogramData(sorted_ideogram_data,
                                    chromosome = c("chrX", "chr3"),
                                    remove_unzoomed = TRUE),
    multi_zoom_removed)
})
