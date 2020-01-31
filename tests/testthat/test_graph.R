context("Graphing functions")
library(gmoviz)
library(grid)
library(ComplexHeatmap)
#### legends ####
genelabels <- data.frame(chr = c("chr1", "chr1"), start = c(100, 300),
                         end = c(150, 350), label = c("a", "b"),
                         type = c("gene", "lncRNA"), colour = c("red", "blue"),
                         stringsAsFactors = FALSE)
test_that("Make Legends is working properly", {
  expect_false(makeLegends())
  expect_error(makeLegends(label_legend = TRUE))
  expect_equivalent(
    class(makeLegends(label_legend = TRUE, label_data = genelabels)),
    "Legends")
    expect_equivalent(class(makeLegends(
    scatterplot_legend = TRUE, point_colour = c("red", "blue"))), "Legends")

  expect_equivalent(class(makeLegends(
      scatterplot_legend = TRUE, point_colour = c("red", "blue"),
      linegraph_legend = TRUE)), "Legends")
})
