context("Various helper functions")
library(gmoviz)

without_1 <- c("1", "2", "3", "4", "5", "6", "X")
without_2 <- c("1", "2", "3", "4", "5", "14", "19")
without_3 <- c("1", "2", "3", "4", "5", "6", "plasmid")
with_1 <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chrX")
with_2 <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr14", "chr19")
with_3 <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "plasmid")

#### addChrToNumber ####
test_that("Adding chr to numbers w/ addChrToNumber works", {
  expect_equivalent(gmoviz:::.addChrToNumber(without_1), as.factor(with_1))
  expect_equivalent(gmoviz:::.addChrToNumber(without_2), as.factor(with_2))
  expect_equivalent(gmoviz:::.addChrToNumber(without_3), as.factor(with_3))
})

#### GRangesToBed ####
library(GenomicRanges)
GRange_1 <- GRanges(c("chr1", "chr2", "chr3"),
               IRanges(start = c(0, 0, 0), end = c(10000, 20000, 14000)))
GRange_2 <- GRanges(c("chr1", "chr2", "chr3"),
               IRanges(start = c(0, 0, 0), end = c(10000, 20000, 14000)),
               val1 = c(133, 144, 155), val2 = c("h", "i", "j"))
dataframe_1 <- data.frame(chr = c("chr1", "chr2", "chr3"), start = c(0, 0, 0),
                  end = c(10000, 20000, 14000), stringsAsFactors = TRUE)
dataframe_2 <- data.frame(chr = c("chr1", "chr2", "chr3"),
                  start = c(0, 0, 0), end = c(10000, 20000, 14000),
                  val1 = c(133, 144, 155), val2 = c("h", "i", "j"),
                  stringsAsFactors = TRUE)

test_that("Converting GRanges to bed format data frame is working", {
  expect_equivalent(gmoviz:::.GRangesToBed(GRange_1, ideogram = TRUE),
                    dataframe_1)
  expect_equivalent(gmoviz:::.GRangesToBed(GRange_2, ideogram = TRUE),
                    dataframe_1)
  expect_equivalent(gmoviz:::.GRangesToBed(GRange_2, ideogram = FALSE),
                    dataframe_2)
})

#### regex maker ####
test_that("constructing regex with makeRegex is working", {
  expect_equal(gmoviz:::.makeRegex("chr1"), "^chr1$")

  expect_equal(gmoviz:::.makeRegex(c("chr1", "chr13", "chr2")),
               "^chr1$|^chr13$|^chr2$")

  expect_equal(gmoviz:::.makeRegex(c("chrX", "chr15", "chr4")),
               "^chrX$|^chr15$|^chr4$")
})

#### getThisSectorCoverage ####
coverage <- data.frame(chr = c("chr1", "chr2", "chr2", "chr3", "chr5", "chr7"),
                       start = rep(0, 6), end = rep(4, 6),
                       coverage = c(1,2,3,4,5,6))
test_that("filtering coverage df with getThisSectorCoverage is working", {
  expect_equivalent(
    gmoviz:::.getThisSectorCoverage(coverage, "chr1"), coverage[1, ])
  expect_equivalent(
    gmoviz:::.getThisSectorCoverage(coverage, c("chr1", "chr2")),
    coverage[1:3, ])
  expect_equivalent(
    gmoviz:::.getThisSectorCoverage(coverage, "chr7"), coverage[6, ])

})

#### insertionsToFeatures ####
insertion_data_1 <- data.frame(
  chr = "chr12", start = 70905597, end = 70917885, name = "plasmid",
  colour = "#7270ea", length = 12000, in_tandem = 11, shape = "forward_arrow",
  stringsAsFactors = FALSE)
insertion_data_2 <- data.frame(
  chr = rep("chr12", 2), start = c(70905597, 70705597),
  end = c(70917885, 70717885), name = c("plasmid1", "plasmid2"),
  colour = c("#7270ea", "#ea7082"), length = c(12000, 10000),
  in_tandem = c(4, 1), shape = c("rectangle", "forward_arrow"),
  stringsAsFactors = FALSE)
feature_data_1 <- data.frame(
  chr = rep("plasmid", 11), start = seq(from = 0, by = 12000, length.out = 11),
  end = seq(from = 12000, by = 12000, length.out = 11),
  colour = rep("#7270ea", 11), label = 1:11, shape = rep("forward_arrow", 11),
  track = rep(1, 11), stringsAsFactors = FALSE)
feature_data_2 <- data.frame(
  chr = c("plasmid1", "plasmid1", "plasmid1", "plasmid1", "plasmid2"),
  start = c(seq(from = 0, by = 12000, length.out = 4),
            seq(from = 0, by = 10000, length.out = 1)),
  end = c(seq(from = 12000, by = 12000, length.out = 4),
          seq(from = 10000, by = 10000, length.out = 1)),
  colour = c("#7270ea", "#7270ea", "#7270ea", "#7270ea", "#ea7082"),
  label = c("1", "2", "3", "4", "plasmid2"),
  shape = c("rectangle", "rectangle", "rectangle", "rectangle",
            "forward_arrow"), track = rep(1, 5),
  stringsAsFactors = FALSE)

test_that("conversion of insertion DF to feature DF with insertionsToFeatures
          is working", {
  expect_equivalent(gmoviz:::.insertionsToFeatures(insertion_data_1),
                    feature_data_1)
  expect_equivalent(gmoviz:::.insertionsToFeatures(insertion_data_2),
                    feature_data_2)

})

#### insertionsToIdeogram ####
ideogram_1 <- data.frame(
  chr = as.factor(c("chr12", "plasmid")), start = c(70905597 - 1843.2, -4800),
  end = c(70917885 + 1843.2, 136800))
ideogram_2 <- data.frame(
  chr = as.factor(c("chr12", "plasmid1", "plasmid2")),
  start = c(70705597 - 2000, -4800, -4000),
  end = c(70917885 + 2000, 52800, 14000))
ideogram_3 <- data.frame(
  chr = as.factor(c("chr12", "plasmid")), start = c(60905597, -4800),
  end = c(80917885, 136800))
GRange_3 <- GRanges("chr12", IRanges(60905597, 80917885))

test_that("conversion of insertion DF to ideo DF with insertionsToIdeogram is
          working", {
  expect_equivalent(
    gmoviz:::.insertionsToIdeogram(insertion_data_1, either_side = "default"),
    ideogram_1)
  expect_equivalent(
    gmoviz:::.insertionsToIdeogram(insertion_data_2, 2000), ideogram_2)
  expect_equivalent(
    gmoviz:::.insertionsToIdeogram(insertion_data_1, c(60905597, 80917885)),
    ideogram_3)
  expect_equivalent(
    gmoviz:::.insertionsToIdeogram(insertion_data_1, GRange_3), ideogram_3)
})

#### simplifyAngle ####
test_that("bringing angles back into the 0-360 degree range with
          simplifyAngle is working", {
            expect_equal(gmoviz:::.simplifyAngle(90), 90)
            expect_equal(gmoviz:::.simplifyAngle(360), 360)
            expect_equal(gmoviz:::.simplifyAngle(490), 130)
            expect_equal(gmoviz:::.simplifyAngle(-45), 315)
})
