context("Data import functions")
library(gmoviz)
library(GenomicRanges)
library(BiocGenerics)
library(GenomeInfoDb)
library(rtracklayer)

#### features from gff ####
coloured_by_type <- GRanges(
  rep("chr4", 5),
  IRanges(start = c(680000, 360000, 1080000, 120000, 840000),
          end = c(760000, 400000, 1200000, 160000, 880000)),
  label = c("geneX", "InsertA", "geneY", "InsertB", "InsertC"),
  colour = c("#8dd3c7", "#e7ec83", "#8dd3c7", "#e7ec83", "#e7ec83"),
  shape = c("forward_arrow", "rectangle", "forward_arrow",
            "rectangle", "rectangle"),
  track = rep(1, 5))

coloured_by_name <- GRanges(
  rep("chr4", 5),
  IRanges(start = c(680000, 360000, 1080000, 120000, 840000),
          end = c(760000, 400000, 1200000, 160000, 880000)),
  label = c("geneX", "InsertA", "geneY", "InsertB", "InsertC"),
  colour = c("#8dd3c7", "#e7ec83", "#bebada", "#fb8072", "#fccde5"),
  shape = c("forward_arrow", "rectangle", "forward_arrow",
            "rectangle", "rectangle"),
  track = rep(1, 5))

path <- system.file("extdata", "example.gff3", package = "gmoviz")
test_that("feature import with getFeatures is working", {
  expect_equivalent(getFeatures(path), coloured_by_type)
  expect_equivalent(getFeatures(path, colour_by_type = FALSE),
                    coloured_by_name)
})

#### labels from gff ####
example <- GRanges(
  rep("chr4", 5),
  IRanges(start = c(680000, 360000, 1080000, 120000, 840000),
          end = c(760000, 400000, 1200000, 160000, 880000)),
  label = c("geneX", "InsertA", "geneY", "InsertB", "InsertC"),
  type = c("Gene", "Insert", "Gene", "Insert", "Insert"),
  colour = c("#ff0000", "#00c8ff", "#ff0000", "#00c8ff", "#00c8ff"))

example_not_colour_by_type <- example
example_not_colour_by_type$colour <- NULL

path <- system.file("extdata", "example.gff3", package = "gmoviz")
test_that("gene label import with getLabels is working", {
  expect_equivalent(getLabels(path), example)
  expect_equivalent(getLabels(path, colour_code = FALSE),
                    example_not_colour_by_type)
})

#### coverage ####
path <- system.file("extdata", "ex1.bam", package="Rsamtools")
range <- GRanges("seq1", IRanges(start = 0, end = 20))

normal <- GRanges(rep("seq1", 21), IRanges(start = seq(0, 20),
                                           end = seq(1, 21)),
                  coverage = c(1, 1, 2, 2, 3, 4, 4, 4, 5, 5, 5, 5, 7, 7, 8, 8,
                               8, 9, 9, 9, 9))
smooth <- GRanges(rep("seq1", 21), IRanges(start = seq(0, 20),
                                           end = seq(1, 21)),
                  coverage = c(1.000000, 1.000000, 1.333333, 1.500000,
                               1.800000, 2.400000, 3.000000, 3.400000,
                               4.000000, 4.400000, 4.600000, 4.800000,
                               5.400000, 5.800000, 6.400000, 7.000000,
                               7.600000, 8.000000, 8.400000, 8.600000,
                               8.800000))
error_range <- GRanges("chr1", IRanges(start = 0, end = 300))

test_that("coverage data import from bam with getCoverage is working", {
  normal_coverage <- getCoverage(range, path)
  names(normal_coverage) <- NULL
  expect_equivalent(normal_coverage, normal)
  smoothed_coverage <- getCoverage(range, path, smoothing_window_size = 5)
  names(smoothed_coverage) <- NULL
  expect_equivalent(smoothed_coverage, smooth)
  expect_error(getCoverage(error_range, path))
})

#### ideogram ####
# bam file
whole <- GRanges(seqnames = c("seq1", "seq2"),
                 ranges = IRanges(start = c(0, 0), end = c(1575, 1584)))

seq1_only <- GRanges("seq1", IRanges(0, 1575))
seq2_only <- GRanges("seq2", IRanges(0, 1584))

path <- system.file("extdata", "ex1.bam", package="Rsamtools")
test_that("ideo data import of .bam with getIdeogramData is working", {
  expect_equivalent(getIdeogramData(path, just_pattern = "^seq"), whole)
  expect_equivalent(getIdeogramData(path, wanted_chr = "seq1"), seq1_only)
  expect_equivalent(getIdeogramData(path, unwanted_chr = "seq1",
                                        just_pattern = "^seq"), seq2_only)
})
# fasta & fasta folder
both <- GRanges(c("HSBGPG Human gene for bone gla protein (BGP)",
                  "HSGLTH1 Human theta 1-globin gene"),
                IRanges(start = c(0, 0), end = c(1231, 1020)))
bgp_only <- keepSeqlevels(both, "HSBGPG Human gene for bone gla protein (BGP)",
                          pruning.mode = "coarse")
tglob_only <- keepSeqlevels(both, "HSGLTH1 Human theta 1-globin gene",
                            pruning.mode = "coarse")
# fasta
path <- system.file("extdata", "fastaFolder/HSBGPG.fasta", package = "gmoviz")
test_that("ideo data import of .fasta with getIdeogramData is working", {
  expect_equivalent(getIdeogramData(fasta_file = path), bgp_only)
})

# fasta_folder
path <- system.file("extdata", "fastaFolder", package = "gmoviz")
test_that("ideo data import of .fasta with getIdeogramData is working", {
  expect_equivalent(getIdeogramData(fasta_folder = path), both)
  expect_equivalent(getIdeogramData(fasta_folder = path,
                                    just_pattern = "HSBGPG"), bgp_only)
  expect_equivalent(getIdeogramData(fasta_folder = path,
                                    just_pattern = "HSGLTH1"), tglob_only)
})
