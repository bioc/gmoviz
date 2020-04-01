#### checkIdeogramData ####
missing_columns <- data.frame(chr = "chr1", start = 1)
nonsense_start_end <- data.frame(chr = "chr1", start = 1000, end = 0)
duplicated_sector_names <- GRanges(c("chr1", "chr2", "chr1"),
                                   IRanges(start = c(0, 0, 1000),
                                           end = c(1000, 1000, 3000)))
passing_GRanges <- GRanges(c("chr1", "chr2", "chr3"),
                           IRanges(start = c(0, 0, 0),
                                   end = c(1000, 2000, 3000)))
passing_dataframe <- data.frame(chr = c("chr1", "chr2", "chr3"),
                                start = c(0, 0, 0),
                                end = c(1000, 2000, 3000), 
                                stringsAsFactors = TRUE)
test_that("checking ideo data with checkIdeogramData is working", {
    expect_error(gmoviz:::.checkIdeogramData("test"))
    expect_error(gmoviz:::.checkIdeogramData(missing_columns))
    expect_error(gmoviz:::.checkIdeogramData(nonsense_start_end))
    expect_error(gmoviz:::.checkIdeogramData(duplicated_sector_names))
    expect_equivalent(gmoviz:::.checkIdeogramData(passing_GRanges),
                      passing_dataframe)
})

#### checkLabelData ####
nonsense_start_end <- data.frame(chr = "chr1", start = 1000, end = 0,
                                 label = "blah")
passing_GRanges <- GRanges(c("chr1", "chr2", "chr3"),
                           IRanges(start = c(900, 1900, 2900),
                                   end = c(1000, 2000, 3000)),
                           label = c("a", "b", "c"))
passing_dataframe <- data.frame(chr = as.factor(c("chr1", "chr2", "chr3")),
                                start = c(900, 1900, 2900),
                                end = c(1000, 2000, 3000),
                                label = c("a", "b", "c"),
                                stringsAsFactors = FALSE)
test_that("checking label data with checkLabelData is working", {
    expect_error(gmoviz:::.checkLabelData("test"))
    expect_error(gmoviz:::.checkLabelData(missing_columns))
    expect_error(gmoviz:::.checkLabelData(nonsense_start_end))
    expect_equivalent(gmoviz:::.checkLabelData(passing_GRanges),
                      passing_dataframe)
})

#### checkFeatureData ####
nonsense_start_end <- data.frame(chr = "chr1", start = 1000, end = 0)
wrong_shapes <- GRanges("chr1", IRanges(1000, 2000), shape = "blah")
wrong_tracks <- GRanges("chr1", IRanges(1000, 2000), track = -4)
default_GRanges <- GRanges("chr1", IRanges(1000, 2000))
default_dataframe <- data.frame(chr = as.factor("chr1"), start = 1000,
                                end = 2000, label = "", shape = "rectangle",
                                track = 1, colour = "#f1d03c",
                                stringsAsFactors = FALSE)
passing_GRanges <- GRanges(c("chr1", "chr2", "chr3"),
                           IRanges(start = c(900, 1900, 2900),
                                   end = c(1000, 2000, 3000)),
                           label = c("a", "b", "c"),
                           colour = c("red", "blue", "yellow"),
                           shape = c("rectangle", "forward_arrow",
                                     "upwards_triangle"),
                           track = c(0, 1, 2))
passing_dataframe <- data.frame(chr = as.factor(c("chr1", "chr2", "chr3")),
                                start = c(900, 1900, 2900),
                                end = c(1000, 2000, 3000),
                                label = c("a", "b", "c"),
                                colour = c("red", "blue", "yellow"),
                                shape = c("rectangle", "forward_arrow",
                                          "upwards_triangle"),
                                track = c(0, 1, 2), stringsAsFactors = FALSE)
test_that("checking feature data with checkFeatureData is working", {
    expect_error(gmoviz:::.checkFeatureData("test"))
    expect_error(gmoviz:::.checkFeatureData(missing_columns))
    expect_error(gmoviz:::.checkFeatureData(nonsense_start_end))
    expect_error(gmoviz:::.checkFeatureData(wrong_shapes))
    expect_error(gmoviz:::.checkFeatureData(wrong_tracks))
    expect_equivalent(gmoviz:::.checkFeatureData(passing_GRanges),
                      passing_dataframe)
    expect_equivalent(gmoviz:::.checkFeatureData(default_GRanges),
                      default_dataframe)
})

#### checkInsertionData ####
nonsense_start_end <- data.frame(chr = "chr1", start = 1000, end = 0,
                                 name = "insA", length = 300)
wrong_shapes <- GRanges("chr1", IRanges(1000, 2000), name = "insA",
                        length = 300, shape = "blah")
wrong_length <- GRanges("chr1", IRanges(1000, 2000), name = "insA", length = 0)
wrong_in_tandem <- GRanges("chr1", IRanges(1000, 2000), name = "insA",
                           length = 100, in_tandem = 0)
duplicated_insertion_names <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
                                      IRanges(start = c(1000, 2000, 3000),
                                              end = c(1100, 2100, 3100)),
                                      name = c("insA", "insB", "insA"),
                                      length = c(100, 200, 100))
default_GRanges <- GRanges("chr1", IRanges(1000, 2000), name = "insA",
                           length = 150)
default_dataframe <- data.frame(chr = as.factor("chr1"), start = 1000,
                                end = 2000, name = "insA", length = 150,
                                shape = "forward_arrow", in_tandem = 1,
                                colour = "#f1d03c", stringsAsFactors = FALSE)
passing_GRanges <- GRanges(
  c("chr1", "chr1", "chr1"),
  IRanges(start = c(900, 1900, 2900), end = c(1000, 2000, 3000)),
  name = c("a", "b", "c"), length = c(100, 150, 200),
  colour = c("red", "blue", "yellow"),
  shape = c("rectangle", "forward_arrow", "reverse_arrow"),
  in_tandem = c(1, 1, 2))
passing_dataframe <- data.frame(
  chr = as.factor(c("chr1", "chr1", "chr1")),
  start = c(900, 1900, 2900), end = c(1000, 2000, 3000),
  name = c("a", "b", "c"), length = c(100, 150, 200),
  colour = c("red", "blue", "yellow"),
  shape = c("rectangle", "forward_arrow", "reverse_arrow"),
  in_tandem = c(1, 1, 2), stringsAsFactors = FALSE)
test_that("checking insertion data with checkInsertionData is working", {
    expect_error(gmoviz:::.checkInsertionData("test"))
    expect_error(gmoviz:::.checkInsertionData(missing_columns))
    expect_error(gmoviz:::.checkInsertionData(nonsense_start_end))
    expect_error(gmoviz:::.checkInsertionData(wrong_shapes))
    expect_error(gmoviz:::.checkInsertionData(wrong_length))
    expect_error(gmoviz:::.checkInsertionData(wrong_in_tandem))
    expect_error(gmoviz:::.checkInsertionData(duplicated_insertion_names))
    expect_equivalent(gmoviz:::.checkInsertionData(passing_GRanges),
                      passing_dataframe)
    expect_equivalent(gmoviz:::.checkInsertionData(default_GRanges),
                      default_dataframe)
})

#### checkCoverageData ####
nonsense_start_end <- data.frame(chr = "chr1", start = 1000, end = 0,
                                 coverage = 300)
wrong_coverage <- data.frame(chr = "chr1", start = 1000, end = 0,
                             coverage = -5)
passing_GRanges <- GRanges(c("chr1", "chr2", "chr3"),
                           IRanges(start = c(0, 0, 0),
                                   end = c(1000, 2000, 3000)),
                           coverage = c(0, 4, 100))
passing_dataframe <- data.frame(chr = c("chr1", "chr2", "chr3"),
                                start = c(0, 0, 0),
                                end = c(1000, 2000, 3000),
                                coverage = c(0, 4, 100),
                                stringsAsFactors = TRUE)
test_that("checking coverage data with checkCoverageData is working", {
    expect_error(gmoviz:::.checkCoverageData("test"))
    expect_error(gmoviz:::.checkCoverageData(missing_columns))
    expect_error(gmoviz:::.checkCoverageData(nonsense_start_end))
    expect_error(gmoviz:::.checkCoverageData(wrong_coverage))
    expect_equivalent(gmoviz:::.checkCoverageData(passing_GRanges),
                      passing_dataframe)
})

#### checkSectorsMatch ####
passing_dataframe <- data.frame(chr = c("a", "b", "c"),
                                start = c(0, 0, 0),
                                end = c(0.5, 0.6, 1),
                                coverage = c(0, 4, 100))
failing_dataframe <- data.frame(chr = c("A", "b", "c"),
                                start = c(0, 0, 0),
                                end = c(0.5, 0.6, 1),
                                coverage = c(0, 4, 100))
circos.initialize(factors = letters[1:10], xlim = c(0, 1))
test_that("checking sectors to plot into with checkSectorsMatch is working", {
    expect_equal(gmoviz:::.checkSectorsMatch(passing_dataframe), TRUE)
    expect_warning(gmoviz:::.checkSectorsMatch(failing_dataframe))
})
