context("MultipleInsertionDiagram helper functions")

#### .equivSlotFinder ####
test_that(
  "finding equivalent slots in other areas with .equivSlotFinder is working",
  {
    expect_equal(gmoviz:::.equivSlotFinder(2, "previous"), 6)
    expect_equal(gmoviz:::.equivSlotFinder(2, "next"), 4)
    expect_equal(gmoviz:::.equivSlotFinder(9, "previous"), 7)
    expect_equal(gmoviz:::.equivSlotFinder(9, "next"), 3)
  })

#### .firstOrLastTwo ####
as_1 <- c("ins1" = 2, "ins2" = 1, "ins3" = 9) # they need to be in order
as_2 <- c("ins1" = 1, "ins2" = 4, "ins3" = 9, "ins4" = 6)
test_that(
  "determining whether plots occupy the first or last two slots of an area with
  .firstOrLastTwo is working",
  {
    expect_equal(
      gmoviz:::.firstOrLastTwo(c("ins1", "ins2"), as_1), "first_two")
    expect_equal(
      gmoviz:::.firstOrLastTwo(c("ins3", "ins4"), as_2), "last_two")
    expect_equal(
      gmoviz:::.firstOrLastTwo(c("ins1", "ins2"), as_2), "last_two")
  })

#### .chooseTwoSlots ####
ea_1 <- c("ins4" = 220, "ins5" = 227)
ea_2 <- c("ins4" = 10, "ins5" = 40)
ea_3 <- c("ins4" = 222, "ins5" = 230)
ea_4 <- c("ins4" = 49, "ins5" = 60)
as_1_1 <- c("ins1" = 2, "ins2" = 1, "ins3" = 9, "ins5" = 7, "ins4" = 4)
as_1_2 <- c("ins1" = 2, "ins2" = 1, "ins3" = 9, "ins5" = 3, "ins4" = 6)
as_1_3 <- c("ins1" = 2, "ins2" = 1, "ins3" = 9, "ins4" = 7, "ins5" = 8)
as_1_4 <- c("ins1" = 2, "ins2" = 1, "ins3" = 9, "ins4" = 3, "ins5" = 2)
test_that(
  "choosing two slots in an area to put plots in with .chooseTwoSlots is
  working",
  {
    expect_equal(gmoviz:::.chooseTwoSlots(ea_1, c(4, 7, 8), 225, as_1), as_1_1)
    expect_equal(gmoviz:::.chooseTwoSlots(ea_2, c(6, 3, 2), 45, as_1), as_1_2)
    expect_equal(gmoviz:::.chooseTwoSlots(ea_3, c(4, 7, 8), 225, as_1), as_1_3)
    expect_equal(gmoviz:::.chooseTwoSlots(ea_4, c(6, 3, 2), 45, as_1), as_1_4)
  })

#### .slotAssigner ####
a1_1 <- c("ins1" = 145)
a2_1 <- c()
a3_1 <- c("ins2" = 300, "ins3" = 320, "ins4" = 350)
a4_1 <- c("ins5" = 35)
a1_2 <- c("ins1" = 145, "ins6" = 130)
a2_2 <- c()
a3_2 <- c("ins2" = 300, "ins3" = 320, "ins4" = 350)
a4_2 <- c("ins5" = 35, "ins7" = 50)
as_1 <- c("ins1" = 1, "ins2" = 8, "ins3" = 9, "ins4" = 6, "ins5" = 3)
as_2 <- c("ins1" = 1, "ins6" = 4, "ins2" = 8, "ins3" = 9, "ins4" = 6,
          "ins5" = 3, "ins7" = 2)
test_that("assigning plots to slots with .slotAssigner is working", {
  expect_equal(gmoviz:::.slotAssigner(a1_1, a2_1, a3_1, a4_1), as_1)
  expect_equal(gmoviz:::.slotAssigner(a1_2, a2_2, a3_2, a4_2), as_2)
})
