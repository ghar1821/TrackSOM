library(data.table)

test_that("split_found_1", {
  dat <- data.table::data.table(
    SCA1=sample(seq(100,400), size=10),
    timepoint=c(rep(1,5), rep(2,5)),
    lineage=c(rep('A',5), rep('A', 3), rep('A|1', 2))
  )
  edge.df <- GetTransitionsAsEdges(dat, c(1,2), 'timepoint', 'lineage')
  expect_true(helper.found.edge(edge.df, '1_A', '2_A'))
  expect_true(helper.found.edge(edge.df, '1_A', '2_A|1'))
})

test_that("split_found_2", {
  dat <- data.table::data.table(
    SCA1=sample(seq(100,400), size=15),
    timepoint=c(rep(1,5), rep(2,5), rep(3,5)),
    lineage=c(rep('A',5), 
              rep('A', 3), rep('A|1', 2), 
              rep('A', 2), rep('A|1', 2), rep('A|2',1))
  )
  edge.df <- GetTransitionsAsEdges(dat, c(1,2,3), 'timepoint', 'lineage')
  expect_true(helper.found.edge(edge.df, '1_A', '2_A'))
  expect_true(helper.found.edge(edge.df, '1_A', '2_A|1'))
  expect_true(helper.found.edge(edge.df, '2_A', '3_A'))
  expect_true(helper.found.edge(edge.df, '2_A', '3_A|2'))
  expect_true(helper.found.edge(edge.df, '2_A|1', '3_A|1'))
})

test_that("merge_found_1", {
  dat <- data.table::data.table(
    SCA1=sample(seq(100,400), size=10),
    timepoint=c(rep(1,5), rep(2,5)),
    lineage=c(rep('A',3), rep('B', 2),
              rep('(A,B)', 5))
  )
  edge.df <- GetTransitionsAsEdges(dat, c(1,2), 'timepoint', 'lineage')
  expect_true(helper.found.edge(edge.df, '1_A', '2_(A,B)'))
  expect_true(helper.found.edge(edge.df, '1_B', '2_(A,B)'))
})

test_that("merge_found_2", {
  dat <- data.table::data.table(
    SCA1=sample(seq(100,400), size=15),
    timepoint=c(rep(1,5), rep(2,5), rep(3,5)),
    lineage=c(rep('A',3), rep('B', 2),
              rep('(A,B)', 3), rep('C', 2),
              rep('((A,B),C)', 2), rep('C|1',3))
  )
  edge.df <- GetTransitionsAsEdges(dat, c(1,2,3), 'timepoint', 'lineage')
  expect_true(helper.found.edge(edge.df, '1_A', '2_(A,B)'))
  expect_true(helper.found.edge(edge.df, '1_B', '2_(A,B)'))
  expect_true(helper.found.edge(edge.df, '2_(A,B)', '3_((A,B),C)'))
  expect_true(helper.found.edge(edge.df, '2_C', '3_C|1'))
})

test_that("merge_split_found_1", {
  dat <- data.table::data.table(
    SCA1=sample(seq(100,400), size=10),
    timepoint=c(rep(1,5), rep(2,5)),
    lineage=c(rep('A',3), rep('B', 2),
              rep('(A,B)', 3),
              rep('A|1', 2))
  )
  edge.df <- GetTransitionsAsEdges(dat, c(1,2), 'timepoint', 'lineage')
  expect_true(helper.found.edge(edge.df, '1_A', '2_(A,B)'))
  expect_true(helper.found.edge(edge.df, '1_B', '2_(A,B)'))
  expect_true(helper.found.edge(edge.df, '1_A', '2_A|1'))
})

test_that("merge_split_found_2", {
  dat <- data.table::data.table(
    SCA1=sample(seq(100,400), size=10),
    timepoint=c(rep(1,5), rep(2,5)),
    lineage=c(rep('A',3), rep('B', 2),
              rep('(A|1,B)', 2),
              rep('A', 3))
  )
  edge.df <- GetTransitionsAsEdges(dat, c(1,2), 'timepoint', 'lineage')
  expect_true(helper.found.edge(edge.df, '1_A', '2_(A|1,B)'))
  expect_true(helper.found.edge(edge.df, '1_B', '2_(A|1,B)'))
  expect_true(helper.found.edge(edge.df, '1_A', '2_A'))
})

test_that("merge_split_found_3", {
  dat <- data.table::data.table(
    SCA1=sample(seq(100,400), size=15),
    timepoint=c(rep(1,5), rep(2,5), rep(3,5)),
    lineage=c(rep('A',3), rep('B', 2),
              rep('(A,B)', 3),
              rep('A|1', 2),
              rep('(A,B)', 2),
              rep('(A,B)|1', 1),
              rep('A|1', 2))
  )
  edge.df <- GetTransitionsAsEdges(dat, c(1,2,3), 'timepoint', 'lineage')
  expect_true(helper.found.edge(edge.df, '1_A', '2_(A,B)'))
  expect_true(helper.found.edge(edge.df, '1_B', '2_(A,B)'))
  expect_true(helper.found.edge(edge.df, '1_A', '2_A|1'))
  expect_true(helper.found.edge(edge.df, '2_(A,B)', '3_(A,B)'))
  expect_true(helper.found.edge(edge.df, '2_(A,B)', '3_(A,B)|1'))
  expect_true(helper.found.edge(edge.df, '2_A|1', '3_A|1'))
})

test_that("merge_split_found_4", {
  dat <- data.table::data.table(
    SCA1=sample(seq(100,400), size=20),
    timepoint=c(rep(1,5), rep(2,5), rep(3,5), rep(4,5)),
    lineage=c(rep('A',5), 
              rep('A', 3), rep('A|1', 2), 
              rep('A', 2), rep('A|1', 2), rep('A|1|1',1),
              rep('A', 3), rep('(A|1,A|1|1)', 2))
  )
  edge.df <- GetTransitionsAsEdges(dat, c(1,2,3,4), 'timepoint', 'lineage')
  expect_true(helper.found.edge(edge.df, '1_A', '2_A'))
  expect_true(helper.found.edge(edge.df, '1_A', '2_A|1'))
  expect_true(helper.found.edge(edge.df, '2_A', '3_A'))
  expect_true(helper.found.edge(edge.df, '2_A|1', '3_A|1'))
  expect_true(helper.found.edge(edge.df, '2_A|1', '3_A|1|1'))
  expect_true(helper.found.edge(edge.df, '3_A', '4_A'))
  expect_true(helper.found.edge(edge.df, '3_A|1', '4_(A|1,A|1|1)'))
  expect_true(helper.found.edge(edge.df, '3_A|1|1', '4_(A|1,A|1|1)'))
})


test_that("repeated_alphabets_1", {
  dat <- data.table::data.table(
    SCA1=sample(seq(100,400), size=10),
    timepoint=c(rep(1,5), rep(2,5)),
    lineage=c(rep('A',5), rep('A', 3), rep('AA', 2))
  )
  edge.df <- GetTransitionsAsEdges(dat, c(1,2), 'timepoint', 'lineage')
  expect_true(helper.found.edge(edge.df, '1_A', '2_A'))
  expect_true(helper.not.found.edge(edge.df, '1_A', '2_AA'))
})

test_that("repeated_alphabets_2", {
  dat <- data.table::data.table(
    SCA1=sample(seq(100,400), size=15),
    timepoint=c(rep(1,5), rep(2,5), rep(3,5)),
    lineage=c(rep('A',5), 
              rep('A', 3), rep('AA', 2),
              rep('A',3), rep('AA', 1), rep('AA|1', 1))
  )
  edge.df <- GetTransitionsAsEdges(dat, c(1,2,3), 'timepoint', 'lineage')
  expect_true(helper.found.edge(edge.df, '1_A', '2_A'))
  expect_true(helper.found.edge(edge.df, '2_A', '3_A'))
  expect_true(helper.found.edge(edge.df, '2_AA', '3_AA'))
  expect_true(helper.found.edge(edge.df, '2_AA', '3_AA|1'))
  
  expect_true(helper.not.found.edge(edge.df, '1_A', '2_AA'))
  expect_true(helper.not.found.edge(edge.df, '2_A', '3_AA'))
  expect_true(helper.not.found.edge(edge.df, '2_A', '3_AA|1'))
  expect_true(helper.not.found.edge(edge.df, '2_AA', '3_A'))
})

test_that("repeated_alphabets_3", {
  dat <- data.table::data.table(
    SCA1=sample(seq(100,400), size=10),
    timepoint=c(rep(1,5), rep(2,5)),
    lineage=c(rep('A',3), rep('AA',2), 
              rep('(A,AA)', 3), rep('AAA', 2))
  )
  edge.df <- GetTransitionsAsEdges(dat, c(1,2), 'timepoint', 'lineage')
  expect_true(helper.found.edge(edge.df, '1_A', '2_(A,AA)'))
  expect_true(helper.found.edge(edge.df, '1_AA', '2_(A,AA)'))
  
  expect_true(helper.not.found.edge(edge.df, '1_AA', '2_AAA'))
  expect_true(helper.not.found.edge(edge.df, '1_A', '2_AAA'))
})