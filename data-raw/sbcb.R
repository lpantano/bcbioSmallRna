library(bcbioSmallRna)

sbcb = loadSmallRnaRun("inst/extra/bcbio/2018-02-21_samples", interestingGroups = "country")
save(sbcb, file = "data/sbcb.rda", compress = "xz")
