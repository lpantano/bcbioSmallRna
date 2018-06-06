library(bcbioSmallRna)

sbcb = loadSmallRnaRun("inst/extra/bcbio", interestingGroups = "country")
save(sbcb, file = "data/sbcb.rda", compress = "xz")
