test_that("report", {
    
    dir.create("report")
    setwd("report")
    data(sbcb)
    save(sbcb, file = "sbcb.rda")
    file.copy(file.path(system.file(package = "bcbioSmallRna"),
                        "rmarkdown",
                        "templates",
                        "srnaseq",
                        "skeleton",
                        "skeleton.Rmd"), "test.Rmd")
    rmarkdown::render("test.Rmd", params = list(bcb = "sbcb.rda"))
    setwd("..")
    unlink("report", recursive = TRUE)
    
})