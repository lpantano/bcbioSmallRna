library(tidyverse)
library(janitor)
left_join(
    read_tsv("inst/extra/PRJEB3365.txt") %>% clean_names() %>% 
        select(secondary_sample_accession, fastq_ftp),
    read_tsv("inst/extra/E-GEUV-2.sdrf.txt") %>% clean_names() %>% 
        select(comment_ena_sample, factor_value_population, factor_value_laboratory),
    by = c("secondary_sample_accession" = "comment_ena_sample")) %>% 
    group_by(factor_value_population) %>% 
    sample_frac(0.02) %>% 
    select(samplenames=fastq_ftp,
           description=secondary_sample_accession,
           population=factor_value_population,
           laboratory=factor_value_laboratory) %>% 
    mutate(samplenames=paste0("ftp://", samplenames)) %>% 
    write_csv("inst/extra/geu_tiny.csv")

library(bcbioSmallRna)
sbcb = loadSmallRnaRun("inst/extra/geu_tiny/final/2018-09-29_geu_tiny", interestingGroups = "population")
save(sbcb, file = "data/sbcb.rda", compress = "xz")
