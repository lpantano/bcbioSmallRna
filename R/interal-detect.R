#' Detect the organism from the genome build name
#'
#' @rdname detect_organism
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param genome_build Genome build.
#' @noRd
#' @return Organism string.
.detect_organism <- function(genome_build) {
    if (str_detect(genome_build, "^(hg|GRCh)\\d+")) {
        return("hsapiens")
    } else if (str_detect(genome_build, "^mm\\d+")) {
        return("mmusculus")
    } else if (str_detect(genome_build, "^rn\\d+")) {
        return("rnorvegicus")
    } else if (str_detect(genome_build, "^WBcel\\d+")) {
        return("celegans")
    } else if (str_detect(genome_build, "^BDGP\\d+")) {
        return("dmelanogaster")
    } else if (str_detect(genome_build, "^Zv\\d+")) {
        return("drerio")
    } else if (str_detect(genome_build, "^ASM\\d+")) {
        return("spombe")
    } else if (str_detect(genome_build, "^(MB|MG)\\d+")) {
        return("ecoli")
    } else {
        warning("Failed to detect organism")
    }
    return("NA")
}