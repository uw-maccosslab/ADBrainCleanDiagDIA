#' Process source data to raw
#'
#' @param dat source peak area output from Skyline
#'
#' @return formatted raw dataframe
#' @export
source_to_raw <- function(dat) {

    colnames(dat) = sapply(strsplit(colnames(dat), " "), "[", 1)
    colnames(dat) = sapply(strsplit(colnames(dat), "-"), "[", 1)
    if ("PeptideModifiedSequence" %in% colnames(dat) &
        "Peptide" %in% colnames(dat))
      dat = dat[,-which(colnames(dat)=="Peptide")]
    if (length(grep("PeptideModifiedSequence", colnames(dat))) > 0)
      colnames(dat)[which(colnames(dat) == "PeptideModifiedSequence")] = "Peptide"
    if (length(grep("ProteinName", colnames(dat))) > 0)
      colnames(dat)[which(colnames(dat) == "ProteinName")] = "Protein"
    uniquePeptide = paste0(dat$Peptide, "@", dat$Protein)
    dat = dat[,-c(1:2)]
    rownames(dat) = uniquePeptide

  return(dat)
}
