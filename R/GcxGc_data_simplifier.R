#' @title Cleaning of the peaktables from the PEG4D Leco instrument
#' @description The function transforms the raw peaktable obtained from the chromatof deconvolution software to a clean pekatable with only the names of the Blanks, QCs, samples, Retention times and compound names.
#' @param x the peaktable obtained from the chromatof deconvolution software
#' @param Blanks should be the peaks delete according to their presence in blanks? If TRUE, every compound present also in the blanks with an average amount lower in the samples lower than the average amount found in the blanks multiplied by the LOD factor. default = TRUE
#' @param limit_of_detection the limit of detection for the compounds found also in the blanks. Default = 3
#' @param sample_name the identificative letter/shortname of the samples. Mandatory to allow the algorithm to recognize the samples from the QCs.
#' @param simplify_matrix if TRUE it will delete the compounds having a maximum value below the min.threshold value. default = TRUE
#' @param minimum.threshold the minimum threshold applied to reduce the matrix according to the maximum intensity found in any of the samples
#' @param unknown should unknown peaks be excluded? default set to TRUE
#' @param column_bleed should common column bleedings be excluded? default set to TRUE
#' @param unlikely should compounds containing unlikely atoms (chlorine, fluorine etc.) be excluded? default set to TRUE
#' @param grouping should compounds having the same name and similar retention times be grouped? Default set to TRUE
#' @param RTspan the maximum distance in seconds between to peaks with the same name to be grouped.
#' @return a processed peaktable having all the real pekas present in your data
#' @export "GcxGc_data_cleaning"
#' @author Luca Narduzzi "nardluca@gmail.com"
#' @examples
#' Gar <- read.csv(system.file("extdata", "Gar.csv", package = "GCxGC.Leco.analyzer"), row.names = 1, stringsAsFactors = FALSE)
#' GcxGCdata <- GcxGc_data_cleaning(Gar)
GcxGc_data_cleaning <- function(x, sample_name, Blanks = TRUE, limit_of_detection = 3, simplify_matrix = TRUE, minimum.threshold = 120000,
                                unknown = TRUE, column_bleed = TRUE, unlikely = TRUE, grouping = TRUE, RTspan = 40) {
  if (missing(sample_name)) {
    samp.nam <- names(table(str_extract(colnames(x), "[a-z]+")))
    } else {samp.nam <- sample_name} 
  Blank <- Blanks
  LOD <- limit_of_detection
  reduce.matrix <- simplify_matrix
  min.threshold <- minimum.threshold
  unknowns <- unknown
  col_bleed <- column_bleed
  unlikel <- unlikely
  group <- grouping
  RTmax <- RTspan
  Gar <- Gar_builder(x, samp.nam)
  if (is.null(Gar$Mass) == FALSE) {ms <- 1} else {ms <- 0}
  leni <- ncol(Gar) - length(which(lapply(Gar, is.numeric) != TRUE)) - ms
  if (Blank == TRUE) {
    Gari <- Blank_subtraction(Gar, Blank, LOD, reduce.matrix)
  } else {Gari <- Gar}
  if (any(grepl("QC", colnames(Gari)) == TRUE) == TRUE) {
    if (any(grep("^QC$", colnames(Gari)))){
      Gari <- QC_sub(Gari)
    } else {
      Gari <- QC_sub_spec(Gari)
    }
  } else {Gari = Gari}
  if (class(Gari) == "data.frame") {
    ad <- GCxGCanalyze(Gari, min.threshold, unknowns, col_bleed, reduce.matrix,
                       unlikel, group, RTmax)
  } else {
    ad <- lapply(seq(1, length(Gari), 1), function(x) GCxGCanalyze(Gari[[x]], min.threshold, unknowns, col_bleed, reduce.matrix, unlikel, group, RTmax))
  }
  return(ad)
}
