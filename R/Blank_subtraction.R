Blank_subtraction <- function(Gar, Blank, LOD, reduce.matrix) {
  nami <- tolower(colnames(Gar))
  if ((length(grep("blank", nami)) == 0)) {Gar <- Gar} else {
  Gar$Peak <- tolower(Gar$Peak)
  Gar[is.na(Gar)] <- 0
  Sub_Gar <- function (Gar) {
    subi <- length(grep("Blank", colnames(Gar))) + length(grep("BLANK", colnames(Gar)))
    if (is.null(Gar$Mass) == FALSE) {
      white <- which((rowMeans(Gar[,5:(5+leni-subi)])-Gar$BLANK)/(Gar$BLANK*LOD) <= 1)
    } else {
      white <- which((rowMeans(Gar[,4:(4+leni-subi)])-Gar$BLANK)/(Gar$BLANK*LOD) <= 1)
    }
    return(white)
  }
  if (Blank == FALSE) {
    Gari <- Gar
    warning("Blank set to false")
  } else {
    if (length(grep("Blank", names(Gar)))) {
      white <- Sub_Gar(Gar)
    } else {
      if (Blank == FALSE){stop("Blank samples are missing")} else {
        #if (is.null(Gar$BLANK) == TRUE) {
        Gar$BLANK <- rowSums(Gar[,grepl("Blank", names(Gar))])
        white <- Sub_Gar(Gar)
      }
      if (length(white) == 0 & Blank != TRUE) {Gar <- Gar} else {Gar <- Gar[-white,]}
    }
  }
  }
  return(Gar)
}