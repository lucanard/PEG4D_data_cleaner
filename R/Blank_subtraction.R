Blank_subtraction <- function(Gar, LOD) {
  nami <- tolower(colnames(Gar))
  Sub_Gar <- function (Gar) {
    subi <- length(grep("Blank", colnames(Gar))) + length(grep("BLANK", colnames(Gar)))
    if (is.null(Gar$Mass) == FALSE) {
      white <- which((rowMeans(Gar[,5:(5+leni-subi)])-Gar$BLANK)/(Gar$BLANK*LOD) <= 1)
    } else {
      white <- which((rowMeans(Gar[,4:(4+leni-subi)])-Gar$BLANK)/(Gar$BLANK*LOD) <= 1)
    }
    return(white)
  }
  if ((length(grep("blank", nami)) == 0)) {
    Gar <- Gar
    warning("Blanks are missing")
  } else {
    #Gar$Peak <- tolower(Gar$Peak)
    #Gar[is.na(Gar)] <- 0
    if (length(grep("Blank", names(Gar))) == 1) {
      white <- Sub_Gar(Gar)
    } else {
      hh <- grepl("Blank", colnames(Gar))
      Gar$BLANK <- apply(Gar[,hh], 1, median)       
      #Gar$BLANK <- rowMeans(Gar[,grepl("Blank", names(Gar))])
      white <- Sub_Gar(Gar)
    }
    if (length(white) == 0) {Gar <- Gar} else {Gar <- Gar[-white,]}
  }
  return(Gar)
}