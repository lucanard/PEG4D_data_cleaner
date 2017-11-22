QC_sub <- function(Gari) {
  Gari_nb <- Gari[,!grepl("Blank", names(Gari))]
  Gari_nb <- Gari_nb[,5: ncol(Gari_nb)]
  mean_QC <- rowMeans(Gari_nb[,grepl("QC", names(Gari_nb))])
  max_sample <- apply(Gari_nb[,!grepl("QC", names(Gari_nb))], 1, max)
  strani <- which(mean_sample - mean_QC <= 0)
  Gari <- Gari[-strani,]
  return(Gari)
}
QC_sub_spec <- function (Gari) {
  Gari_nb <- Gari[,!grepl("Blank", names(Gari))]
  Gari_nb <- Gari_nb[,5: ncol(Gari_nb)]
  Gari_nbq <- Gari_nb[,!grepl("QC", names(Gari_nb))]
  level1 <- levels(as.factor(str_sub(colnames(Gari[, grepl("QC", names(Gari))]), start = 1, end = 4)))
  meany <- sapply(level1, function (x) rowMeans(Gari_nb[,grep(x, names(Gari_nb))]))
  meany1 <- sapply(str_sub(level1, start = nchar(level1)), function (x) apply(Gari_nbq[,grep(x, names(Gari_nbq))], 1, max))
  diffy <- meany1 - meany
  diffy1 <- diffy <= 0
  liste <- lapply(seq(1, ncol(diffy1), 1), function(x) Gari[diffy1[,x] == FALSE,])
  return(liste)
}