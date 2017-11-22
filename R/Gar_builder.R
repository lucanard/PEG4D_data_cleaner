Gar_builder <- function(prova, samp.nam) {
  options(warn=-1)
  if ((colnames(prova[1]) != "Peak") == TRUE) {
    hiti <- prova
    colnames(hiti) <- prova[1,]
    hiti <- hiti[-1,]
    heiti <- hiti[,grep("Area", colnames(hiti))]
    names <- colnames(prova)[grep(samp.nam, colnames(prova))]
    colnames(heiti) <- c(names, rep("Blank", (length(heiti) - length(names))))
    heitix <- as.data.frame(lapply(heiti, function(x) as.numeric(as.character(x))))
    rt1s <- hiti[,grep("1st", colnames(hiti))]
    rt2s <- hiti[,grep("2nd", colnames(hiti))]
    RTis <- function(rts) {
      rts1 <- lapply(rts, function(x) as.numeric(as.period(ms(x), unit = "sec")))
      RT <- as.Date(seconds_to_period(rowMeans(as.data.frame(rts1), na.rm = T)),
                    format = "%S/%oS", origin = Sys.Date())
      RT <- gsub(RT, pattern = Sys.Date(), replacement = "", fixed = T)
      RT <- gsub(RT, pattern = "UTC", replacement = "", fixed = T)
      RT <- gsub(RT, pattern = " ", replacement = "", fixed = T)
      return(RT)
    }
    rt1 <- RTis(rt1s)
    rt2 <- RTis(rt2s)
    Peak <- as.character(unlist(hiti[,1]))
    Gar <- cbind(Peak, rt1, rt2, heitix)
  } else {Gar <- prova}
  return(Gar)
}