GCxGCanalyze <- function(Gari, min.threshold, unknowns, col_bleed, reduce.matrix, unlikel, 
                         group, RTmax) {options(warn=-1)
  if (group != TRUE) {ad <- Gari} else {
    ab <- interaction(Gari$Peak)
    ad <- cbind(Gari,id=factor(ab,levels=1:length(unique(ab))))
    ad$Cartridge <- NULL
    ac <- chron(times = Gari$rt1)
    ac <- chron(times = ad$rt1)
    for (i in 1:length(unique(ad$id))) {
      if (length(which(ad$id == i)) >= 2) {
        #ll <- rownames(ad[ad$id == i,])
        #ah <- which(ll == "NA")
        #ll <- ll[-ah]
        ff <- ad[ad$id == i,]
        hh <- ac[ad$id == i]
        d <- abs(outer(hh,hh,"-"))
        #diag(d) <- NA
        d[lower.tri(d, diag = TRUE)] <- 1
        d[is.na(d)] <- 1
        if (any(d <= (RTmax/84600))) {
          gh <- as.numeric(which(d <= (RTmax/86400), arr.ind=TRUE))
          jj <- colSums(ff[gh, 5:(ncol(ff)-1)])
          ff[-gh,1] <- paste("X", ff[-gh,1], sep = "_")
          fff <- as.character(ff[gh[1],])
          names(fff) <- colnames(ff[gh[1],])
          jjj <- c(fff[1:4], jj, fff[ncol(fff)])
          #jjj <- as.data.frame(t(do.call(rbind, jjj)))
          if (length(jjj) == 13) {jjj <- c(jjj, 0)}
          names(jjj[14]) <- "id"
          ff[gh,5:(ncol(ff)-1)] <- data.frame(t(jj))
          ad[ad$id == i,] <- ff
          ad$id[is.na(ad$id)] <- i
        } else {next()}
      }
    }
    if (!is.null(Gari$Cartridge)) {
      Cartridge <- Gari$Cartridge
      ad <- cbind(ad, Cartridge)
    } else {ad <- ad}
    ad$Peak <- tolower(ad$Peak)
    ad <- ad[!duplicated(ad$Peak, ad$id),]
    ad$Peak <- gsub("x_", "", ad$Peak)
  }
  excl_list1 <- c("SILOXANE", "glycol", "Silane", "silane", "Glycol",
                  "Xylil", "Xylyl", "xylil", "xylyl")
  excl_list2 <- c("CHLOR", "chlor", "FLUOR", "fluor", "BROMO", "bromo",
                  "METHANE, TRICHLORO-", "Tungsten", "tungsten",
                  "acetonitrile", "ACETONITRILE")
  excl_list3 <- c("Analyte", "analyte")
  if (unknowns == TRUE & col_bleed == TRUE & unlikel == TRUE) {
    excl_list <- c(excl_list1, excl_list2, excl_list3)
  }
  if (unknowns == TRUE & col_bleed == TRUE & unlikel == FALSE) {
    excl_list <- c(excl_list1, excl_list3)
  }
  if (unknowns == TRUE & col_bleed == FALSE & unlikel == FALSE) {
    excl_list <- c(excl_list3)
  }
  if (any(c(unknowns, col_bleed, unlikel) == TRUE)) {
    ut <- paste(excl_list, collapse="|")
    matches <- grep(ut, ad$Peak)
    ad <- ad[-matches,]
  } else {ad <- ad}
  if (reduce.matrix == TRUE) {
    subi <- length(grep("Blank", colnames(Gar))) + length(grep("BLANK", colnames(Gar)))
    oj <- length(which(lapply(ad, is.numeric) != TRUE))
    maxi <- apply(ad[,oj:(oj+leni-subi-1)], 1, max)
    ad <- ad[-which(maxi <= min.threshold),]
  }
  return(ad)
}
