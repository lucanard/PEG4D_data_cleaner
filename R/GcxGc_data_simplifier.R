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
#' @return a processed peaktable having all the real pekas present in your data in csv file in your working directory
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
  Gar_builder <- function(x, samp.nam) {
    options(warn=-1)
    if ((colnames(x[1]) != "Peak") == TRUE) {
      hiti <- x
      colnames(hiti) <- x[1,]
      hiti <- hiti[-1,]
      heiti <- hiti[,grep("Area", colnames(hiti))]
      names <- colnames(x)[grep(samp.nam, colnames(x))]
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
    } else {Gar <- x}
    return(Gar)
  }
  missing_eliminator <- function(x, samp.nam) {
    ab <- grep(samp.nam, colnames(x))
    ac <- which(apply(x[,ab], 2, median) == 0)
    
    if ((length(ac) != 0) == TRUE) {x1 <- x[,-ab[ac]]} else {x1 <- x}
    return(x1)
  }
  Sub_Gar <- function (Gar) {
    subi <- length(grep("Blank", colnames(Gar))) + length(grep("BLANK", colnames(Gar)))
    if (is.null(Gar$Mass) == FALSE) {
      white <- which((apply(Gar[,5:(5+leni-subi)], 1, median)-Gar$BLANK)/(Gar$BLANK*LOD) <= 1)
    } else {
      white <- which((apply(Gar[,4:(4+leni-subi)], 1, median)-Gar$BLANK)/(Gar$BLANK*LOD) <= 1)
    }
    return(white)
  }
  Blank_subtraction <- function(Gar, LOD) {
    nami <- tolower(colnames(Gar))
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
      subi <- length(grep("Blank", colnames(Gari))) + length(grep("BLANK", colnames(Gari)))
      oj <- length(which(lapply(ad, is.numeric) != TRUE))
      maxi <- apply(ad[,oj:(oj+leni-subi-1)], 1, max)
      ad <- ad[-which(maxi <= min.threshold),]
    }
    return(ad)
  }
  Gar <- Gar_builder(x, samp.nam)
  Gar$Peak <- tolower(Gar$Peak)
  Gar[is.na(Gar)] <- 0
  Gar <- missing_eliminator(Gar, samp.nam)
  #clean$Peak <- NULL
  #no_Gar <- Gar[,-grep(samp.nam, colnames(Gar))]
  #blank <- no_Gar[,grep("Blank", colnames(no_Gar))]
  #no_blank <- no_Gar[,-grep("Blank", colnames(no_Gar))]
  if (is.null(Gar$Mass) == FALSE) {ms <- 1} else {ms <- 0}
  leni <- ncol(Gar) - length(which(lapply(Gar, is.numeric) != TRUE)) - ms
  if (Blank == TRUE) {
    Gari <- Blank_subtraction(Gar, LOD)
  } else {
    Gari <- Gar
    warning("Blank set to false")
    }
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
  write.csv(ad, file = "simple_GCxGC_table.csv")
  return(ad)
}
