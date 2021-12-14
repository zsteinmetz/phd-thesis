read_openchrom <- function(
  path = ".", sep = "\t", dec = ".", encoding = getOption("encoding"),
  verbose = FALSE) {

  files <- list.files(path = path, pattern = "\\.txt$")

  ret <- data.frame()
  for (cfile in files) {
    if (verbose) message("reading ", cfile)
    
    con <- file(file.path(path, cfile), encoding = encoding)
    lines <- readLines(con)
    close(con)

    pls <- which(lines == "PEAK LIST OVERVIEW") + 3
    nl <- which(lines == "")
    ple <- nl[nl > pls][1] - 1

    peaks <- read.table(text = lines[pls:ple], sep = sep, dec = dec,
                        header = F)
    names(peaks) <- strsplit(lines[pls-1], "\t")[[1]]

    peaks <- peaks[!duplicated(names(peaks))]

    ids <- which(lines == "PEAK IDENTIFICATION RESULTS") + 3
    nl <- which(lines == "")
    ide <- nl[nl > ids][1] - 1
    
    sds <- which(lines == "SCAN IDENTIFICATION RESULTS") + 3
    nl <- which(lines == "")
    sde <- nl[nl > sds][1] - 1
    
    txt <- c(lines[ids:ide][!(lines[ids:ide] %in% c("---", ""))],
             lines[sds:sde][!(lines[sds:sde] %in% c("---", ""))])
    
    if (length(txt)) {
      ident <- read.table(text = txt, sep = sep, dec = dec, header = F)
      names(ident) <- strsplit(lines[ids-1], "\t")[[1]]
      
      peaks <- peaks[!duplicated(names(peaks))]
      
      fname <- gsub("^Name: ", "", unique(lines[grepl("^Name: ", lines)]))
      
      comps <- base::merge(peaks, ident, by = "[#]", all = T)
      
      comps <- cbind(`File Name` = fname, comps)
      ret <- rbind(ret, comps)
    } else {
     warning("skipping ", cfile, "; empty file")
     next
    }
  }
  return(ret)
}
