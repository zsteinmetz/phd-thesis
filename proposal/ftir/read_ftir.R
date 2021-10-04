read_ftir <- function(path = ".", dec = ".", encoding = getOption("encoding")) {
  files <- list.files(path = path, pattern = "\\.asp$")
  
  ret <- data.frame()
  for (cfile in files) {
    con <- file(file.path(path, cfile), encoding = encoding)
    lines <- as.numeric(readLines(con))
    close(con)
    
    sig <- lines[-c(1:6)]
    wav <- seq(lines[2], lines[3], length.out = lines[1])
    dat <- data.frame(Measurement = sub("\\.asp$", "", cfile),
                      Wavenumber = wav, Absorbance = sig)
    ret <- rbind(ret, dat)
  }
  return(ret)
}