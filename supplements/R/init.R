wd <- getwd()
setwd(dirname(sys.frame(1)$ofile))

source("packages.R")
source("variables.R")
source("wrappers.R")
source("read_openchrom.R")
source("read_targets.R")

setwd(wd)
rm(wd)
