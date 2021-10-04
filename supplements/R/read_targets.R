read_targets <- function(...) {
  library(data.table)
  targets <- fread(...)
  names(targets) <-
    c("RTstart", "RTstop", "Name", "CAS", "Comment", "Contributor", "Reference",
      "Traces", "Reference Identifier")[1:length(targets)]
  return(targets)
}
