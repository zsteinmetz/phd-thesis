batchcal <- function(...) {
  tryCatch(
    {
      cal <- calibration(...)
      # plot(cal)
      list(Int = cal$model$coefficients[1],
           Slope = cal$model$coefficients[2],
           Radj = cal$adj.r.squared,
           LOD = cal$lod[1], LOQ = cal$loq[1],
           BlankMean = mean(cal$blanks),
           BlankSD = sd(cal$blanks))
    },
    error = function(e) {
      warning(e)
      list(Int = NaN,
           Slope = NaN,
           Radj = NaN,
           LOD = NaN, LOQ = NaN,
           BlankMean = NaN,
           BlankSD = NaN)
    }
  )
}

batchpred <- function(...) {
  cal <- calibration(...)
  
  model <- cal$model
  
  conc <- model$model[,2]
  new <- data.frame(conc = seq(min(conc), max(conc), length.out = 100 * length(conc)))
  names(new) <- all.vars(model$formula)[2]
  
  data.frame(new, predict(cal$model, new, interval = "conf"))
}

multirow <- function(x) {
  u <- sapply(unique(x), match, x)
  x[-u] <- ""
  return(x)
}

theme_red <- function() {
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
}

scientific_10 <- function(x) {
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
}