# Tufte Latex
marginwidth <- 1.94525 #in
textwidth <- 4.21342
pagewidth <- 6.29707
font <- "Source Sans Pro"

# ggplot scalings and point sizes
scal <- 1.5
pt <- 2
sz <- .25

pyr <- c("1,13-Tetradecadiene" = "14:2(1,13)",
         "1,14-Pentadecadiene" = "15:2(1,14)",
         "1,15-Hexadecadiene" = "16:2(1,15)",
         "1,16-Heptadecadiene" = "17:2(1,16)",
         "1,17-Octadecadiene" = "18:2(1,17)",
         "1,18-Nonadecadiene" = "19:2(1,18)",
         "1,19-Icosadiene" = "20:2(1,19)",
         "1,20-Henicosadiene" = "21:2(1,20)",
         "1,21-Docosadiene" = "22:2(1,21)",
         "1,22-Tricosadiene" = "23:2(1,22)",
         "2,4-Dimethyl-1-heptene" = "2,4Me9:1(1)",
         "Styrene" = "Sty",
         "α-Methylstyrene" = parse(text = "alpha*'MeSty'"))

poly <- c("Polyethylene" = "PE",
          "Polypropylene" = "PP",
          "Polystyrene" = "PS")

markers <- c("1,21-Docosadiene", "2,4-Dimethyl-1-heptene", "Styrene")

polymer_labeller <- as_labeller(function(string) {
  ifelse(is.na(poly[string]), string, poly[string])
})
