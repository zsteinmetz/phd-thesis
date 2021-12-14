#### Load packages and functions ####
source("func/init.R")

pyrograms <- rbind(
  cbind(`Sample Type` = "Blank", fread("../data/chromatograms/Pyro_Plast-015.csv")),
  cbind(`Sample Type` = "Polymer mixture", fread("../data/chromatograms/Pyro_Plast-026.csv")),
  cbind(`Sample Type` = "Pure PE", fread("../data/chromatograms/Pyro_Plast-027.csv")),
  cbind(`Sample Type` = "Pure PP", fread("../data/chromatograms/Pyro_Plast-028.csv")),
  cbind(`Sample Type` = "Pure PS", fread("../data/chromatograms/Pyro_Plast-029.csv"))
)

chrom <- data.table(read_openchrom("../data/reports/sample/"))
chrom <- chrom[!is.na(Name), .(Name, RI, `Integrated Area`)]

abbr <- fread("../../../tables/py-products.tex", sep = "&")
abbr <- abbr[, .(V2, V3)]
names(abbr) <- c("Label", "Name")
abbr[Name == "$\\alpha$-Methylstyrene", Name := "α-Methylstyrene"]
abbr[Label == "$\\alpha$MeSty", Label := "αMeSty"]
abbr[Label == "2,4,6Me12:1(1)\\textit{i}", Label := "2,4,6Me12:1(1)i"]
abbr[Label == "2,4,6Me12:1(1)\\textit{h}", Label := "2,4,6Me12:1(1)h"]

labels <- merge(chrom[RI %between% c(820, 2000)], abbr, by = "Name")
labels[, rTIC := `Integrated Area`]
labels[RI %between% c(885, 960), rTIC := `Integrated Area`/60]
labels[, RI := RI + 2] # correct for column aging

pyrograms <- pyrograms[RI %between% c(820, 2000)]
# Calculate total ion current
pyrograms[, TIC := Reduce(`+`, .SD), .SDcols = 5:length(pyrograms)][]

pyrograms[, rTIC := TIC]
# Scale down styrene
pyrograms[RI %between% c(885, 960), rTIC := TIC/60]

pyrograms[RI %between% c(833, 855), Label := "2,4Me9:1(1)"]
pyrograms[RI %between% c(880, 902), Label := "Sty"]
pyrograms[RI %between% c(976, 998), Label := "αMeSty"]
pyrograms[RI %between% c(1480, 1502), Label := "15:2(1,14)"]
pyrograms[RI %between% c(1680, 1702), Label := "17:2(1,16)"]
pyrograms[RI %between% c(1780, 1802), Label := "18:2(1,17)"]

pyrograms[, Rel := TIC/max(TIC), by = Label]
pyrograms[, RI := RI + 2] # correct for column aging

pyrograms[, `Sample Type`:= factor(`Sample Type`, levels = unique(`Sample Type`))]
pyrograms[, Label := factor(Label,
                            levels = c("15:2(1,14)", "17:2(1,16)", "18:2(1,17)",
                                      "2,4Me9:1(1)", "Sty", "αMeSty", NA))]

ggplot(pyrograms[`Sample Type` == "Polymer mixture"], aes(RI, rTIC)) +
  geom_text_repel(data = labels, aes(label = Label), size = 3.5, alpha = .66,
                  family = font,
                  nudge_y = 950000, 
                  point.padding = 0.25,
                  segment.size = 0.25) +
  geom_line(size = 0.4, color = viridis(1)) +
  xlab("Retention index (RI)") +
  ylim(c(0, 1.2 * max(pyrograms$rTIC))) +
  theme_publish(base_family = font) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave('../../../figures/py-sample.pdf', scale = 1.5, width = pagewidth,
       height = 3, unit = 'in', device = cairo_pdf)

ggplot(pyrograms[!is.na(Label)], aes(RI, Rel*100)) +
  geom_line(aes(color = `Sample Type`, linetype = `Sample Type`)) +
  geom_text_repel(data = labels[Label %in% unique(pyrograms$Label)],
                  aes(label = Label), y = 100, size = 3.5, alpha = .66,
                  family = font,
                  nudge_y = 30, nudge_x = -1,
                  point.padding = 0.25,
                  segment.size = 0.25) +
  facet_wrap(~ Label, scales = "free", labeller =
               as_labeller(function(string){c("PE", "PE", "PE", "PP", "PS", "PS")})) +
  xlab("Retention index (RI)") +
  scale_color_manual(values = c("black", "black", viridis(3))) +
  scale_linetype_manual(values = c("dotted", "dashed", rep("solid", 3))) +
  scale_y_continuous(limits = c(0,150)) +
  theme_publish(base_family = font) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave('../../../figures/py-selectivity1.pdf', scale = 1.5, width = pagewidth,
       height = 3, unit = 'in', device = cairo_pdf)

