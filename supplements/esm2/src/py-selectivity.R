#### Load packages and functions ####
source("func/init.R")

## Read data
seqtable <- fread("../data/seqtable.csv")

allcompounds <- rbind(
  cbind(read_targets("../data/targets/PE.txt", header = F)),
  cbind(read_targets("../data/targets/PP.txt", header = F)),
  cbind(read_targets("../data/targets/PS.txt", header = F))
)
names(allcompounds)[names(allcompounds) == "Contributor"] <- "Polymer"

comp <- list()

comp$reports <- rbind(
  data.table(read_openchrom("../data/reports/calibration/")),
  data.table(read_openchrom("../data/reports/selectivity/"))
)

## Data handling
comp$data <- merge(seqtable, comp$reports, by = "File Name")

comp$data[, Conc := as.numeric(gsub(" ug/mL.*", "", `Sample Name`))]

comp$data[Conc == 0, `Sample Type` := "Blank"]
comp$data[Conc == 150 & `Sample ID` != 21, `Sample Type` := paste("Pure", gsub(".*ug/mL ", "", `Sample Name`))]
comp$data[`Sample Type` == "Pure PlastMix", `Sample Type` := "Polymer mixture"]

comp$targets <- merge(comp$data[Comments == "quantifiable" & `Sample Type` != "Unknown"],
                      allcompounds, by = "Name")

comp$targets[, .N, by = .(`Sample Type`, Name)]
comp$targets[Name == "Styrene", .(`Sample ID`, `Sample Type`, `Integrated Area`)]

comp$targets[, `Relative Area` := `Integrated Area`/max(`Integrated Area`), by = Polymer]

## Stats
comp$targets[, .(Mean = mean(`Relative Area`),
                 RSD = sd(`Relative Area`)/mean(`Relative Area`)), by = .(`Sample Type`, Name)]

comp$targets[, SampleType := as.factor(`Sample Type`)]

comp$aov <- sapply(unique(comp$targets$Name), function(x) {
  mod <- aov(`Integrated Area` ~ SampleType, data = comp$targets[Name == x])
  op <- par(mfrow = c(1,2))
  plot(mod, which = 1:2)
  par(op)
  mod
}, simplify = FALSE)

sapply(comp$aov, summary, simplify = F)

comp$glht <-sapply(comp$aov, function(obj) {
  glht(obj, linfct = mcp(SampleType = "Tukey"))
}, simplify = FALSE)
  
sapply(comp$glht, summary, test = adjusted("bonferroni"), simplify = F)
comp$cld <- sapply(comp$glht, function(x) {
  cld(x, test = adjusted("bonferroni"), level = 0.05)$mcletters$Letters
  }, simplify = F)

comp$letters <- as.data.table(t(as.data.table(comp$cld)))
names(comp$letters) <- names(comp$cld[[1]])
comp$letters[, Name := names(comp$cld)]
comp$letters <- merge(
  melt(comp$letters, id.vars = "Name", value.name = "Letter", variable.name = "Sample Type"),
  comp$targets[, .(`Integrated Area` = max(`Integrated Area`)), by = .(Polymer, Name, `Sample Type`)],
  by = c("Name", "Sample Type"))

## Plotting
comp$plot <- comp$targets[!(Name %in% c("1,13-Tetradecadiene", "1,15-Hexadecadiene"))]
comp$letters <- comp$letters[!(Name %in% c("1,13-Tetradecadiene", "1,15-Hexadecadiene"))]

levs <- comp$plot[, unique(`Sample Type`)]
comp$plot[, `Sample Type`:= factor(`Sample Type`, levels = levs)]
comp$letters[, `Sample Type`:= factor(`Sample Type`, levels = levs)]

pd <- position_dodge(.7)

comp$gg <- list()
(comp$gg$pe <-
  ggplot(data = comp$plot[Polymer == "Polyethylene"],
                     aes(Name, `Integrated Area`*10^-6)) +
  geom_point(size = pt, position = pd,
             aes(color = `Sample Type`, shape = `Sample Type`, fill = `Sample Type`)) +
  geom_text(data = comp$letters[Polymer == "Polyethylene"], position = pd, size = 3,
            aes(Name, `Integrated Area`*10^-6 + 0.1 * max(`Integrated Area`*10^-6),
                label = paste0("italic(", Letter, ")"), group = `Sample Type`), parse = T) +
  xlab("") + ylab(expression("Peak area"~"["%.%~10^6*"]")) +
  facet_wrap(. ~ Polymer, scales = "free", labeller = polymer_labeller) +
  scale_color_manual(values = c("black", "black", viridis(3))) +
  scale_fill_manual(values = c("black", "black", viridis(3))) +
  scale_shape_manual(values = c(1, 19, 22, 23, 24)) +
  scale_x_discrete(labels = pyr) +
  ylim(0, NA) +
  theme_publish(base_family = font))
(comp$gg$pp <-
  ggplot(data = comp$plot[Polymer == "Polypropylene"],
         aes(Name, `Integrated Area`*10^-6)) +
  geom_point(size = pt, position = pd,
             aes(color = `Sample Type`, shape = `Sample Type`, fill = `Sample Type`)) +
  geom_text(data = comp$letters[Polymer == "Polypropylene"], position = pd, size = 3,
            aes(Name, `Integrated Area`*10^-6 + 0.1 * max(`Integrated Area`*10^-6),
                label = paste0("italic(", Letter, ")"), group = `Sample Type`), parse = T) +
  xlab("") + ylab(expression("Peak area"~"["%.%~10^6*"]")) +
  facet_wrap(. ~ Polymer, scales = "free", labeller = polymer_labeller) +
  scale_color_manual(values = c("black", "black", viridis(3))) +
  scale_fill_manual(values = c("black", "black", viridis(3))) +
  scale_shape_manual(values = c(1, 19, 22, 23, 24)) +
  scale_x_discrete(labels = pyr) +
  ylim(0, NA) +
  theme_publish(base_family = font) + theme(axis.title.y = element_blank()))
(comp$gg$ps <-
  ggplot(data = comp$plot[Polymer == "Polystyrene"],
         aes(Name, `Integrated Area`*10^-6)) +
  geom_point(size = pt, position = pd,
             aes(color = `Sample Type`, shape = `Sample Type`, fill = `Sample Type`)) +
  geom_text(data = comp$letters[Polymer == "Polystyrene"], position = pd, size = 3,
            aes(Name, `Integrated Area`*10^-6 + 0.1 * max(`Integrated Area`*10^-6),
                label = paste0("italic(", Letter, ")"), group = `Sample Type`), parse = T) +
  xlab("") + ylab(expression("Peak area"~"["%.%~10^6*"]")) +
  facet_wrap(. ~ Polymer, scales = "free", labeller = polymer_labeller) +
  scale_color_manual(values = c("black", "black", viridis(3))) +
  scale_fill_manual(values = c("black", "black", viridis(3))) +
  scale_shape_manual(values = c(1, 19, 22, 23, 24)) +
  scale_x_discrete(labels = pyr) +
  ylim(0, NA) +
  theme_publish(base_family = font) + theme(axis.title.y = element_blank()))

combined <- comp$gg$pe + comp$gg$pp + comp$gg$ps &
  theme(legend.position = "bottom")
combined + plot_layout(guides = "collect", widths = c(3, 1, 2)) +
  plot_annotation(theme = theme(plot.margin = margin()))
ggsave('../../../figures/py-selectivity2.pdf', scale = 1.5, width = pagewidth,
       height = 2.2, unit = 'in', device = cairo_pdf)
ggsave('../../../defense/py-selectivity2.png', scale = 1, width = 8,
       height = 3, unit = 'in', bg = "white")
