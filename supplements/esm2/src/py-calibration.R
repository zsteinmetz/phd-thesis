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

plast <- list()

## Data handling
plast$reports <- data.table(read_openchrom("../data/reports/calibration/"))

plast$data <- merge(seqtable, plast$reports, by = "File Name")
plast$data[, Conc := as.numeric(gsub(" ug/mL.*", "", `Sample Name`))]

## Calibrate
plast$all <- merge(plast$data, allcompounds, by = "Name")
plast$all[, batchcal(`Integrated Area` ~ Conc), by = .(Polymer, Name)]

plast$targets <- merge(plast$data[Comments == "quantifiable"],
                       allcompounds, by = "Name")

(plast$cal <- plast$targets[, batchcal(`Integrated Area` ~ Conc), by = .(Polymer, Name)])

plast$plot <- plast$targets[Comments == "quantifiable" &
                              !(Name %in% c("1,13-Tetradecadiene", "1,15-Hexadecadiene"))]
plast$pred <- plast$plot[, batchpred(`Integrated Area` ~ Conc), by = .(Polymer, Name)]

## Plotting
ggplot(data = plast$plot, aes(Conc, `Integrated Area`*10^-6)) +
  geom_ribbon(data = plast$pred, alpha = .25,
              aes(y = fit*10^-6, ymin = lwr*10^-6, ymax = upr*10^-6, fill = Name)) +
  geom_line(data = plast$pred, aes(y = fit*10^-6, color = Name)) +
  geom_point(size = pt, aes(color = Name, shape = Name)) +
  facet_wrap(~ Polymer, scales = "free", labeller = polymer_labeller) +
  xlab(expression("Polymer concentration"~"["*"\U003BCg"~mL^-1*"]")) +
  ylab(expression("Peak area"~"["%.%~10^6*"]")) +
  scale_shape_manual(name = "Pyrolysates", values = c(1, 17, 15, 6, 23, 19), labels = pyr) +
  scale_color_viridis(discrete = T, name = "Pyrolysates", labels = pyr) +
  scale_fill_viridis(discrete = T, name = "Pyrolysates", labels = pyr) +
  theme_publish(base_family = font) + theme(legend.text.align = 0)
ggsave('../../../figures/py-calibration.pdf', scale = 1.5, width = pagewidth,
       height = 2.5, unit = 'in', device = cairo_pdf)

ggplot(data = plast$plot[Name %in% c("Styrene", "2,4-Dimethyl-1-heptene",
                                     "1,17-Octadecadiene")],
       aes(Conc, `Integrated Area`*10^-6)) +
  geom_ribbon(data = plast$pred[Name %in% c("Styrene", "2,4-Dimethyl-1-heptene",
                                            "1,17-Octadecadiene")], alpha = .25,
              aes(y = fit*10^-6, ymin = lwr*10^-6, ymax = upr*10^-6, fill = Name)) +
  geom_line(data = plast$pred[Name %in% c("Styrene", "2,4-Dimethyl-1-heptene",
                                          "1,17-Octadecadiene")],
            aes(y = fit*10^-6, color = Name)) +
  geom_point(size = pt, aes(color = Name, shape = Name)) +
  facet_wrap(~ Polymer, scales = "free", labeller = polymer_labeller) +
  xlab(expression("Polymer concentration"~"["*"\U003BCg"~mL^-1*"]")) +
  ylab(expression("Peak area"~"["%.%~10^6*"]")) +
  scale_shape_manual(name = "Pyrolysates", values = c(17, 15, 19)) +
  scale_color_viridis(discrete = T, name = "Pyrolysates") +
  scale_fill_viridis(discrete = T, name = "Pyrolysates") +
  theme_publish(base_family = font) + theme(legend.position = "none")
ggsave('../../../defense/py-calibration.png', scale = 1, width = 8,
       height = 3, unit = 'in', bg = "white")
