#### Load packages and global vars ####
source("func/init.R")

#### Read data ####
seqtable <- fread("../data/py-gc-ms/seqtable.csv")

allcompounds <- rbind(
  read_targets("../data/py-gc-ms/targets/Ref.txt", header = F),
  read_targets("../data/py-gc-ms/targets/PS-d5.txt", header = F),
  read_targets("../data/py-gc-ms/targets/PS.txt", header = F),
  read_targets("../data/py-gc-ms/targets/PP.txt", header = F),
  read_targets("../data/py-gc-ms/targets/PE.txt", header = F)
)
names(allcompounds)[names(allcompounds) == "Contributor"] <- "Polymer"

mama <- list()

mama$reports <- read_openchrom("../data/py-gc-ms/reports/matrix-match/") %>% 
  data.table() %>% 
  rbind()

names(mama$reports)[names(mama$reports) == "Comments"] <- "Extra"

#### Arrange measurement data ####
mama$data <- merge(seqtable, mama$reports, by = "File Name")

mama$data[grep("ug/mL", `Sample Name`, ignore.case = T),
          `Sample Type` := "Standard"]
mama$data[`Sample Type` == "Standard",
               Conc := as.numeric(gsub(" ug/mL.*", "",
                                       `Sample Name`, ignore.case = T))]
mama$data[grepl("PS-d5", `Sample Name`), `Cal Type` := "In solvent mixture"]
mama$data[grepl("LUFA", `Sample Name`), `Cal Type` := "In LUFA 2.2 matrix"]

mama$data[, Acquisition := as.POSIXct(Date)]
mama$data[, CW := isoweek(Acquisition)]

mama$all <- merge(mama$data,
                  allcompounds[, .(Name, Polymer, `Reference Identifier`)],
                  by = "Name")

# Check for contamination and other problems
mama$all[!(Comment %in% c("", "Repeatability")) & Name == "Styrene-d5",
         .(`Sample ID`, `Sample Name`, Comment)]

# Remove measurements
mama$targets <- mama$all[Extra == "quantifiable" &
                           !(`Sample ID` %in% c(899))]

#### Repeatability ####
# Check for internal repeatability per day
mama$targets[grepl("*-d5$", Name),
             .(RSD = sd(`Integrated Area`) / mean(`Integrated Area`)),
             by = .(Acquisition, Name, `Cal Type`)]

ggplot(mama$targets[grepl("*-d5$", Name)],
       aes(`Sample ID`, `Integrated Area`)) +
  geom_point() +
  geom_text(aes(label = `Sample ID`), size = 2.5, nudge_x = 2) +
  facet_wrap(. ~ Name, scales = "free_y") + theme_publish()

# Internal standard correction?
mama$targets[, `IS corrected Area` := `Integrated Area` /
               `Integrated Area`[Name == "Styrene-d5"],
               by = .(`Sample ID`)]

#### Calibration ####
# Linear model
caldates <- mama$data[(`Sample Type` == "Standard" & Conc == 5),
                      unique(Acquisition)]

mama$markers <- mama$targets[Acquisition %in% caldates & Name %in% markers &
                               !(`Sample ID` == 919 & Polymer == "Polypropylene") &
                               Polymer != "Polystyrene-d5"]

mama$markers[Polymer == "Polyethylene", Label := "PE via 22:2(1,21)"]
mama$markers[Polymer == "Polypropylene", Label := "PP via 2,4Me9:1(1)"]
mama$markers[Polymer == "Polystyrene", Label := "PS via Sty"]

(mama$cal <- mama$markers[, batchcal(`Integrated Area` ~ Conc),
                          by = .(CW, `Cal Type`, Polymer, Label, Name)])

mama$pred <- mama$markers[, batchpred(`Integrated Area` ~ Conc),
                          by = .(CW, `Cal Type`, Polymer, Label, Name)]

ggplot(data = mama$markers, aes(Conc, `Integrated Area` * 10^-6)) +
  geom_ribbon(data = mama$pred, alpha = .25,
              aes(y = fit * 10^-6, ymin = lwr * 10^-6, ymax = upr * 10^-6,
                  fill = `Cal Type`)) +
  geom_line(data = mama$pred, aes(y = fit * 10^-6, color = `Cal Type`)) +
  geom_point(size = pt, aes(color = `Cal Type`, shape = `Cal Type`)) +
  facet_wrap(~ Label, scales = "free") +
  xlab(expression("Polymer concentration"~"["*µg~mL^-1*"]")) +
  ylab(expression("Peak area"~"["%.%~10^6*"]")) +
  scale_shape_discrete(name = "Calibration") +
  scale_color_viridis(discrete = T, name = "Calibration",
                      option = "inferno", begin = .1, end = .9) +
  scale_fill_viridis(discrete = T, name = "Calibration",
                     option = "inferno", begin = .1, end = .9) +
  # theme_publish(base_size = 12) +
  theme_publish(base_family = font)
ggsave('../../../figures/matrix-match.pdf', scale = 1.5, width = pagewidth,
       height = 2.5, unit = 'in', device = cairo_pdf)

mama$pe <- mama$markers[Polymer == "Polyethylene"]
matrix_effect(calibration(`Integrated Area` ~ Conc,
                          data = mama$pe[grepl("TCB", `Cal Type`)]),
              calibration(`Integrated Area` ~ Conc,
                          data = mama$pe[grepl("matrix", `Cal Type`)]))

mama$pp <- mama$markers[Polymer == "Polypropylene"]
matrix_effect(calibration(`Integrated Area` ~ Conc,
                          data = mama$pp[grepl("TCB", `Cal Type`)]),
              calibration(`Integrated Area` ~ Conc,
                          data = mama$pp[grepl("matrix", `Cal Type`)]))

mama$ps <- mama$markers[Polymer == "Polystyrene"]
matrix_effect(calibration(`Integrated Area` ~ Conc,
                          data = mama$ps[grepl("TCB", `Cal Type`)]),
              calibration(`Integrated Area` ~ Conc,
                          data = mama$ps[grepl("matrix", `Cal Type`)]))
