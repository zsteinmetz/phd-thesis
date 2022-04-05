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

fisa <- list()

fisa$reports <- rbind(
  data.table(read_openchrom("../data/py-gc-ms/reports/calibration/")),
  data.table(read_openchrom("../data/py-gc-ms/reports/field-samples/"))
)
names(fisa$reports)[names(fisa$reports) == "Comments"] <- "Extra"

samtable <- fread("../data/py-gc-ms/field-samples.csv")

#### Arrange measurement data ####
fisa$data <- merge(seqtable, fisa$reports, by = "File Name")

fisa$data[is.na(`Integrated Area`), `Integrated Area` := 0]
fisa$data[grepl("ZS-", `Sample Name`), `Sample Type` := "Sample"]
fisa$data[grep("ug/mL", `Sample Name`, ignore.case = T),
          `Sample Type` := "Standard"]
fisa$data[`Sample Type` == "Standard",
               Conc := as.numeric(gsub(" ug/mL.*", "",
                                       `Sample Name`, ignore.case = T))]
fisa$data[, Acquisition := as.POSIXct(Date)]
fisa$data[, CW := isoweek(Acquisition)]

fisa$all <- merge(fisa$data,
                  allcompounds[, .(Name, Polymer, `Reference Identifier`)],
                  by = "Name")

# Check for contamination and other problems
fisa$all[!(Comment %in% c("") | is.na(Comment)) & Name == "Styrene-d5",
         .(`Sample ID`, `Sample Name`, Comment)]

# Remove measurements
fisa$targets <- fisa$all[Extra == "quantifiable" &
                           !(`Sample ID` %in% c(
                             # Calibration outliers
                             479, 542, 605, 668, 740, 810, 816,
                             # Bracketing outliers
                             562,
                             # Empty tube
                             597,
                             # Conc too high
                             645,
                             # Dilution too high / repetition
                             891, 896, 898,
                             # Contamination
                             618, 771,
                             # Incomplete pyrolysis (IS low in intensity)
                             568, 569, 594, 599, 616, 636, 637, 638, 647,
                             793, 795, 842,
                             # Vial open
                             798))]

#### Repeatability ####
# Check for internal repeatability per day
fisa$targets[grepl("*-d5$", Name) &
               !(`Sample ID` %in% c(819, 861, 870, 871, 874, 891)),
             .(RSD = sd(`Integrated Area`) / mean(`Integrated Area`)),
             by = .(Acquisition, Name)]
# IDs 682, 747, 819, 861, 870, 871, 874, 891 are (far) below LOD anyway
# off IS thus negligible

ggplot(fisa$targets[Name == "Styrene-d5"],
       aes(`Sample ID`, `Integrated Area`)) +
  # geom_point() +
  geom_text(aes(label = `Sample ID`), size = 2.5) +
  theme_publish()

# Internal standard correction?
fisa$targets[, `IS corrected Area` := `Integrated Area` /
               `Integrated Area`[Name == "Styrene-d5"],
               by = .(`Sample ID`)]

#### Calibration ####
# Linear model
caldates <- fisa$data[(`Sample Type` == "Standard" & Conc %in% c(5, 10, 20, 50)),
                      unique(Acquisition)]

# User Conc == 2 for automatic LOD calculations
fisa$targets[Conc == 0, Conc := NA]
fisa$targets[Conc == 2, Conc := 0]

fisa$cal <- fisa$targets[Acquisition %in% caldates &
                            Polymer != "Polystyrene-d5",
                          batchcal(`Integrated Area` ~ Conc),
                          by = .(CW, Polymer, Name)]
fisa$cal[Name %in% markers]

fisa$plot <- fisa$targets[Acquisition %in% caldates &
                            Name %in% markers & Polymer != "Polystyrene-d5"]
fisa$pred <- fisa$plot[, batchpred(`Integrated Area` ~ Conc, check_assumptions = F),
                       by = .(CW, Polymer, Name)]

ggplot(data = fisa$plot, aes(Conc, `Integrated Area` * 10^-6)) +
  geom_ribbon(data = fisa$pred, alpha = .25,
              aes(y = fit * 10^-6, ymin = lwr * 10^-6, ymax = upr * 10^-6,
                  fill = Name)) +
  geom_line(data = fisa$pred, aes(y = fit * 10^-6, color = Name)) +
  geom_point(size = pt, aes(color = Name, shape = Name)) +
  facet_grid(Polymer ~ CW, scales = "free_y") +
  xlab(expression("Polymer concentration"~"["*µg~mL^-1*"]")) +
  ylab(expression("Peak area"~"["%.%~10^6*"]")) +
  scale_shape_discrete(name = "Pyrolysates") +
  scale_color_viridis(discrete = T, name = "Pyrolysates",
                      option = "inferno", begin = .1, end = .9) +
  scale_fill_viridis(discrete = T, name = "Pyrolysates",
                     option = "inferno", begin = .1, end = .9) +
  theme_publish(base_family = font, base_size = 12)

# Standard bracketing
ggplot(fisa$targets[Conc == 100 & Name %in% markers],
       aes(Acquisition, `Integrated Area`)) +
  geom_text(aes(label = `Sample ID`), size = 2.5) +
  facet_wrap(. ~ Name, scales = "free_y") + theme_publish()

fisa$quant <- merge(fisa$targets,
                    fisa$cal[, .(CW, Polymer, Name, Int, Slope, Radj)],
                    by = c("Polymer", "Name", "CW"))

fisa$quant[`Sample Type` == "Standard" & Conc == 100,
           `Bracketing Factor` := `Integrated Area` / (Conc * Slope + Int)]

# Daily factor
.brackets <- fisa$quant[`Sample Type` == "Standard" & Conc == 100,
                        .(`Bracketing Factor` = mean(`Bracketing Factor`)),
                        by = .(Acquisition, Name)]
fisa$quant[`Sample Type` == "Sample",
           `Bracketing Factor` :=
             setDT(.brackets)[fisa$quant[`Sample Type` == "Sample"],
                              `Bracketing Factor`, roll = "nearest",
                              on = c("Name", "Acquisition")]
           ]
rm(.brackets)

fisa$quant[, `Bracketed Area` := `Integrated Area` / `Bracketing Factor`]

#### Calculate concentrations and contents ####
fisa$quant[`Sample Type` == "Sample", Conc := (`Bracketed Area` - Int)/Slope]
fisa$quant[Conc < 0, Conc := 0]
fisa$quant[Conc < 5, Info := "conc too small"]
fisa$quant[Conc > 200, Info := "conc too high"]

fisa$quant[Info == "conc too high" & Name %in% markers,
           .(`Sample Name`, Polymer, Conc)]

fisa$final <- merge(fisa$quant[Name %in% markers], samtable, by = "Sample Name")

fisa$final[, CleanConc := Conc - mean(Conc[Location == "Blank"], na.rm = T),
           by = .(Polymer, Name, `Separation CW`)]
fisa$final[CleanConc < 0, CleanConc := 0]
fisa$final[Location != "Blank",
           Content := CleanConc * (Extract / `Soil mass`) * Dilution * `Dil Factor`]

# Check for missing/duplicate samples
fisa$final[Site != "B",
           .(.N,
             Measured = paste(`Sample Name`, collapse = ", ")),
           by = .(Name, Site, Transect, Row)] %>%
  .[N != 5]

fisa$sum <- fisa$final[Location != "Blank"]
fisa$sum[, Site := paste("Site", Site)]

fisa$labels <- fisa$sum[Polymer == "Polyethylene",
                        .(Label = unique(Cover)), by = .(Polymer, Site)]

fisa$mean <- fisa$sum[, .(Content = mean(Content),
                          SD = sd(Content)),
                      by = .(Polymer, Site, Transect, Row)]

ggplot(data = fisa$sum, aes(Transect, Content)) +
  geom_col(data = fisa$mean, aes(fill = Row, color = Row),
           position = position_dodge(1)) +
  geom_point(aes(color = Row, shape = Row), fill = "white",
             position = position_dodge(1), size = pt) +
  facet_grid(Polymer ~ Site, scales = "free", labeller = polymer_labeller) +
  scale_y_continuous(name = expression("Polymer content"~"[mg kg"^-1*"]"),
                     breaks = c(1, 2.5, 5, 10, 25, 50),
                     trans = scales::pseudo_log_trans(base = 10)) +
  scale_x_discrete(labels = c("FC", "FE", "FM", "P")) +
  scale_shape_manual(breaks = c("Planting", "Track", "None"),
                     values = c(22, 23, 21)) +
  scale_color_manual(breaks = c("Planting", "Track", "None"),
                     values = c(viridis(3)[c(2,1)], "gray20")) +
  scale_fill_manual(breaks = c("Planting", "Track", "None"),
                    values = c(viridis(3)[c(2,1)], "gray20")) +
  theme_publish(base_family = font) +
  theme(strip.text.y = element_text(angle = 0))
ggsave('../../../figures/py-screening.pdf', scale = 1.5, width = pagewidth,
       height = 5, unit = 'in', device = cairo_pdf)
