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

nacl <- list()

nacl$reports <- rbind(
  data.table(read_openchrom("../data/py-gc-ms/reports/calibration/")),
  data.table(read_openchrom("../data/py-gc-ms/reports/recovery/")),
  data.table(read_openchrom("../data/py-gc-ms/reports/interferences/")),
  data.table(read_openchrom("../data/py-gc-ms/reports/repeatability/"))
)
names(nacl$reports)[names(nacl$reports) == "Comments"] <- "Extra"

recovtable <- fread("../data/py-gc-ms/recovery.csv")

inftable <- fread("../data/py-gc-ms/interferences.csv")

#### Arrange measurement data ####
nacl$data <- merge(seqtable, nacl$reports, by = "File Name")

nacl$data[grepl("(ZS-)|(Interfer*)", `Sample Name`),
          `Sample Type` := "Sample"]
nacl$data[grep("ug/mL", `Sample Name`, ignore.case = T),
          `Sample Type` := "Standard"]
nacl$data[`Sample Type` == "Standard",
               Conc := as.numeric(gsub(" ug/mL.*", "",
                                       `Sample Name`, ignore.case = T))]
nacl$data[, Acquisition := as.POSIXct(Date)]
nacl$data[, CW := isoweek(Acquisition)]

nacl$all <- merge(nacl$data,
                  allcompounds[, .(Name, Polymer, `Reference Identifier`)],
                  by = "Name")

# Check for contamination and other problems
nacl$all[!(Comment %in% c("", "Repeatability")) & Name == "Styrene-d5",
         .(`Sample ID`, `Sample Name`, Comment)]

# Remove measurements
nacl$targets <- nacl$all[Extra == "quantifiable" &
                           !(`Sample ID` %in% c(479, 487, 496))]
# 479: Calibration outlier (Conc == 200)
# 487 & 496: Displaced quartz tube

nacl$all[`Sample ID` %in% c(479, 487, 496) &
           Name == "Styrene", .(`Sample Name`, Comment)]

#### Repeatability ####
# Check for internal repeatability per day
nacl$targets[grepl("*-d5$", Name),
             .(RSD = sd(`Integrated Area`) / mean(`Integrated Area`)),
             by = .(Acquisition, Name)]

ggplot(nacl$targets[grepl("*-d5$", Name)],
       aes(`Sample ID`, `Integrated Area`)) +
  geom_point() +
  geom_text(aes(label = `Sample ID`), size = 2.5, nudge_x = 2) +
  facet_wrap(. ~ Name, scales = "free_y") + theme_publish()

# Internal standard correction?
nacl$targets[, `IS corrected Area` := `Integrated Area` /
               `Integrated Area`[Name == "Styrene-d5"],
               by = .(`Sample ID`)]

nacl$targets[Comment == "Repeatability",
             .(AreaRSD = sd(`Integrated Area`) / mean(`Integrated Area`),
               IScorRSD = sd(`IS corrected Area`) / mean(`IS corrected Area`)),
             by = .(Name)]

ggplot(nacl$targets[Comment == "Repeatability" & Name %in% markers],
       aes(`Sample ID`, `Integrated Area`)) +
  geom_point() +
  facet_wrap(. ~ Name, scales = "free_y") + theme_publish()

#### Calibration ####
# Linear model
caldates <- nacl$data[(`Sample Type` == "Standard" & Conc == 5),
                      unique(Acquisition)]

# User Conc == 2 for automatic LOD calculations
nacl$targets[Conc == 0, Conc := NA]
nacl$targets[Conc == 2, Conc := 0]

(nacl$cal <- nacl$targets[Acquisition %in% caldates &
                            Polymer != "Polystyrene-d5",
                          batchcal(`Integrated Area` ~ Conc),
                          by = .(CW, Polymer, Name)])

nacl$plot <- nacl$targets[Acquisition %in% caldates &
                            Name %in% markers & Polymer != "Polystyrene-d5"]
nacl$pred <- nacl$plot[, batchpred(`Integrated Area` ~ Conc),
                       by = .(CW, Polymer, Name)]

ggplot(data = nacl$plot, aes(Conc, `Integrated Area` * 10^-6)) +
  geom_ribbon(data = nacl$pred, alpha = .25,
              aes(y = fit * 10^-6, ymin = lwr * 10^-6, ymax = upr * 10^-6,
                  fill = Name)) +
  geom_line(data = nacl$pred, aes(y = fit * 10^-6, color = Name)) +
  geom_point(size = pt, aes(color = Name, shape = Name)) +
  facet_wrap(~ Polymer, scales = "free") +
  xlab(expression("Polymer concentration"~"["*µg~mL^-1*"]")) +
  ylab(expression("Peak area"~"["%.%~10^6*"]")) +
  scale_shape_discrete(name = "Pyrolysates") +
  scale_color_viridis(discrete = T, name = "Pyrolysates",
                      option = "inferno", begin = .1, end = .9) +
  scale_fill_viridis(discrete = T, name = "Pyrolysates",
                     option = "inferno", begin = .1, end = .9) +
  theme_publish(base_family = font, base_size = 12)

# Standard bracketing
ggplot(nacl$targets[Conc == 100 & Name %in% markers],
       aes(Acquisition, `Integrated Area`)) +
  geom_text(aes(label = `Sample ID`), size = 2.5) +
  facet_wrap(. ~ Name, scales = "free_y") + theme_publish()

nacl$quant <- merge(nacl$targets, nacl$cal, by = c("Polymer", "Name", "CW"))

nacl$quant[`Sample Type` == "Standard" & Conc == 100,
           `Bracketing Factor` := `Integrated Area` / (Conc * Slope + Int)]

# Daily factor
.brackets <- nacl$quant[`Sample Type` == "Standard" & Conc == 100,
                        .(`Bracketing Factor` = mean(`Bracketing Factor`)),
                        by = .(Acquisition, Name)]
nacl$quant[`Sample Type` == "Sample",
           `Bracketing Factor` :=
             setDT(.brackets)[nacl$quant[`Sample Type` == "Sample"],
                              `Bracketing Factor`, roll = "nearest",
                              on = c("Name", "Acquisition")]
           ]
rm(.brackets)

nacl$quant[, `Bracketed Area` := `Integrated Area` / `Bracketing Factor`]

#### Calculate concentrations and contents ####
nacl$quant[`Sample Type` == "Sample", Conc := (`Bracketed Area` - Int)/Slope]
nacl$quant[Conc < 0, Conc := 0]
nacl$quant[Conc < LOD, Info := "below LOD"]

nacl$recov <- merge(nacl$quant, recovtable, by = c("Polymer", "Sample Name"))

nacl$recov[, CleanConc := Conc - mean(Conc[Soil == "Blank"], na.rm = T),
           by = .(Polymer, Name)]
nacl$recov[CleanConc < 0, CleanConc := 0]
nacl$recov[, Content := CleanConc * (Extract / `Soil mass`) * Dilution]

nacl$inf <- merge(nacl$quant, inftable, by = c("Polymer", "Sample Name"))

nacl$inf[, CleanConc := Conc - mean(Conc[Soil == "Blank"], na.rm = T),
         by = .(Polymer, Name)]
nacl$inf[CleanConc < 0, CleanConc := 0]
nacl$inf[, Content := CleanConc * (Extract / `Soil mass`) * Dilution]

#### Background contents and interferences ####
(nacl$ref <- nacl$recov[Soil == "Reference" & Name %in% markers,
                       .(`Actual Conc` = `Actual weight` / Extract,
                         CleanConc),
                       by = .(Polymer, Name, `Target weight`)]
 )

nacl$backgr <- nacl$recov[`Target content` == 0 &
                            !(Soil %in% c("Blank", "Reference")),
                          .(Mean = mean(Conc * (Extract / `Soil mass`) * Dilution),
                            SD = sd(Conc * (Extract / `Soil mass`) * Dilution),
                            mLOD = sd(Conc * (Extract / `Soil mass`) * Dilution)
                            * -qt(0.01, 6 - 1) * sqrt(1/.N + 1/6)),
                          by = .(Polymer, Name, Soil)]
nacl$backgr[, .(mLOD = mean(mLOD)), by = .(Polymer, Name)]

nacl$inf[Soil != "Blank" & Name %in% markers,
         .(Mean = mean(Content),
           SD = sd(Content)),
         by = .(Polymer, Name, Soil)]

#### Recovery ####
nacl$recov[, Recovery := Content / (`Actual weight` / `Soil mass`)]

nacl$sum <- nacl$recov[`Target content` != 0,
                       .(RecMean = mean(Recovery),
                         RecSD = sd(Recovery)),
                       by = .(Polymer, Name, `Target content`, Soil)]

nacl$sum[, Rec := signifig(RecMean * 100, as.numeric(RecSD * 100)),
         by = Soil]
nacl$sum[Name %in% markers]

instr <-
  nacl$cal[, .(Polymer, Name, Radj, LOD)] %>%
  merge(nacl$targets[Comment == "Repeatability",
                     .(RSD = sd(`Integrated Area`) / mean(`Integrated Area`)),
                     by = .(Name)],
        by = "Name")

instr[, .(Polymer, Pyrolysate = factor(Name, labels = pyr[Name]),
          Radj = round(Radj, 4),
          LOD = round(LOD, 1) * 2, RSD = round(100 * RSD, 1))]

meth <-
  nacl$backgr[, .(mLOD = round(mLOD, 1)),
            by = .(Polymer, Name, Soil)] %>%
  merge(nacl$inf[Soil != "Blank",
                 .(Selectivity = signifig(mean(Content), sd(Content))),
                 by = .(Polymer, Name, Soil)],
        by = c("Polymer", "Name", "Soil"), all.x = T) %>%
  merge(
    dcast(nacl$sum, Polymer + Name + Soil ~ `Target content`,
          value.var = "Rec"),
    by = c("Polymer", "Name", "Soil")
  )

meth[order(Soil), .(Polymer, Pyrolysate = factor(Name, labels = pyr[unique(Name)]),
                    mLOD, Selectivity, `2`, `20`)]
