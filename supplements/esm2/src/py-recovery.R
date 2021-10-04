#### Load packages and functions ####

setwd("~/Documents/PhD/Thesis/supplements/")
source("R/init.R")

## Read data
seqtable <- fread("esm2/data/seqtable.csv")

allcompounds <- rbind(
  cbind(read_targets("Targets/Refs.txt", header = F)),
  cbind(read_targets("Targets/PE.txt", header = F)),
  cbind(read_targets("Targets/PP.txt", header = F)),
  cbind(read_targets("Targets/PS.txt", header = F))
)
names(allcompounds)[names(allcompounds) == "Contributor"] <- "Polymer"

rec <- list()
rec$reports <- data.table(read_openchrom("esm2/data/reports/recovery/"))
rec$samples <- fread("esm2/data/samples.csv")

rec$data <- merge(seqtable, rec$reports, by = "File Name")

# Naming and sorting
rec$data[grep("ZS-", `Sample Name`), `Sample Type` := "Sample"]
rec$data[grep("ug/mL", `Sample Name`), `Sample Type`  := "Standard"]
rec$data[`Sample Type` == "Standard",
         Conc := as.numeric(gsub(" ug/mL.*", "", `Sample Name`))]
rec$data[, Acquisition := as.POSIXct(Date)]
rec$data[, CW := week(Acquisition)]

rec$targets <- merge(rec$data[grep("Recovery", Comment) & Comments == "quantifiable"],
                     allcompounds, by = "Name")
# Remove standards with impurities
rec$clean <- rec$targets[!(Name == "2,4-Dimethyl-1-heptene" & `Sample ID` %in% c(49, 52, 150))]

## Calibration
caldates <- rec$data[`Sample Type` == "Standard" & Conc == 0, unique(Acquisition)]
(rec$cal <- rec$clean[`Sample Type` == "Standard" & Acquisition %in% caldates, 
                      batchcal(`Integrated Area` ~ Conc), by = .(CW, Polymer, Name)])

rec$pred <- rec$clean[`Sample Type` == "Standard" & Acquisition %in% caldates, 
                      batchpred(`Integrated Area` ~ Conc), by = .(CW, Polymer, Name)]

ggplot(data = rec$clean[Comments == "quantifiable"], aes(Conc, `Integrated Area`*10^-6)) +
  geom_ribbon(data = rec$pred, alpha = .2,
              aes(y = fit*10^-6, ymin = lwr*10^-6, ymax = upr*10^-6, fill = Name)) +
  geom_line(data = rec$pred, aes(y = fit*10^-6, color = Name)) +
  geom_point(size = pt, aes(color = Name)) +
  facet_wrap(~ interaction(CW, Polymer), scales = "free") +
  xlab(expression("Polymer concentration"~"["*µg~mL^-1*"]")) +
  ylab(expression("Peak area"~"["%.%~10^6*"]")) +
  theme_publish() + theme(legend.box = "vertical")

## Standard bracketing
rec$quant <- merge(rec$clean, rec$cal, by = c("Polymer", "Name", "CW"))

rec$quant[Conc == 100, `Bracketing Factor` := `Integrated Area` / (100 * Slope + Int)]
rec$quant[, `Corrected Area` := `Integrated Area` /
            mean(`Bracketing Factor`, na.rm = T), by = .(CW, Name, Acquisition)]

rec$quant[`Sample Type` == "Sample", Conc := (`Corrected Area` - Int)/Slope]
rec$quant[`Sample Type` == "Sample" & (Conc < 0 | is.na(Conc)), Conc := 0]

rec$final <- merge(rec$quant, rec$samples, by = c("Polymer", "Sample Name"))

# Blank correction
rec$final[, CleanConc := Conc - mean(Conc[Soil == "Blank"], na.rm = T),
          by = .(Polymer, Name, `Extraction procedure`)]
rec$final[is.na(CleanConc), CleanConc := Conc]
rec$final[CleanConc < 0, CleanConc := 0]

# Calculate plastic contents
rec$final[, Content := CleanConc * (Extract / `Soil mass`) * Dilution]
rec$blanks <- rec$final[Soil %in% c("Blank", "Interference blank") &
                        Name %in% c("1,16-Heptadecadiene", "2,4-Dimethyl-1-heptene",
                                    "α-Methylstyrene"),
                        .(Blank = mean(CleanConc * (Extract / 4) * Dilution)),
                        by = .(Name, Soil, `Extraction procedure`)]

rec$final[, CleanContent := Content - mean(Content[`Target content` == 0], na.rm = T),
           by = .(Soil, Polymer, `Extraction procedure`)]
rec$final[is.na(CleanContent), CleanContent := Content]
rec$final[CleanContent < 0, CleanContent := 0]

rec$final <- rec$final[Soil != "Blank"]

## Recovery
rec$final[, Recovery := CleanContent / `Actual content`]

rec$sum <- rec$final[Name %in% c("1,16-Heptadecadiene", "2,4-Dimethyl-1-heptene",
                                 "α-Methylstyrene"),
                     .(ContentMean = mean(Content),
                       ContentSD = sd(Content),
                       RecMean = mean(Recovery),
                       RecSD = sd(Recovery),
                       N = .N),
                     by = .(Polymer, Soil, `Target content`, `Extraction procedure`, `Interference target`)]
rec$sum[`Interference target` == 0.2 & Soil == "LUFA 2.2", Soil := "LUFA 2.2*"]

rec$sum[, `Extraction procedure` := factor(`Extraction procedure`,
                                           levels = unique(`Extraction procedure`))]

pd <- position_dodge(.4)

rec$sum[Soil == "LUFA 2.2*", Soil := "LUFA 2.2 with non-\ntarget polymers*"]

## Plotting
ggplot(rec$sum[`Target content` != 0 & RecMean > .05],
       aes(as.factor(`Target content`), RecMean*100)) +
  geom_rect(xmin = 0, ymin = 70, xmax = Inf, ymax = 130, fill = "snow2", alpha = .25) +
  geom_hline(yintercept = 100, size = sz) +
  geom_errorbar(aes(ymin = (RecMean - RecSD)*100, ymax = (RecMean + RecSD)*100, color = `Extraction procedure`,
                    group = interaction(`Extraction procedure`, Soil)), width = .25, position = pd) +
  geom_point(size = pt, aes(shape = Soil, color = `Extraction procedure`), position = pd) +
  facet_grid(. ~ Polymer, labeller = polymer_labeller) +
  xlab(bquote("Nominal polymer content"~"["*"\U003BCg"~g^-1*"]")) +
  scale_y_continuous(name = "Matrix-corrected recovery [%]", breaks = c(0, 50, 100, 150, 200, 300)) +
  scale_color_viridis(discrete = T, option = "inferno", begin = .1, end = .9,
                      labels = c("TCB only", "Methanol\ncleanup",
                                 bquote(KAl(SO[4])[2]~"flocculation"),
                                 "Fenton digestion")) +
  scale_shape_manual(values = c(19, 1, 17, 15)) +
  theme_publish(base_family = font) + theme(legend.box = "vertical")
ggsave('../figures/py-recovery.pdf', scale = 1.5, width = pagewidth,
       height = 3, unit = 'in', device = cairo_pdf)
