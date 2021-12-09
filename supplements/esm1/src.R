#### Load packages and functions ####

setwd("~/Documents/PhD/Thesis/supplements/")
source("R/init.R")

#### Examplary TGA/DTG curve ####

example <- list()
example$data <- fread("esm1/data/tga-dtg-sample.csv")
example$data[, Sample := factor(Sample, c("LUFA soil", "PET in soil", "PET dust"))]

lty <- c("solid", "dashed", "solid")

example$tga <- 
  ggplot(example$data, aes(Temp, TGA)) +
  geom_line(aes(color = Sample)) +
  ylab("TGA mass [%]") +
  scale_color_viridis(option = "inferno", begin = .1, end = .9, direction = -1,
                      discrete = T) +
  theme_publish(base_family = font) + theme_red()
example$dtg <- 
  ggplot(example$data[Temp > 40], aes(Temp, DTG)) +
  geom_line(aes(color = Sample)) +
  ylab(expression("DTG [% K"^-1*"]")) + 
  scale_color_viridis(option = "inferno", begin = .1, end = .9, direction = -1,
                      discrete = T) +
  theme_publish(base_family = font) + theme_red()
example$mass <- 
  ggplot(example$data, aes(Temp, mz105*10^9)) +
  geom_line(aes(color = Sample)) +
  xlab("Temperature [°C]") + ylab(bquote(105~italic("m/z")~"["%.% 10^-9*"]")) +
  scale_color_viridis(option = "inferno", begin = .1, end = .9, direction = -1,
                      discrete = T,
                      labels = c("Blank soil", "PET debris in soil",
                                 "Pure PET")) +
  theme_publish(base_family = font) + theme(legend.direction = "vertical")
example$tga / example$dtg / example$mass + plot_annotation(theme = theme(plot.margin = margin()))
ggsave('../figures/tga-dtg-sample.pdf', scale = 1.5, width = marginwidth,
       height = 3.5, unit = 'in', device = cairo_pdf)

#### TGA/MS curves ####

tga <- list()
tga$data <- fread("esm1/data/tga-ms.csv")
  
tga$data[, Label := paste0(toupper(substring(IS, 1, 1)), substring(IS, 2), " internal standard")]

tga$tga <-
  ggplot(tga$data, aes(Temp, TGA)) +
  geom_line(aes(color = PET, group = Sample)) +
  ylab("TGA mass [%]") +
  facet_wrap(~ Label, scales = "free") +
  scale_color_viridis(option = "inferno", begin = .1, end = .9, direction = -1) +
  theme_publish(base_family = font) + theme_red()
tga$mz105 <- 
  ggplot(tga$data, aes(Temp, mz105*10^9)) +
  geom_line(aes(color = PET, group = Sample)) +
  xlab("Temperature [°C]") + ylab(bquote(105~italic("m/z")~"["%.% 10^-9*"]")) +
  facet_wrap(~ Label, scales = "free") +
  scale_color_viridis(option = "inferno", begin = .1, end = .9, direction = -1) +
  theme_publish(base_family = font) + theme_red() + theme(strip.text.x = element_text(size=0))
tga$mz154 <- 
  ggplot(tga$data, aes(Temp, mz154*10^9)) +
  geom_line(aes(color = PET, group = Sample)) +
  facet_wrap(~ Label, scales = "free") +
  xlab("Temperature [°C]") + ylab(bquote(154~italic("m/z")~"["%.% 10^-9*"]")) +
  scale_color_viridis(option = "inferno", begin = .1, end = .9, direction = -1,
                      name = bquote(bold("PET content [g kg"^-1*"]"))) +
  theme_publish(base_family = font) + theme(strip.text.x = element_text(size=0)) +
  guides(color = guide_colorbar(title.vjust = 1),
         shape = guide_legend(title.vjust = 1))
tga$tga / tga$mz105 / tga$mz154 + plot_annotation(theme = theme(plot.margin = margin()))
ggsave('../figures/tga-ms.pdf', scale = 1.5, width = pagewidth,
       height = 4, unit = 'in', device = cairo_pdf)

#### Calibration ####

cal <- list()
cal$data <- fread("esm1/data/tga-calibration.csv")

cal$data[, Label := paste0(toupper(substring(IS, 1, 1)), substring(IS, 2), " internal standard")]

(cal$cal105 <- cal$data[, batchcal(mz105norm ~ PET), by = .(Label)])
(cal$cal154 <- cal$data[, batchcal(mz154norm ~ PET), by = .(Label)])

cal$pred105 <- cal$data[, batchpred(mz105norm ~ PET), by = .(Label)]
cal$pred154 <- cal$data[, batchpred(mz154norm ~ PET), by = .(Label)]

cal$plot105 <-
  ggplot(data = cal$data, aes(PET, mz105norm*10^9)) +
  geom_ribbon(data = cal$pred105, alpha = .25, fill = inferno(1, begin = .1, end = .9),
              aes(y = fit*10^9, ymin = lwr*10^9, ymax = upr*10^9)) +
  geom_line(data = cal$pred105, aes(y = fit*10^9), color = inferno(1, begin = .1, end = .9)) +
  geom_point(size = pt, color = inferno(1, begin = .1, end = .9)) +
  facet_wrap(~ Label, scales = "free") +
  xlab(bquote("PET content"~"[g kg"^-1*"]")) +
  ylab(bquote(105~italic("m/z")~"["%.% 10^-9*"]")) +
  scale_y_continuous(label=scientific_10) +
  theme_publish(base_family = font) + theme_red()
cal$plot154 <-
  ggplot(data = cal$data, aes(PET, mz154norm*10^9)) +
  geom_ribbon(data = cal$pred154, alpha = .25, fill = inferno(1, begin = .1, end = .9),
              aes(y = fit*10^9, ymin = lwr*10^9, ymax = upr*10^9)) +
  geom_line(data = cal$pred154, aes(y = fit*10^9), color = inferno(1, begin = .1, end = .9)) +
  geom_point(size = pt, color = inferno(1, begin = .1, end = .9)) +
  facet_wrap(~ Label, scales = "free") +
  xlab(bquote("PET content"~"[g kg"^-1*"]")) +
  ylab(bquote(154~italic("m/z")~"["%.% 10^-9*"]")) +
  scale_y_continuous(label=scientific_10) +
  theme_publish(base_family = font) + theme(strip.text.x = element_text(size=0))
cal$plot105 / cal$plot154 + plot_annotation(theme = theme(plot.margin = margin()))
ggsave('../figures/tga-calibration.pdf', scale = 1.5, width = pagewidth,
       height = 3, unit = 'in', device = cairo_pdf)

#### Capillary checks ####

cap <- list()
cap$data <- fread("esm1/data/tga-capillary-checks.csv")

cap$data[, Label := paste0(toupper(substring(IS, 1, 1)), substring(IS, 2), " internal standard")]
cap$data[PET > 30, Case := "40–50"]
cap$data[PET < 10, Case := "2"]

lty <- c("dotted", "dashed", "solid")

ggplot(cap$data, aes(Temp, mz105*10^9, group = interaction(PET, Case, Type))) +
  geom_line(aes(color = PET, linetype = Type)) +
  facet_wrap(~ Label, scales = "free") +
  xlab("Temperature [°C]") +   ylab(bquote(105~italic("m/z")~"["%.% 10^-9*"]")) +
  scale_color_viridis(option = "inferno", begin = .1, end = .9, direction = -1,
                      name = bquote(bold("PET content [g kg"^-1*"]"))) +
  scale_linetype_manual(values = lty) +
  theme_publish(base_family = font) +
  guides(color = guide_colorbar(title.vjust = 1),
         shape = guide_legend(title.vjust = 1))
ggsave('../figures/tga-capillary-checks.pdf', scale = 1.5, width = pagewidth,
       height = 2.2, unit = 'in', device = cairo_pdf)

#### Capillary pressure ####

cap$pressure <- fread("esm1/data/tga-capillary-pressure.csv")
cap$pressure[, Date := as.Date(Date)]

ggplot(cap$pressure, aes(x = Date, y = Pressure * 10^6)) +
  geom_line(color = cividis(1)) +
  geom_point(color = cividis(1)) +
  geom_hline(yintercept = 8, linetype = "dashed") +
  geom_text(x = -Inf, y = 7, label = "'capillary blocked' < 8 %.% 10^-9~'bar'",
            hjust = -0.5, parse = T, size = 3, family = font) + 
  geom_hline(yintercept = 20, linetype = "dashed") +
  geom_text(x = -Inf, y = 21, label = "'capillary broken' > 20 %.% 10^-9~'bar'",
            hjust = -0.5, parse = T, size = 3, family = font) + 
  ylab(expression("Pressure [10"^-9~"bar]")) +
  scale_y_continuous(limits = c(7, 21)) +
  theme_publish(base_family = font)
ggsave('../figures/tga-capillary-pressure.pdf', scale = 1.5, width = textwidth,
       height = 2, unit = 'in', device = cairo_pdf)

#### Polymer DSC ####

dsc <- list()
dsc$data <- fread("esm1/data/dsc.csv")

ggplot(dsc$data, aes(Temp, AvgHF)) +
  geom_line(aes(color = Sample)) +
  xlab("Temperature [°C]") + ylab(expression("Average heat flow [W g"^-1*"]")) +
  scale_color_viridis(option = "cividis", direction = -1,
                      discrete = T, name = "PET sample",
                      labels = c("Bottle (EDEKA)",
                                 "Dust (PETKA CZ)",
                                 "Reference (PlasticsEurope)")) +
  theme_publish(base_family = font)
ggsave('../figures/dsc.pdf', scale = 1.5, width = textwidth,
       height = 2, unit = 'in', device = cairo_pdf)

#### Cysteine comparison ####

cys <- list()
cys$comp <- fread("esm1/data/cysteine-comparison.csv")
cys$comp[, Cys := as.factor(Cys)]

cys$comptga <-
  ggplot(cys$comp, aes(Temp, TGA)) +
  geom_line(aes(color = Cys, group = Cys)) +
  ylab("TGA mass [%]") +
  scale_color_viridis(option = "cividis", direction = -1, discrete = T) +
  theme_publish(base_family = font) + theme_red()
cys$compmz33 <- 
  ggplot(cys$comp, aes(Temp, mz33*10^9)) +
  geom_line(aes(color = Cys, group = Cys)) +
  xlab("Temperature [°C]") + ylab(bquote(33~italic("m/z")~"["%.% 10^-9*"]")) +
  scale_color_viridis(option = "cividis", direction = -1, discrete = T,
                      name = "Cysteine content [%]") +
  theme_publish(base_family = font) + theme(strip.text.x = element_text(size=0))
cys$comptga / cys$compmz33 +
  plot_annotation(theme = theme(plot.margin = margin()))
ggsave('../figures/cysteine-comparison.pdf', scale = 1.5, width = textwidth,
       height = 3, unit = 'in', device = cairo_pdf)

#### Cysteine tests ####

cys$lin <- fread("esm1/data/cysteine-linearity.csv")
cys$deg <- fread("esm1/data/cysteine-degradation.csv")

cys$cal <- calibration(mz33 ~ Cys, data = cys$lin)
cys$conc <- cys$cal$model$model[,2]
cys$new <- data.frame(conc = seq(min(cys$conc), max(cys$conc),
                                 length.out = 100 * length(cys$conc)))
names(cys$new) <- all.vars(cys$cal$model$formula)[2]
cys$pred <- data.frame(cys$new, predict(cys$cal$model, cys$new, interval = "conf"))

cys$cal
summary(cys$cal$model)

cys$calcurve <-
  ggplot(data = cys$lin, aes(Cys, mz33*10^6)) +
  geom_ribbon(data = cys$pred, alpha = .25, fill = viridis(1, begin = .1, end = .9),
              aes(y = fit*10^6, ymin = lwr*10^6, ymax = upr*10^6)) +
  geom_line(data = cys$pred, aes(y = fit*10^6), color = viridis(1, begin = .1, end = .9)) +
  geom_point(color = cividis(1)) +
  xlab("Cysteine content [%]") +
  ylab(bquote(33~italic("m/z")~"["%.% 10^-6*"]")) +
  theme_publish(base_family = font) +
  ggtitle("a")

cys$tga <-
  ggplot(cys$deg, aes(Temp, TGA)) +
  geom_line(aes(color = Cys, group = Sample)) +
  ylab("TGA mass [%]") +
  scale_color_viridis(option = "cividis", direction = -1) +
  theme_publish(base_family = font) + theme_red() +
  ggtitle("b")
cys$mz33 <- 
  ggplot(cys$deg, aes(Temp, mz33*10^9)) +
  geom_line(aes(color = Cys, group = Sample)) +
  xlab("Temperature [°C]") + ylab(bquote(33~italic("m/z")~"["%.% 10^-9*"]")) +
  scale_color_viridis(option = "cividis", direction = -1,
                      name = "Cysteine content [%]") +
  theme_publish(base_family = font) + theme(strip.text.x = element_text(size=0)) +
  guides(color = guide_colorbar(title.vjust = 1),
         shape = guide_legend(title.vjust = 1))
cys$calcurve / cys$tga / cys$mz33 +
  plot_annotation(theme = theme(plot.margin = margin()))
ggsave('../figures/cysteine-tests.pdf', scale = 1.5, width = textwidth,
       height = 4, unit = 'in', device = cairo_pdf)

#### Control measurements ####
contr <- list()
contr$data <- fread("esm1/data/control-measurements.csv")
contr$melt <- melt(contr$data, id = c("PET", "IS"),
                   measure = c("mz18", "mz28", "mz44", "SigRatio"),
                   variable.name = "Mass", value.name = "Signal")

contr$melt[, Label := paste0(toupper(substring(IS, 1, 1)), substring(IS, 2), " internal standard")]

contr$sigratio <- 
  ggplot(contr$melt[Mass == "SigRatio",], aes(PET, Signal)) +
  geom_line() +
  geom_point() +
  ylab("Signal ratio") +
  facet_wrap(~ Label) +
  theme_publish(base_family = font) + theme_red()
contr$masses <- 
  ggplot(contr$melt[Mass != "SigRatio",], aes(PET, Signal * 10^9)) +
  geom_line(aes(color = Mass), linetype = "dashed") +
  geom_point(aes(color = Mass)) +
  ylab(bquote(105~italic("m/z")~"["%.% 10^-9*"]")) +
  scale_color_viridis(option = "cividis", direction = -1,
                      discrete = T, name = bquote(bolditalic("m/z")),
                      labels = c("18", "28", "44")) +
  xlab(bquote("PET content"~"[g kg"^-1*"]")) +
  ylab(bquote(italic("m/z")~"["%.% 10^-9*"]")) +
  theme_publish(base_family = font) + theme(strip.text.x = element_text(size=0)) +
  facet_wrap(~ Label)
contr$sigratio / contr$masses +
  plot_annotation(theme = theme(plot.margin = margin()))
ggsave('../figures/control-measurements.pdf', scale = 1.5, width = pagewidth,
       height = 3.5, unit = 'in', device = cairo_pdf)
