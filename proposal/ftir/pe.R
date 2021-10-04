setwd("~/Documents/PhD/Thesis/proposal/ftir")

source("read_ftir.R")
require(data.table)
require(ggplot2)
require(cowplot)
require(envalysis)

pe <- data.table(read_ftir())

ggplot(pe, aes(Wavenumber, Absorbance)) +
  geom_line(aes(linetype = Measurement)) +
  labs(x = expression("Wavenumber [cm"^-1*"]")) +
  scale_linetype_discrete(name = "") +
  theme_publish(base_family = "Source Sans Pro")
# ggsave('../aging.pdf', scale = 1.5, width = 4.21342,
#        height = 2.2, unit = 'in', device = cairo_pdf)

write.csv(pe[Measurement == "Virgin PE",
             .(Wavelength = Wavenumber, Absorbance = Absorbance)],
          file = "Virgin PE.csv", row.names = F)

# res <- data.table(Waveno = pe[Measurement == "Aged PE", Waveno], 
#                   Signal = pe[Measurement == "Aged PE", Signal] -
#                     pe[Measurement == "Pure PE", Signal])
# 
# ggplot(res, aes(Waveno, Signal)) +
#   geom_line() +
#   labs(x = expression("Wavenumber [cm"^-1*"]"), y = "Absorbance") +
#   theme_publish()

                  