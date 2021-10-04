setwd("~/Documents/PhD/Thesis/figures")

library(ggplot2)
library(envalysis)
library(data.table)

breakthrough <- rbind(
  data.table(
    Compound = "Tracer",
    PoreVol = seq(0, 2.5, length=1000),
    Conc = dnorm(seq(0, 2.5, length=1000), mean=1, sd=.12)/1.1
  ),
  data.table(
    Compound = "Plastic particles",
    PoreVol = seq(0, 2.5, length=1000),
    Conc = dnorm(seq(0, 2.5, length=1000), mean=.6, sd=.07)/6
  )
)

ggplot(breakthrough, aes(PoreVol, Conc)) +
  geom_line(aes(linetype = Compound)) +
  labs(x = "Pore volumes", y = "Normalized concentration") +
  scale_linetype_manual(values = c("solid", "dashed")) +
  theme_publish(base_family = "Source Sans Pro") +
  theme(legend.position = "none")
ggsave('breakthrough.pdf', scale = 1.5, width = 1.94525,
       height = 1.7, unit = 'in', device = cairo_pdf)
