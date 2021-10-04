setwd("~/Documents/PhD/Thesis/")

library(ggplot2)
library(cowplot)
library(envalysis)
library(data.table)

trace <- fread("proposal/pyrogram.csv")
trace[, Rel := Intensity / max(Intensity)]


compl <- ggplot(trace, aes(Time, Rel * 100)) +
  geom_line(size = .2) +
  xlab('Time [min]') + ylab('Relative intensity [%]') +
  theme_publish(base_family = "Source Sans Pro")
pip <- ggplot(trace[Time %between% c(8.8,9.3)], aes(Time, Rel * 100)) +
  geom_line(size = .2) +
  ylab("") +
  scale_x_continuous(breaks = c(9,9.2), name = "") +
  theme_publish(base_family = "Source Sans Pro")

ggdraw() +
  draw_plot(compl) +
  draw_plot(pip, 0.20, 0.45, 0.2, 0.5)
ggsave('proposal/pyrogram.pdf', scale = 1.5, width = 4.21342,
       height = 2, unit = 'in', device = cairo_pdf)

