#### Load packages and global vars ####
source("func/init.R")

specs <- rbind(
  data.table(Name = "bold('Oleonitrile')", spec = "NIST08",
             fread("../data/py-gc-ms/nist/oleonitrile_nist.csv")),
  data.table(Name = "bold('Oleonitrile')", spec = "Measured",
             fread("../data/py-gc-ms/nist/oleonitrile_td.csv")),
  data.table(Name = "bold('9-Octadecenamide')", spec = "NIST08",
             fread("../data/py-gc-ms/nist/9-octadecenamide_nist.csv")),
  data.table(Name = "bold('9-Octadecenamide')", spec = "Measured",
             fread("../data/py-gc-ms/nist/9-octadecenamide_td.csv")),
  data.table(Name = "bold('Di-'*bolditalic('tert')*'-butylphenol')", spec = "NIST08",
             fread("../data/py-gc-ms/nist/tert-butylphenol_nist.csv")),
  data.table(Name = "bold('Di-'*bolditalic('tert')*'-butylphenol')", spec = "Measured",
             fread("../data/py-gc-ms/nist/tert-butylphenol_td.csv")),
  data.table(Name = "bold('Propyl dodecanoate')", spec = "NIST08",
             fread("../data/py-gc-ms/nist/propyl dodecanoate_nist.csv")),
  data.table(Name = "bold('Propyl dodecanoate')", spec = "Measured",
             fread("../data/py-gc-ms/nist/9-octadecenamide_td.csv"))
)

specs <- specs[`Ion [m/z]` %between% c(50,280)]
specs[spec == "NIST08", `Intensity [%]` := `Intensity [%]`/10]
specs[spec == "Measured", `Ion [m/z]` := round(`Ion [m/z]`)]

top <- specs[, .(`Intensity [%]` = max(`Intensity [%]`)),
              by = .(Name, `Ion [m/z]`)]
top <- top[order(top$`Intensity [%]`, decreasing = TRUE),] 

ggplot(specs, aes(x = `Ion [m/z]`, ymin = 0, ymax = `Intensity [%]`)) +
  geom_linerange(aes(color = spec), position = position_dodge(.5)) +
  geom_text(data = top[ , head(.SD, 5), by = Name],
            nudge_y = 8, size = 3, aes(y = `Intensity [%]`, label = `Ion [m/z]`)) +
  facet_wrap(~ Name, scales = "free", ncol = 2, labeller = "label_parsed") +
  xlab(expression(italic("m/z"))) +
  scale_color_manual(values = c("gray20", inferno(3, begin = .1, end = .9)[2]),
                     name = "Spectrum") +
  # theme_publish(base_size = 12) +
  theme_publish(base_family = font) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave('../../../figures/ms-covers.pdf', scale = 1.5, width = pagewidth,
       height = 4, unit = 'in', device = cairo_pdf)
# ggsave("../Manuscript/figures/covers-specs.pdf",
#        scale = scal, width = 170, height = 100, unit = "mm", device = cairo_pdf)
# ggsave("../Manuscript/figures/covers-specs.jpg", bg = "white",
#        scale = scal, width = 170, height = 100, unit = "mm")
