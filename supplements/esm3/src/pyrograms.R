#### Load packages and global vars ####
source("func/init.R")

#### Read data ####
seqtable <- fread("../data/py-gc-ms/seqtable.csv")

reports <- data.table(read_openchrom("../data/py-gc-ms/reports/scans/"))
reports <- reports[!is.na(Name), .(`File Name`, Name, RT = `RT (Min)`,
                                   TIC = `Integrated Area`)]

csvs <- list.files("../data/py-gc-ms/chromatograms", full.names = T,
                   pattern = "*Pyro_Plast-.*\\.csv$") %>%
  lapply(fread) %>% 
  rbindlist(idcol = "fileno") %>% 
  setkey(fileno, `RT(milliseconds)`, `RT(minutes) - NOT USED BY IMPORT`, `RI`)

csvs[, `File Name` := factor(
  fileno, labels = list.files("../data/py-gc-ms/chromatograms",
                              pattern = "*Pyro_Plast-.*\\.csv$") %>%
       gsub("\\.csv$", "", .))]

names(csvs)[names(csvs) == "RT(minutes) - NOT USED BY IMPORT"] <- "RT"

pyrograms <- base::merge(seqtable, csvs, by = "File Name")

pyrograms[, TIC := Reduce(`+`, .SD), .SDcols = 26:length(pyrograms)][]

pyrograms[grep("300", Comment), Compounds := "Volatile additives"]
pyrograms[grep("750", Comment), Compounds := "Polymer"]

lets <- c("(a)", "", "(b)", "(c)", "(d)", "(e)")

pyrograms[, Facet := factor(`File Name`, labels = rep(lets, each = 2))]
pyrograms[, Facet := factor(Facet, levels = sort(levels(Facet)))]

pyrograms[grepl("(\\(a\\))|(\\(c\\))", Facet), Cover := "Fleece"]
pyrograms[grepl("(\\(d\\))|(\\(e\\))", Facet), Cover := "Perforated foil"]
pyrograms[grepl("\\(b\\)", Facet), Cover := "Mulch"]

reports[Name == "Oleanitrile", Component := "'Oleonitrile'"]
reports[grepl("dodecan", Name), Component := "'Propyl dodecanoate'"]
reports[grepl("Octadecenamide", Name), Component := "'9-Octadecenamide'"]
reports[grepl("Phenol", Name), Component := "'Di-'*italic('tert')*'-butylphenol'"]

labels <- reports[!is.na(Component),
                  .(Compounds = "Volatile additives",
                    `File Name`, RT, TIC, Component)] %>% 
  unique()
labels[, Facet := factor(`File Name`, labels = lets)]

labels[RT %between% c(28,30) & Component == "'9-Octadecenamide'",
       Component := NA]


ggplot(pyrograms[Facet != ""], aes(RT, TIC)) +
  geom_text_repel(data = labels[Facet != ""], aes(label = Component),
                  size = 3.5, alpha = .66, parse = T,
                  nudge_y = 20000000,
                  nudge_x = .5,
                  point.padding = 0.25,
                  min.segment.length = 0,
                  segment.size = 0.25) +
  geom_line(size = 0.4, color = "gray20") +
  facet_wrap(~ paste(Facet, Compounds), ncol = 2, scales = "free",
             labeller = as_labeller(function(string){
               c("(a) Polymer pyrolysis", "Volatile additives", "(b)", "", "(c)", "",
                     "(d)", "", "(e)", "")
             })) +
  xlab("Retention time") +
  # theme_publish(base_size = 12) +
  theme_publish(base_family = font) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave('../../../figures/py-covers.pdf', scale = 1.5, width = pagewidth,
       height = 7, unit = 'in', device = cairo_pdf)
# ggsave("../Manuscript/figures/covers-pyrograms.pdf",
#        scale = scal, width = 170, height = 180, unit = "mm", device = cairo_pdf)


