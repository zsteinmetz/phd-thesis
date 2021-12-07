#### Load packages and global vars ####
source("func/init.R")

ftir <- covers <- debris <- list()
ftir$asp <- data.table(batch_asp("../data/imaging/ftir/"))

ftir$samlist <- fread("../data/imaging/samples.csv")
ftir$samlist[, `File Name` := paste0(`Image file`, "_FTIR") ]

# get_lib()
spec_lib <- load_lib()

#### Covers ####
# Adjust spectral intensity
covers$asp <- ftir$asp[`File Name` %in% paste0("Image0", 46:52, "_FTIR")]
covers$lib <- covers$asp[,
                         adj_intens(wavenumber, intensity) %>% 
                           smooth_intens() %>% 
                           subtr_bg() %>% 
                           match_spec(library = spec_lib, which = "ftir",
                                      type = "peaks", top_n = 4), 
                         by = `File Name`]
covers$take <- covers$lib[c(1,5,12,15,17,21,25)]

covers$full <- rbind(
  find_spec(sample_name == 494, library = spec_lib, which = "ftir", type = "full"),
  find_spec(sample_name == 584, library = spec_lib, which = "ftir", type = "full"),
  find_spec(sample_name == 607, library = spec_lib, which = "ftir", type = "full")
)

covers$match <- merge(covers$take, covers$full, by = "sample_name",
                      allow.cartesian = T)

covers$specs <- merge(covers$asp,
                      ftir$samlist[, .(Site, `Transect/Cover`, `Field replicate`,
                                       Row, `File Name`)],
                      by = "File Name", allow.cartesian = T)

covers$specs[Site == 1 & `Transect/Cover` == "Fleece", Label := "(a)"]
covers$specs[Site == 2 & `Transect/Cover` == "Mulch", Label := "(b)"]
covers$specs[Site == 5 & `Transect/Cover` == "Fleece", Label := "(c)"]
covers$specs[Site == 5 & `Transect/Cover` == "Perforated foil", Label := "(d)"]
covers$specs[Site == 8 & `Transect/Cover` == "Perforated foil", Label := "(e)"]

covers$specs[, Label := factor(Label)]

covers$match[wavenumber %between% c(650,4000) &
               `File Name` %in% covers$specs[!is.na(Label), unique(`File Name`)],
             Label := factor(`File Name`, labels = c("(a)", "(b)", "(c)", "(d)",
                                                     "(e)"))]

ggplot(covers$specs[Label != ""], aes(wavenumber, intensity)) +
  geom_line(data = covers$match[Label != ""], aes(color = spectrum_identity),
            alpha = .4) +
  geom_area(data = covers$match[Label != ""], aes(fill = spectrum_identity),
            alpha = .4) +
  geom_line(size = 0.4, color = "gray20") +
  facet_wrap(~ Label, scales = "free", ncol = 2) +
  scale_x_reverse() +
  xlab(expression("Wavenumber"~"["*cm^-1*"]")) +
  scale_color_viridis(discrete = T, direction = -1,
                      option = "inferno", begin = .1, end = .9,
                      name = "Spectrum match",
                      labels = c("PP", "PE", "Low-density PE")) +
  scale_fill_viridis(discrete = T, direction = -1,
                     option = "inferno", begin = .1, end = .9,
                     name = "Spectrum match",
                     labels = c("PP", "PE", "Low-density PE")) +
  # theme_publish(base_size = 12) + 
  theme_publish(base_family = font) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave('../../../figures/ftir-covers.pdf', scale = 1.5, width = pagewidth,
       height = 4.5, unit = 'in', device = cairo_pdf)
# ggsave("../Manuscript/figures/covers-ftir.pdf",
#        scale = scal, width = 170, height = 120, unit = "mm", device = cairo_pdf)
# ggsave("../Manuscript/figures/covers-ftir.jpg", bg = "white",
#        scale = scal, width = 170, height = 120, unit = "mm")


#### Debris ####
debris$asp <- ftir$asp[`File Name` %in% paste0("Image0", c("05", "01", 19, 29, 36,
                                                           40), "_FTIR")]
debris$specs <- debris$asp
debris$specs[, Label := factor(`File Name`, labels = c("(a)", "(c)", "(d)", "(b)",
                                                      "(e)", "(f)"))]
debris$specs[, Label := factor(Label, levels = sort(levels(Label)))]

debris$lib <- debris$specs[,
                           adj_intens(wavenumber, intensity) %>% 
                             smooth_intens() %>% 
                             subtr_bg() %>% 
                             match_spec(library = spec_lib, which = "ftir",
                                        type = "peaks", top_n = 4), 
                           by = `File Name`]
debris$take <- debris$lib[`File Name` %in% debris$specs[Label != "",
                                                        unique(`File Name`)]] %>% 
  .[c(4,5,10,13,17,23)]
debris$take[,unique(sample_name)]

debris$full <- rbind(
  find_spec(sample_name == 602, library = spec_lib, which = "ftir", type = "full"),
  find_spec(sample_name == 308, library = spec_lib, which = "ftir", type = "full"),
  find_spec(sample_name == 600, library = spec_lib, which = "ftir", type = "full"),
  find_spec(sample_name == 554, library = spec_lib, which = "ftir", type = "full")
)

debris$match <- merge(debris$take, debris$full, by = "sample_name",
                      allow.cartesian = T)

debris$match[wavenumber %between% c(650,4000) &
               `File Name` %in% debris$specs[!is.na(Label), unique(`File Name`)],
             Label := factor(`File Name`, labels = c("(a)", "(c)", "(d)", "(b)",
                                                     "(e)", "(f)"))]
debris$match[, Label := factor(Label, levels = sort(levels(Label)))]

ggplot(debris$specs[Label != ""], aes(wavenumber, intensity)) +
  geom_line(data = debris$match[Label != ""], aes(color = spectrum_identity),
            alpha = .4) +
  geom_area(data = debris$match[Label != ""], aes(fill = spectrum_identity),
            alpha = .4) +
  geom_line(size = 0.4, color = "gray20") +
  facet_wrap(~ Label, scales = "free", ncol = 2) +
  scale_x_reverse() +
  xlab(expression("Wavenumber"~"["*cm^-1*"]")) +
  scale_color_viridis(discrete = T, direction = -1,
                      option = "inferno", begin = .1, end = .9,
                      name = "Spectrum match",
                      labels = c("Chitin/cotton", "PE", "PS", "(Natural) resin")) +
  scale_fill_viridis(discrete = T, direction = -1,
                     option = "inferno", begin = .1, end = .9,
                     name = "Spectrum match",
                     labels = c("Chitin/cotton", "PE", "PS", "(Natural) resin")) +
  # theme_publish(base_size = 12) + 
  theme_publish(base_family = font) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave('../../../figures/ftir-debris.pdf', scale = 1.5, width = pagewidth,
       height = 4.5, unit = 'in', device = cairo_pdf)
# ggsave("../Manuscript/figures/debris-ftir.pdf",
#        scale = scal, width = 170, height = 110, unit = "mm", device = cairo_pdf)
# ggsave("../Manuscript/figures/debris-ftir.jpg", bg = "white",
#        scale = scal, width = 170, height = 110, unit = "mm")

