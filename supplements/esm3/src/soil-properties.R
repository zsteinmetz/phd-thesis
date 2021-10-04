#### Load packages and global vars ####
source("func/init.R")

#### pH and EC
ph <- ec <- cn <- samples <- list()
ph$raw <- fread("../data/soil-properties/pH.csv")
ec$raw <- fread("../data/soil-properties/EC.csv")
cn$raw <- fread("../data/soil-properties/CN.csv")

ph$sum <- ph$raw[, .(pH = round(-log(mean(exp(-pH))), 1)), by = Site]

samples$raw <- fread("../data/py-gc-ms/field-samples.csv")

#### Site info ####
samples$sum <- samples$raw[Location != "Blank",
                           lapply(.SD, unique), 
                           .SDcols = c("Cover", "Cultivation",
                                       "Location"),
                           by = Site]
samples$sum[, Site := as.numeric(Site)]

#### Soil texture ####
tex <- list()
tex$raw <- fread("../data/soil-properties/texture.csv")

# Data complete?
tex$raw[, .(No = length(unique(Replicate))), by = .(Site)]
tex$raw[, .N, by = .(Site, Replicate)]

tex$dist <- tex$raw[, .(Class = c("Clay", "Silt", "Sand"),
                        Fraction = 
                          texture(Hydrometer, Blank, Time, Temperature,
                                  plot = T)$din["Estimate",] * 100),
                    by = .(Site)]

tex$cast <- tex$dist %>% dcast(Site ~ Class, value.var = 'Fraction')

# Prepare data
tex$st <- tex$cast
names(tex$st)[names(tex$st) != "Site"] <- toupper(names(tex$st)[names(tex$st) != "Site"])

# Get texture class, for example, in accordance with the German
# "Bodenartendiagramm" (DE.BK94.TT)
require(soiltexture)
tex$din <- data.table(Site = tex$st$Site,
                      TT.points.in.classes(tex$st, class.sys = "DE.BK94.TT")) %>%
  melt(id.vars = "Site", variable.name = "Class") %>% 
  .[value != 0, .(Site, Class)] %>% 
  rbind(data.table(Site = c(3, 6), Class = c("Tu3", "Ut4")))

# Summarize data
tex$sum <- tex$cast[, lapply(.SD, round), 
                    .SDcols = c("Clay", "Silt", "Sand"), by = Site] %>% 
  merge(tex$din, by = "Site")
tex$sum[grep("Tu", Class), Texture := "Silty clay"]
tex$sum[grep("Ut", Class), Texture := "Clayey silt"]

#### Carbon and nitrogen ####
cn$melt <- melt(cn$raw, id = c("Site", "Replicate", "Type"), meas = c("C", "N"),
                var = "Element", val = "Percent")

cn$sum <- cn$melt[, .(Mean = mean(Percent), SD = sd(Percent)),
                  by = .(Site, Type, Element)]
cn$all <- cn$sum[, .(Type = "inorganic",
                     Mean = Mean[Type == "total"] - Mean[Type == "organic"],
                     SD = SD[Type == "total"] + SD[Type == "organic"]),
                 by = .(Site, Element)] %>% 
  rbind(cn$sum)

cn$all[Mean < 0, Mean := 0]

cn$fin <- cn$all[Type != "inorganic"] %>% 
  .[Type == "organic", Type := "org"] %>% 
  .[, Var := paste0(Element, Type)] %>% 
  .[Var != "Norg"] %>%
  .[, Mean := round(Mean, 1)] %>% 
  dcast(Site ~ Var, value.var = "Mean")

#### Aggregate all ####
samples$sum[, -c("Cultivation"), with = F] %>% 
  merge(tex$sum[, -c("Texture"), with = F], by = "Site") %>% 
  merge(cn$fin, by = "Site") %>% 
  merge(ph$sum, by = "Site") %>% 
  merge(ec$raw[, .(Site, EC)], by = "Site")
