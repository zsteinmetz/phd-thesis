#### Load packages and functions ####
source("func/init.R")

## Read data
seqtable <- fread("../data/seqtable.csv")

allcompounds <- rbind(
  cbind(read_targets("../data/targets/PE.txt", header = F)),
  cbind(read_targets("../data/targets/PP.txt", header = F)),
  cbind(read_targets("../data/targets/PS.txt", header = F))
)
names(allcompounds)[names(allcompounds) == "Contributor"] <- "Polymer"

repl <- list()

repl$reports <- data.table(read_openchrom("../data/reports/replication/"))

## Data handling
repl$data <- merge(seqtable, repl$reports, by = "File Name")
repl$data[, Conc := as.numeric(gsub(" ug/mL.*", "", `Sample Name`))]

repl$targets <- merge(repl$data[Comments == "quantifiable"],
                      allcompounds, by = "Name")

repl$sum <- repl$targets[, .(Mean = mean(`Integrated Area`),
                             RSD = sd(`Integrated Area`)/mean(`Integrated Area`),
                             `Sample ID` = range(`Sample ID`)),
                         by = .(Polymer, Name)]
repl$sum[, range(RSD*100), by = Polymer]

repl$plot <- repl$targets[Comments == "quantifiable" &
                            !(Name %in% c("1,13-Tetradecadiene", "1,15-Hexadecadiene"))]
repl$sumplot <- repl$sum[!(Name %in% c("1,13-Tetradecadiene", "1,15-Hexadecadiene"))]

## Plotting
ggplot(data = repl$plot, aes(as.numeric(`Sample ID`)-33, `Integrated Area`*10^-6)) +
  geom_ribbon(data = repl$sumplot, alpha = .25,
              aes(y = Mean*10^-6, ymin = 1.05*Mean*10^-6, ymax = 0.95*Mean*10^-6, fill = Name)) +
  geom_hline(data = repl$sumplot, aes(yintercept = Mean*10^-6, color = Name),
             linetype = "dashed") +
  geom_point(size = pt, position = position_dodge(.2),
             aes(color = Name, shape = Name)) +
  facet_wrap(~ Polymer, scales = "free", labeller = polymer_labeller) +
  xlab("Run") +
  ylab(expression("Peak area"~"["%.%~10^6*"]")) +
  ylim(0, NA) +
  scale_x_continuous(expand = c(0,0), breaks = seq(2, 10, 2)) +
  scale_shape_manual(name = "Pyrolysates", values = c(1, 17, 15, 6, 23, 19), labels = pyr) +
  scale_color_viridis(discrete = T, name = "Pyrolysates", labels = pyr) +
  scale_fill_viridis(discrete = T, name = "Pyrolysates", labels = pyr) +
  theme_publish(base_family = font) + theme(legend.text.align = 0)
ggsave('../../../figures/py-repeatability.pdf', scale = 1.5, width = pagewidth,
       height = 2.2, unit = 'in', device = cairo_pdf)

## Check for linear trends
sapply(unique(repl$plot$Name), function(x) {
  mod <- lm(`Integrated Area` ~ I(as.numeric(`Sample ID`)-33), data = repl$targets[Name == x])
  op <- par(mfrow = c(2,2))
  plot(mod)
  par(op)
  summary(mod)
}, simplify = FALSE)

