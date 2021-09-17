# set working directory
setwd('~/Documents/awed_trial_modeling/')

# clear existing workspace
rm(list = ls())

# install necessary packages
if(!require(seqinr)){install.packages('seqinr'); library(seqinr)}

# load output
load('./fig_2.RData')

# specify palette 
palette <- c('#222222',
             '#F8CF2C',
             '#9F2424',
             '#4FA9BB')

# generate plot 
jpeg(filename = './fig_2_comb.jpg', width = 8, height = 3.25, units = 'in', res = 500)
## layout(mat = matrix(1:2, nrow = 2))
par(mar = c(3.3,3.6,1.1,1.1))
frac.bias.mosquito <-  (efficacy.bestcase + efficacy.human + efficacy.humsupp - efficacy.mosquito - efficacy.hummoz - efficacy.fullmodel) / 3 / (efficacy.bestcase - efficacy.fullmodel)
frac.bias.suppression <- (2 * efficacy.hummoz + efficacy.human - efficacy.humsupp - 2 * efficacy.fullmodel) / 3 / (efficacy.bestcase - efficacy.fullmodel)
frac.bias.human <- 1 - frac.bias.suppression - frac.bias.mosquito
plot(NA, NA, xlim = c(0.01,1), ylim = c(0,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
polygon(c(epsilon.vec, rev(epsilon.vec)), c(c(NA,rep(0,length(epsilon.vec)-1)), rev(frac.bias.human)),
        border = NA, col = palette[2])
polygon(c(epsilon.vec, rev(epsilon.vec)), c(frac.bias.human, rev(frac.bias.mosquito + frac.bias.human)),
        border = NA, col = palette[3])
polygon(c(epsilon.vec, rev(epsilon.vec)), c(frac.bias.mosquito + frac.bias.human, rev(frac.bias.human + frac.bias.mosquito + frac.bias.suppression)),
        border = NA, col = palette[4])
box()
axis(side = 2, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20), las = 1)
mtext(side = 2, line = 2.3, 'Contribution to observed bias (%)')
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
mtext(side = 1, line = 2.3, expression('Reduction in ' * 'R'[0] * ' (%), ' * epsilon))
text(epsilon.vec[2],0.09,"Human\nmovement",pos=4)
text(epsilon.vec[2],0.47,"Mosquito\ncontamination",
     pos=4,col="white")
text(epsilon.vec[2],0.86,"Transmission\ncoupling",
     pos=4,col="white")


dev.off()
