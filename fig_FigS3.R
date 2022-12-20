# set working directory
setwd('~/Documents/awed_trial_modeling/')

# clear existing workspace
rm(list = ls())

# install necessary packages
if(!require(seqinr)){install.packages('seqinr'); library(seqinr)}

# load output
load('./fig_2_tempSIRvarC.RData')

# specify palette 
palette <- c('#222222',
             '#F8CF2C',
             '#9F2424',
             '#4FA9BB')

# generate plot 
jpeg(filename = './fig_S3.jpg', width = 8, height = 9.5, units = 'in', res = 500)
layout(mat = matrix(1:3, nrow = 3))
par(mar = c(3.3,3.6,1.1,1.1))
frac.bias.mosquito <-  (efficacy.bestcase - efficacy.mosquito) / (efficacy.bestcase - efficacy.fullmodel)
frac.bias.human <- (efficacy.mosquito - efficacy.hummoz) / (efficacy.bestcase - efficacy.fullmodel)
frac.bias.suppression <- 1 - frac.bias.human - frac.bias.mosquito
plot(NA, NA, xlim = c(0.01,1), ylim = c(0,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
polygon(c(epsilon.vec, rev(epsilon.vec)), c(c(NA,rep(0,length(epsilon.vec)-1)), rev(frac.bias.human)),
        border = NA, col = palette[2])
polygon(c(epsilon.vec, rev(epsilon.vec)), c(frac.bias.human, rev(frac.bias.mosquito + frac.bias.human)),
        border = NA, col = palette[3])
polygon(c(epsilon.vec, rev(epsilon.vec)), c(frac.bias.mosquito + frac.bias.human, rev(frac.bias.human + frac.bias.mosquito + frac.bias.suppression)),
        border = NA, col = palette[4])
box()
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = NA)#seq(from = 0, to = 100, by = 20))
axis(side = 2, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20), las = 1)
## mtext(side = 1, line = 2.3, expression('Reduction in ' * 'R'[0] * ' (%), ' * epsilon))
mtext(side = 2, line = 2.3, 'Contribution to bias (%)')
mtext(side = 3, line = 0, adj = 0, 'A', font = 2)
mid.index <- as.integer(length(epsilon.vec)/2 + 0.5)
text(epsilon.vec[mid.index],
     frac.bias.human[mid.index]/2,
     "Human movement",adj=0.5)
text(epsilon.vec[mid.index],
     frac.bias.human[mid.index] + frac.bias.mosquito[mid.index]/2,
     "Mosquito movement",
     adj=0.5,col="white")
text(epsilon.vec[101],
     frac.bias.human[mid.index] + frac.bias.mosquito[mid.index] + frac.bias.suppression[mid.index]/2,
     "Transmission coupling",
     adj=0.5,col="white")

frac.bias.human <-  (efficacy.bestcase - efficacy.human) / (efficacy.bestcase - efficacy.fullmodel)
frac.bias.mosquito <- (efficacy.human - efficacy.hummoz) / (efficacy.bestcase - efficacy.fullmodel)
frac.bias.suppression <- 1 - frac.bias.mosquito - frac.bias.human
plot(NA, NA, xlim = c(0.01,1), ylim = c(0,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
polygon(c(epsilon.vec, rev(epsilon.vec)), c(c(NA,rep(0,length(epsilon.vec)-1)), rev(frac.bias.human)),
        border = NA, col = palette[2])
polygon(c(epsilon.vec, rev(epsilon.vec)), c(frac.bias.human, rev(frac.bias.mosquito + frac.bias.human)),
        border = NA, col = palette[3])
polygon(c(epsilon.vec, rev(epsilon.vec)), c(frac.bias.mosquito + frac.bias.human, rev(frac.bias.human + frac.bias.mosquito + frac.bias.suppression)),
        border = NA, col = palette[4])
box()
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
axis(side = 2, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20), las = 1)
## mtext(side = 1, line = 2.3, expression('Reduction in ' * 'R'[0] * ' (%), ' * epsilon))
mtext(side = 2, line = 2.3, 'Contribution to bias (%)')
mtext(side = 3, line = 0, adj = 0, 'B', font = 2)
mid.index <- as.integer(length(epsilon.vec)/2 + 0.5)
text(epsilon.vec[mid.index],
     frac.bias.human[mid.index]/2,
     "Human movement",adj=0.5)
text(epsilon.vec[mid.index],
     frac.bias.human[mid.index] + frac.bias.mosquito[mid.index]/2,
     "Mosquito movement",
     adj=0.5,col="white")
text(epsilon.vec[101],
     frac.bias.human[mid.index] + frac.bias.mosquito[mid.index] + frac.bias.suppression[mid.index]/2,
     "Transmission coupling",
     adj=0.5,col="white")

frac.bias.human <-  (efficacy.bestcase - efficacy.human) / (efficacy.bestcase - efficacy.fullmodel)
frac.bias.suppression <- (efficacy.human - efficacy.humsupp) / (efficacy.bestcase - efficacy.fullmodel)
frac.bias.mosquito <- 1 - frac.bias.suppression - frac.bias.human
plot(NA, NA, xlim = c(0.01,1), ylim = c(0,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
polygon(c(epsilon.vec, rev(epsilon.vec)), c(c(NA,rep(0,length(epsilon.vec)-1)), rev(frac.bias.human)),
        border = NA, col = palette[2])
polygon(c(epsilon.vec, rev(epsilon.vec)), c(frac.bias.human, rev(frac.bias.mosquito + frac.bias.human)),
        border = NA, col = palette[3])
polygon(c(epsilon.vec, rev(epsilon.vec)), c(frac.bias.mosquito + frac.bias.human, rev(frac.bias.human + frac.bias.mosquito + frac.bias.suppression)),
        border = NA, col = palette[4])
box()
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
axis(side = 2, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20), las = 1)
mtext(side = 1, line = 2.3, expression('Reduction in ' * 'R'[0] * ' (%), ' * epsilon))
mtext(side = 2, line = 2.3, 'Contribution to bias (%)')
mtext(side = 3, line = 0, adj = 0, 'C', font = 2)
mid.index <- as.integer(length(epsilon.vec)/2 + 0.5)
text(epsilon.vec[mid.index],
     frac.bias.human[mid.index]/2,
     "Human movement",adj=0.5)
text(epsilon.vec[mid.index],
     frac.bias.human[mid.index] + frac.bias.mosquito[mid.index]/2,
     "Mosquito movement",
     adj=0.5,col="white")
text(epsilon.vec[101],
     frac.bias.human[mid.index] + frac.bias.mosquito[mid.index] + frac.bias.suppression[mid.index]/2,
     "Transmission coupling",
     adj=0.5,col="white")

dev.off()
