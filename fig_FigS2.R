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
linewidth <- 2.5

# generate plot 
jpeg(filename = './fig_S2.jpg', width = 8, height = 9, units = 'in', res = 500)
layout(mat = matrix(1:3, nrow = 3))
par(mar = c(3.3,3.6,1.1,1.1))
plot(NA, NA, xlim = c(0.01,1), ylim = c(-1,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
lines(epsilon.vec, -(efficacy.bestcase - efficacy.fullmodel), lwd = linewidth, col = palette[1])
lines(epsilon.vec, efficacy.human - efficacy.bestcase, lwd = linewidth*2, col = palette[2],lty=2)
lines(epsilon.vec, efficacy.hummoz - efficacy.mosquito, lwd = linewidth, col = palette[2])
box()
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = NA)#seq(from = 0, to = 100, by = 20))
axis(side = 2, at = seq(from = -1, to = 1, by = 0.5), labels = seq(from = -100, to = 100, by = 50),
     las = 1)
abline(h=0,col="gray",lty=3)
## mtext(side = 1, line = 2.3, expression('Reduction in ' * 'R'[0] * ' (%), ' * epsilon))
mtext(side = 2, line = 2.3, 'Total bias (%)')
mtext(side = 3, line = 0, adj = 0, 'A. Human movement', font = 2)

plot(NA, NA, xlim = c(0.01,1), ylim = c(-1,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
lines(epsilon.vec, -(efficacy.bestcase - efficacy.fullmodel), lwd = linewidth, col = palette[1])
lines(epsilon.vec, efficacy.mosquito - efficacy.bestcase, lwd = linewidth, col = palette[3])
lines(epsilon.vec, efficacy.hummoz - efficacy.human, lwd = linewidth, col = palette[3],lty=2)
lines(epsilon.vec, efficacy.fullmodel - efficacy.humsupp, lwd = linewidth, col = palette[3],lty=3)
box()
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = NA)#seq(from = 0, to = 100, by = 20))
axis(side = 2, at = seq(from = -1, to = 1, by = 0.5), labels = seq(from = -100, to = 100, by = 50),
     las = 1)
abline(h=0,col="gray",lty=3)
## mtext(side = 1, line = 2.3, expression('Reduction in ' * 'R'[0] * ' (%), ' * epsilon))
mtext(side = 2, line = 2.3, 'Total bias (%)')
mtext(side = 3, line = 0, adj = 0, 'B. Mosquito contamination', font = 2)

plot(NA, NA, xlim = c(0.01,1), ylim = c(-1,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
lines(epsilon.vec, -(efficacy.bestcase - efficacy.fullmodel), lwd = linewidth, col = palette[1])
lines(epsilon.vec, efficacy.humsupp - efficacy.human, lwd = linewidth, col = palette[4])
lines(epsilon.vec, efficacy.fullmodel - efficacy.hummoz, lwd = linewidth*2, col = palette[4],lty=2)
box()
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
axis(side = 2, at = seq(from = -1, to = 1, by = 0.5), labels = seq(from = -100, to = 100, by = 50),
     las = 1)
abline(h=0,col="gray",lty=3)
## mtext(side = 1, line = 2.3, expression('Reduction in ' * 'R'[0] * ' (%), ' * epsilon))
mtext(side = 2, line = 2.3, 'Total bias (%)')
mtext(side = 3, line = 0, adj = 0, 'C. Transmission coupling', font = 2)
mtext(side = 1, line = 2.3, expression('Reduction in ' * 'R'[0] * ' (%), ' * epsilon))

dev.off()
