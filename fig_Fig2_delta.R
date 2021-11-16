# set working directory
setwd('~/Documents/awed_trial_modeling/')

# clear existing workspace
rm(list = ls())

# install necessary packages
if(!require(seqinr)){install.packages('seqinr'); library(seqinr)}

# load output
load('./fig_2_delta.RData')

# specify palette 
palette <- c('#222222',
             '#F8CF2C',
             '#9F2424',
             '#4FA9BB')
linewidth <- 2.5

# generate plot 
jpeg(filename = './fig_2_delta.jpg', width = 8, height = 6.5, units = 'in', res = 500)

layout(mat = matrix(1:2, nrow = 2))
par(mar = c(3.3,3.6,1.1,1.1))
plot(NA, NA, xlim = c(0,1e4), ylim = c(0,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
abline(v = seq(from = 0, to = 1e4, by = 2e3),
       h = seq(from = 0, to = 1, by = 0.2),
       col = col2alpha('gray', alpha = 0.5), lwd = 1.25, lty = 3)
lines(delta.vec, efficacy.bestcase - efficacy.fullmodel, lwd = linewidth, col = palette[1])
lines(delta.vec, (2*(efficacy.bestcase - efficacy.human) + efficacy.mosquito - efficacy.hummoz)/3,
      lwd = linewidth, col = palette[2])
lines(delta.vec, (efficacy.human - efficacy.hummoz + efficacy.bestcase - efficacy.mosquito + efficacy.humsupp - efficacy.fullmodel)/3, lwd = linewidth, col = palette[3])
lines(delta.vec, (efficacy.human - efficacy.humsupp + 2*(efficacy.hummoz - efficacy.fullmodel))/3,
      lwd = linewidth, col = palette[4])
print(approxfun(delta.vec, efficacy.bestcase)(1e3))
print(approxfun(delta.vec, efficacy.mosquito)(1e3))
print(approxfun(delta.vec, efficacy.human)(1e3))
print(approxfun(delta.vec, efficacy.hummoz)(1e3))
print(approxfun(delta.vec, efficacy.humsupp)(1e3))
print(approxfun(delta.vec, efficacy.fullmodel)(1e3))
box()
axis(side = 1, at = seq(from = 0, to = 1e4, by = 1e3), labels = seq(from = 0, to = 1e4, by = 1e3))
axis(side = 2, las = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 1, by = 0.20))
mtext(side = 2, line = 2.3, 'Total bias')
mtext(side = 3, line = 0, adj = 0, 'A', font = 2)
legend('topright', #pch = c(15,15,15,15), pt.cex = 1.5,
       col = palette[c(1,2,3,4)],
       legend = c('None', 'Human movement', 'Mosquito contamination',
                  'Transmission coupling'),
       lty = c(1,1,1,1), lwd=linewidth,
       title = expression(underline('Sources of bias')), bty = 'n')

frac.bias.mosquito <-  (efficacy.bestcase + efficacy.human + efficacy.humsupp - efficacy.mosquito - efficacy.hummoz - efficacy.fullmodel) / 3 / (efficacy.bestcase - efficacy.fullmodel)
frac.bias.suppression <- (2 * efficacy.hummoz + efficacy.human - efficacy.humsupp - 2 * efficacy.fullmodel) / 3 / (efficacy.bestcase - efficacy.fullmodel)
frac.bias.human <- 1 - frac.bias.suppression - frac.bias.mosquito
plot(NA, NA, xlim = c(0,1e4), ylim = c(0,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
polygon(c(delta.vec, rev(delta.vec)), c(c(NA,rep(0,length(delta.vec)-1)), rev(frac.bias.human)),
        border = NA, col = palette[2])
polygon(c(delta.vec, rev(delta.vec)), c(frac.bias.human, rev(frac.bias.human + frac.bias.mosquito)),
        border = NA, col = palette[3])
polygon(c(delta.vec, rev(delta.vec)), c(frac.bias.human + frac.bias.mosquito, rev(frac.bias.mosquito + frac.bias.human + frac.bias.suppression)),
        border = NA, col = palette[4])
box()
print(approxfun(delta.vec, frac.bias.mosquito)(1e3))
print(approxfun(delta.vec, frac.bias.human)(1e3))
print(approxfun(delta.vec, frac.bias.suppression)(1e3))
print(approxfun(frac.bias.fullmodel,delta.vec)(0.1))
print(approxfun(frac.bias.fullmodel,delta.vec)(0.2))
axis(side = 1, at = seq(from = 0, to = 1e4, by = 1e3), labels = seq(from = 0, to = 1e4, by = 1e3))
axis(side = 2, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20), las = 1)
mtext(side = 1, line = 2.3, expression('Width of cluster (m), ' * delta))
mtext(side = 2, line = 2.3, 'Contribution to bias (%)')
mtext(side = 3, line = 0, adj = 0, 'B', font = 2)
text(delta.vec[2],0.09,"Human\nmovement",pos=4)
text(delta.vec[2],0.47,"Mosquito\ncontamination",
     pos=4,col="white")
text(delta.vec[2],0.86,"Transmission\ncoupling",
     pos=4,col="white")

dev.off()
