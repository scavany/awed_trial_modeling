# set working directory
setwd('~/Documents/awed_trial_modeling/')

# clear existing workspace
rm(list = ls())

# install necessary packages
if(!require(seqinr)){install.packages('seqinr'); library(seqinr)}

# load output
load('./fig_2_HSM.RData')

# specify palette 
palette <- c('#222222',
             '#F8CF2C',
             '#9F2424',
             '#4FA9BB')

# generate plot 
jpeg(filename = './fig_2_HSM.jpg', width = 8, height = 6.5, units = 'in', res = 500)
layout(mat = matrix(1:2, nrow = 2))
par(mar = c(3.3,3.6,1.1,1.1))
plot(NA, NA, xlim = c(0.01,1), ylim = c(0,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
abline(h = seq(from = 0, to = 1, by = 0.2),
       v = seq(from = 0, to = 1, by = 0.2),
       col = col2alpha('gray', alpha = 0.5), lwd = 1.25, lty = 3)
lines(epsilon.vec, efficacy.bestcase, lwd = 1.5, col = palette[1])
lines(epsilon.vec, efficacy.mosquito, lwd = 1.5, col = palette[2])
lines(epsilon.vec, efficacy.human, lwd = 1.5, col = palette[3])
lines(epsilon.vec, efficacy.fullmodel, lwd = 1.5, col = palette[4])
print(approxfun(epsilon.vec, efficacy.bestcase)(0.849))
print(approxfun(epsilon.vec, efficacy.mosquito)(0.849))
print(approxfun(epsilon.vec, efficacy.human)(0.849))
print(approxfun(epsilon.vec, efficacy.fullmodel)(0.849))
box()
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
axis(side = 2, las = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
mtext(side = 2, line = 2.3, 'Observed efficacy (%)')
mtext(side = 3, line = 0, adj = 0, 'A', font = 2)

legend('bottomright', pch = 15, pt.cex = 1.5, col = palette,
       legend = c('None', 'Mosquito contamination', 'Human movement', 'Pathogen suppression'),
       title = expression(underline('Source of bias')), bty = 'n')

plot(NA, NA, xlim = c(0.01,1), ylim = c(0,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
polygon(c(epsilon.vec, rev(epsilon.vec)), c(c(NA,rep(0,length(epsilon.vec)-1)), rev(frac.bias.mosquito)),
        border = NA, col = palette[2])
polygon(c(epsilon.vec, rev(epsilon.vec)), c(frac.bias.mosquito, rev(frac.bias.mosquito + frac.bias.human)),
        border = NA, col = palette[3])
polygon(c(epsilon.vec, rev(epsilon.vec)), c(frac.bias.mosquito + frac.bias.human, rev(frac.bias.mosquito + frac.bias.human + frac.bias.fullmodel)),
        border = NA, col = palette[4])
box()
print(approxfun(epsilon.vec, frac.bias.mosquito)(0.849))
print(approxfun(epsilon.vec, frac.bias.human)(0.849))
print(approxfun(epsilon.vec, frac.bias.fullmodel)(0.849))
print(approxfun(frac.bias.fullmodel,epsilon.vec)(0.1))
print(approxfun(frac.bias.fullmodel,epsilon.vec)(0.2))
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
axis(side = 2, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20), las = 1)
mtext(side = 1, line = 2.3, expression('Reduction in ' * 'R'[0] * ' (%), ' * epsilon))
mtext(side = 2, line = 2.3, 'Contribution to observed bias (%)')
mtext(side = 3, line = 0, adj = 0, 'B', font = 2)

dev.off()
