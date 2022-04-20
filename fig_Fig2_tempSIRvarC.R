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
linewidth <- 2.5

# generate plot 
jpeg(filename = './fig_2_tempSIRvarC.jpg', width = 8, height = 6.5, units = 'in', res = 500)

layout(mat = matrix(1:2, nrow = 2))
par(mar = c(3.3,3.6,1.1,1.1))
plot(NA, NA, xlim = c(0.01,1), ylim = c(0,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
abline(h = seq(from = 0, to = 1, by = 0.2),
       v = seq(from = 0, to = 1, by = 0.2),
       col = col2alpha('gray', alpha = 0.5), lwd = 1.25, lty = 3)
lines(epsilon.vec, efficacy.bestcase, lwd = linewidth, col = palette[1])
lines(epsilon.vec, efficacy.human, lwd = linewidth, col = palette[2])
lines(epsilon.vec, efficacy.mosquito, lwd = linewidth, col = palette[3],lty="dotted")
lines(epsilon.vec, efficacy.hummoz, lwd = linewidth, col = palette[2])
lines(epsilon.vec, efficacy.hummoz, lwd = linewidth, col = palette[3],lty="dotted")
lines(epsilon.vec, efficacy.humsupp, lwd = linewidth * 2.5, col = palette[4])
lines(epsilon.vec, efficacy.humsupp, lwd = linewidth, col = palette[2])
lines(epsilon.vec, efficacy.fullmodel, lwd = linewidth * 2.5, col = palette[4])
lines(epsilon.vec, efficacy.fullmodel, lwd = linewidth, col = palette[2])
lines(epsilon.vec, efficacy.fullmodel, lwd = linewidth, col = palette[3],lty="dotted")
print(approxfun(epsilon.vec, efficacy.bestcase)(0.703))
print(approxfun(epsilon.vec, efficacy.mosquito)(0.703))
print(approxfun(epsilon.vec, efficacy.human)(0.703))
print(approxfun(epsilon.vec, efficacy.hummoz)(0.703))
print(approxfun(epsilon.vec, efficacy.humsupp)(0.703))
print(approxfun(epsilon.vec, efficacy.fullmodel)(0.703))
box()
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
axis(side = 2, las = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
mtext(side = 2, line = 2.3, 'Estimated efficacy (%)')
mtext(side = 3, line = 0, adj = 0, 'A', font = 2)

legend('bottomright', #pch = c(15,15,15,15), pt.cex = 1.5,
       col = palette[c(1,2,3,2,4,4)],
       legend = c(expression(italic(Eff^(0))), expression(italic(Eff^(h))),
                  expression(italic(Eff^(m))), expression(italic(Eff^(hm))),
                  expression(italic(Eff^(ht))), expression(italic(Eff^(hmt)))),
       lty = c(1,1,3,1,1,1), lwd=linewidth*c(1,1,1,1,2.5,2.5),
       title = expression(underline('Estimate of efficacy')), bty = 'n')
legend('bottomright', #pch = c(15,15,15,15), pt.cex = 1.5,
       col = palette[c(1,2,3,2,2,2)],
       legend = c(expression(italic(Eff^(0))), expression(italic(Eff^(h))),
                  expression(italic(Eff^(m))), expression(italic(Eff^(hm))),
                  expression(italic(Eff^(ht))), expression(italic(Eff^(hmt)))),
       lty = c(1,1,3,1,1,1), lwd=linewidth*c(1,1,1,1,1,1),
       title = expression(underline('Estimate of efficacy')), bty = 'n')
legend('bottomright', #pch = c(15,15,15,15), pt.cex = 1.5,
       col = palette[c(1,2,3,3,2,3)],
       legend = c(expression(italic(Eff^(0))), expression(italic(Eff^(h))),
                  expression(italic(Eff^(m))), expression(italic(Eff^(hm))),
                  expression(italic(Eff^(ht))), expression(italic(Eff^(hmt)))),
       lty = c(1,1,3,3,1,3), lwd=linewidth*c(1,1,1,1,1,1),
       title = expression(underline('Estimate of efficacy')), bty = 'n')

frac.bias.mosquito <-  (efficacy.bestcase + efficacy.human + efficacy.humsupp - efficacy.mosquito - efficacy.hummoz - efficacy.fullmodel) / 3 / (efficacy.bestcase - efficacy.fullmodel)
frac.bias.suppression <- (2 * efficacy.hummoz + efficacy.human - efficacy.humsupp - 2 * efficacy.fullmodel) / 3 / (efficacy.bestcase - efficacy.fullmodel)
frac.bias.human <- 1 - frac.bias.suppression - frac.bias.mosquito
plot(NA, NA, xlim = c(0.01,1), ylim = c(0,1), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
polygon(c(epsilon.vec, rev(epsilon.vec)), c(c(NA,rep(0,length(epsilon.vec)-1)), rev(frac.bias.human)),
        border = NA, col = palette[2])
polygon(c(epsilon.vec, rev(epsilon.vec)), c(frac.bias.human, rev(frac.bias.human + frac.bias.mosquito)),
        border = NA, col = palette[3])
polygon(c(epsilon.vec, rev(epsilon.vec)), c(frac.bias.human + frac.bias.mosquito, rev(frac.bias.mosquito + frac.bias.human + frac.bias.suppression)),
        border = NA, col = palette[4])
box()
print(approxfun(epsilon.vec, frac.bias.mosquito)(0.703))
print(approxfun(epsilon.vec, frac.bias.human)(0.703))
print(approxfun(epsilon.vec, frac.bias.suppression)(0.703))
print(approxfun(frac.bias.fullmodel,epsilon.vec)(0.1))
print(approxfun(frac.bias.fullmodel,epsilon.vec)(0.2))
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
axis(side = 2, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20), las = 1)
mtext(side = 1, line = 2.3, expression('Reduction in ' * 'R'[0] * ' (%), ' * epsilon))
mtext(side = 2, line = 2.3, 'Contribution to bias (%)')
mtext(side = 3, line = 0, adj = 0, 'B', font = 2)
text(epsilon.vec[1],0.07,"Human\nmovement",pos=4)
text(0.18,0.247,"Mosquito movement",
     pos=4,col="white",srt=10.4)
text(epsilon.vec[1],0.6,"Transmission\ncoupling",
     pos=4,col="white")

dev.off()
