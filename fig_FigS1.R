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
jpeg(filename = './fig_S1.jpg', width = 8, height = 6.5, units = 'in', res = 500)
layout(mat = matrix(1:6, nrow = 3))
par(mar = c(3.3,3.6,1.1,1.1))
plot(NA, NA, xlim = c(0.01,1), ylim = c(0,min(1.5*max(IAR.bestcase),1)), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
lines(epsilon.vec, IAR.bestcase[1,], lwd = linewidth, col = palette[1])
lines(epsilon.vec, IAR.bestcase[2,], lwd = linewidth, col = palette[1],lty=2)
box()
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = NA)#seq(from = 0, to = 100, by = 20))
axis(side = 2, at = seq(from = 0, to = 1, by = 0.1), labels = seq(from = 0, to = 100, by = 10), las = 1)
## mtext(side = 1, line = 2.3, expression('Reduction in ' * 'R'[0] * ' (%), ' * epsilon))
mtext(side = 2, line = 2.3, 'Infection attack rate (%)')
mtext(side = 3, line = 0, adj = 0, 'A. No biases', font = 2, cex=0.75)
legend('bottomright',
       legend = c('Control arm', 'Treatment arm'),
       lty = c(2,1), lwd=linewidth*c(1,1),
       bty = 'n')


plot(NA, NA, xlim = c(0.01,1), ylim = c(0,min(1.5*max(IAR.bestcase),1)), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
lines(epsilon.vec, IAR.mosquito[1,], lwd = linewidth, col = palette[3])
lines(epsilon.vec, IAR.mosquito[2,], lwd = linewidth, col = palette[3],lty=2)
box()
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = NA)#seq(from = 0, to = 100, by = 20))
axis(side = 2, at = seq(from = 0, to = 1, by = 0.1), labels = seq(from = 0, to = 100, by = 10), las = 1)
## mtext(side = 1, line = 2.3, expression('Reduction in ' * 'R'[0] * ' (%), ' * epsilon))
mtext(side = 2, line = 2.3, 'Infection attack rate (%)')
mtext(side = 3, line = 0, adj = 0, 'B. Mosquito contamination', font = 2, cex=0.75)

plot(NA, NA, xlim = c(0.01,1), ylim = c(0,min(1.5*max(IAR.bestcase),1)), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
lines(epsilon.vec, IAR.human[1,], lwd = linewidth, col = palette[2])
lines(epsilon.vec, IAR.human[2,], lwd = linewidth, col = palette[2],lty=2)
box()
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
axis(side = 2, at = seq(from = 0, to = 1, by = 0.1), labels = seq(from = 0, to = 100, by = 10), las = 1)
mtext(side = 1, line = 2.3, expression('Reduction in ' * 'R'[0] * ' (%), ' * epsilon))
mtext(side = 2, line = 2.3, 'Infection attack rate (%)')
mtext(side = 3, line = 0, adj = 0, 'C. Human movement', font = 2, cex=0.75)

plot(NA, NA, xlim = c(0.01,1), ylim = c(0,min(1.5*max(IAR.bestcase),1)), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
lines(epsilon.vec, IAR.hummoz[1,], lwd = linewidth, col = palette[2])
lines(epsilon.vec, IAR.hummoz[1,], lwd = linewidth, col = palette[3],lty=2)
lines(epsilon.vec, IAR.hummoz[2,], lwd = linewidth, col = palette[2])
lines(epsilon.vec, IAR.hummoz[2,], lwd = linewidth, col = palette[3],lty=2)
lines(epsilon.vec, IAR.hummoz[2,], lwd = linewidth*1.5, col = "white",lty=3)
box()
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = NA)#seq(from = 0, to = 100, by = 20))
axis(side = 2, at = seq(from = 0, to = 1, by = 0.1), labels = seq(from = 0, to = 100, by = 10), las = 1)
## mtext(side = 1, line = 2.3, expression('Reduction in ' * 'R'[0] * ' (%), ' * epsilon))
mtext(side = 2, line = 2.3, 'Infection attack rate (%)')
mtext(side = 3, line = 0, adj = 0, 'D. Human movement and mosquito contamination', font = 2, cex=0.75)

plot(NA, NA, xlim = c(0.01,1), ylim = c(0,min(1.5*max(IAR.bestcase),1)), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
lines(epsilon.vec, IAR.humsupp[1,], lwd = linewidth, col = palette[2])
lines(epsilon.vec, IAR.humsupp[1,], lwd = linewidth, col = palette[4],lty=2)
lines(epsilon.vec, IAR.humsupp[2,], lwd = linewidth, col = palette[2])
lines(epsilon.vec, IAR.humsupp[2,], lwd = linewidth, col = palette[4],lty=2)
lines(epsilon.vec, IAR.humsupp[2,], lwd = linewidth*1.5, col = "white",lty=3)
box()
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = NA)#seq(from = 0, to = 100, by = 20))
axis(side = 2, at = seq(from = 0, to = 1, by = 0.1), labels = seq(from = 0, to = 100, by = 10), las = 1)
## mtext(side = 1, line = 2.3, expression('Reduction in ' * 'R'[0] * ' (%), ' * epsilon))
mtext(side = 2, line = 2.3, 'Infection attack rate (%)')
mtext(side = 3, line = 0, adj = 0, 'E. Human movement and transmission coupling', font = 2, cex=0.75)

plot(NA, NA, xlim = c(0.01,1), ylim = c(0,min(1.5*max(IAR.bestcase),1)), axes = F,
     xlab = '', ylab = '', xaxs = 'i', yaxs = 'i')
lines(epsilon.vec, IAR.fullmodel[1,], lwd = linewidth, col = palette[2])
lines(epsilon.vec, IAR.fullmodel[1,], lwd = linewidth, col = palette[3],lty=2)
lines(epsilon.vec, IAR.fullmodel[1,], lwd = linewidth, col = palette[4],lty=3)
lines(epsilon.vec, IAR.fullmodel[2,], lwd = linewidth, col = palette[2])
lines(epsilon.vec, IAR.fullmodel[2,], lwd = linewidth, col = palette[3],lty=2)
lines(epsilon.vec, IAR.fullmodel[2,], lwd = linewidth, col = palette[4],lty=3)
lines(epsilon.vec, IAR.fullmodel[2,], lwd = linewidth*1.5, col = "white",lty=3)
box()
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
axis(side = 2, at = seq(from = 0, to = 1, by = 0.1), labels = seq(from = 0, to = 100, by = 10), las = 1)
mtext(side = 1, line = 2.3, expression('Reduction in ' * 'R'[0] * ' (%), ' * epsilon))
mtext(side = 2, line = 2.3, 'Infection attack rate (%)')
mtext(side = 3, line = 0, adj = 0, 'F. All biases', font = 2, cex=0.75)

dev.off()
