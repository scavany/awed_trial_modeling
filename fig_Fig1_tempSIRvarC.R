# set working directory
setwd('~/Documents/awed_trial_modeling/')

# clear existing workspace
rm(list = ls())

# install necessary packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(seqinr, contoureR, imager, fields)

# load the data 
load('./fig_1_tempSIRvarC.RData',verbose=TRUE)

palette <- c('#4FA9BB',
             'paleturquoise')

efficacy.implied.contour <- expand.grid(x = 1:nrow(efficacy.implied),
                                        y = 1:ncol(efficacy.implied))
efficacy.implied.contour$z <- apply(efficacy.implied.contour,1,function(xx){ efficacy.implied[ xx[1],xx[2] ]} )
efficacy.implied.contour$x <- rho.implied.vec[efficacy.implied.contour$x]
efficacy.implied.contour$y <- epsilon.implied.vec[efficacy.implied.contour$y]

efficacy.implied.contour.df <- getContourLines(efficacy.implied.contour, levels = c(0.653,0.771,0.849))

# generate plot 
##pdf(file = './fig_1.pdf', width = 10 * (5/6), height = 5)
jpeg(filename = './fig_1_tempSIRvarC.jpg', width = 10 * (5/6), height = 5,units="in",res=1000)
layout(mat = matrix(1:6, nrow = 2, byrow = T))
par(mar = c(2,0.1,1.3,0.1))
checkerboard.im <- load.image("./checkerboard_diagram.png")
plot(checkerboard.im,axes=FALSE)
mtext(side = 3, line = 0, adj = 0.04, 'A', font = 2)
##text(200.6,-70,"A",cex=1.5)

#plot.new()
#mtext(side = 3, line = 0, adj = 0, 'B', font = 2)
par(mar = c(3.3,3.6,1.3,1.2))
plot(NA, NA, xlim = c(0,max(b.vec)), ylim = c(0.5,1), axes = F,
     xaxs = 'i', yaxs = 'i', xlab = '', ylab = '')
abline(h = seq(from = 0.5, to = 1, by = 0.1),
       v = seq(from = 0, to = 1000, by = 200),
       col = col2alpha('gray', alpha = 0.5), lwd = 1.25, lty = 3)
lines(b.vec, rho.vec.by.b, lwd = 1.5, col = '#222222')
box()
axis(side = 1)
axis(side = 2, las = 1, at = seq(from = 0.5, to = 1, by = 0.1), labels = seq(from = 50, to = 100, by = 10))
mtext(side = 1, line = 2.3, expression('Scale of human movement (m), ' * italic('b')))
mtext(side = 2, line = 2.3, expression('Time in allocated arm (%), ' * rho))
mtext(side = 3, line = 0, adj = 0, 'B', font = 2)

plot(NA, NA, xlim = c(0.75,1), ylim = c(0.5,1), axes = F,
     xaxs = 'i', yaxs = 'i', xlab = '', ylab = '')
abline(h = seq(from = 0.5, to = 1, by = 0.1),
       v = seq(from = 0.75, to = 1, by = 0.05),
       col = col2alpha('gray', alpha = 0.5), lwd = 1.25, lty = 3)
polygon(c(1,rev(efficacy.implied.contour.df$x[efficacy.implied.contour.df$z == 0.653]),0.75,
          efficacy.implied.contour.df$x[efficacy.implied.contour.df$z == 0.849]),
        c(0.5,rev(efficacy.implied.contour.df$y[efficacy.implied.contour.df$z == 0.653]),1,
          efficacy.implied.contour.df$y[efficacy.implied.contour.df$z == 0.849]),
        col = col2alpha(palette[2]), border = NA)
contour(rho.implied.vec, epsilon.implied.vec,
        efficacy.implied,levels = c(0.653,0.771,0.849),
        lty = 1, lwd = 1.5, add=  T, drawlabels = F,
        col = c(palette[2], palette[1], palette[2]))
print(min(rho.implied.vec[efficacy.implied[,length(epsilon.implied.vec)] > 0.771])) # minimum rho to possibly observe efficacy in trial
box()
axis(side = 1, at = seq(from = 0.75, to = 1, by = 0.05), labels = seq(from = 75, to = 100, by = 5))
axis(side = 2, las = 1, at = seq(from = 0.5, to = 1, by = 0.1), labels = seq(from = 50, to = 100, by = 10))
mtext(side = 1, line = 2.3, expression('Time in allocated arm (%), ' * rho))
mtext(side = 2, line = 2.3, expression(epsilon * ' needed for observed efficacy'))
legend('topright', pch = 15, col = palette,
       legend = c('Mean', '95% CI'), pt.cex = 1.5, bty = 'n')
mtext(side = 3, line = 0, adj = 0, 'C', font = 2)

plot(NA, NA, xlim = c(0,1.01), ylim = c(0,1.0), axes = F,
     xaxs = 'i', yaxs = 'i', xlab = '', ylab = '')
abline(h = seq(from = 0, to = 1, by = 0.2), v = seq(from = 0, to = 1, by = 0.2),
       col = col2alpha('gray', alpha = 0.5), lwd = 1.25, lty = 3)
lines(epsilon.vec, efficacy.epsilon.bestcase, lwd = 1.5, col = '#222222')
lines(epsilon.vec, efficacy.epsilon.fullmodel, lwd = 1.5, col = palette[1])
segments(x0 = epsilon.trial, y0 = rep(0,3), y1 = efficacy.trial, lty = 3, col = '#222222', lwd = 1.5)
segments(x0  = rep(0,3), x1 = epsilon.trial, y0 = efficacy.trial, lty = 3, col = '#222222', lwd = 1.5)
points(epsilon.trial, efficacy.trial, pch = 22, bg = c(palette[2], palette[1], palette[2]), col = '#222222', cex = 1.5)
print(approxfun(efficacy.epsilon.bestcase,epsilon.vec)(efficacy.trial)) # values of epsilon in bestcase to produce trial efficacy
print(approxfun(efficacy.epsilon.fullmodel,epsilon.vec)(efficacy.trial)) # values of epsilon in fullmodel to produce trial efficacy
box()
axis(side = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
axis(side = 2, las = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
mtext(side = 1, line = 2.3, expression('Reduction in ' * 'R'[0] * ' (%), ' * epsilon))
mtext(side = 2, line = 2.3, 'Estimated efficacy (%)')
mtext(side = 3, line = 0, adj = 0, 'D', font = 2)

rho.vec = seq(0.5,1,by=0.001)
plot(NA, NA, xlim = c(0.5,1), ylim = c(0,1), axes = F,
     xaxs = 'i', yaxs = 'i', xlab = '', ylab = '')
abline(h = seq(from = 0, to = 1, by = 0.2), v = seq(from = 0.5, to = 1, by = 0.1),
       col = col2alpha('gray', alpha = 0.5), lwd = 1.25, lty = 3)
polygon(c(rho.vec, rev(rho.vec)), c(efficacy.rho[,1], rev(efficacy.rho[,3])), col = col2alpha(palette[2]), border = NA)
lines(rho.vec, efficacy.rho[,3], lwd = 1.5, col = palette[2])
lines(rho.vec, efficacy.rho[,1], lwd = 1.5, col = palette[2])
lines(rho.vec, efficacy.rho[,2], lwd = 1.5, col = palette[1])
print(efficacy.rho[abs(rho.vec - rho.vec.by.b[70]) < 0.0005,2]) # observed efficacy if b=120m
print((0.771 - efficacy.rho[abs(rho.vec - rho.vec.by.b[70]) < 0.0005,2]) / 0.771) # proportional reduction
box()
axis(side = 1, at = seq(from = 0.5, to = 1, by = 0.1), labels = seq(from = 50, to = 100, by = 10))
axis(side = 2, las = 1, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20))
mtext(side = 1, line = 2.3, expression('Time in allocated arm (%), ' * rho))
mtext(side = 3, line = 0, adj = 0, 'E', font = 2)

deltafn <- approxfun(x = delta.vec, y = efficacy.delta[1,])
efficacy.delta.low <- deltafn(1e3)
deltafn <- approxfun(x = delta.vec, y = efficacy.delta[2,])
efficacy.delta.mean <- deltafn(1e3)
deltafn <- approxfun(x = delta.vec, y = efficacy.delta[3,])
efficacy.delta.high <- deltafn(1e3)
plot(NA, NA, xlim = c(0,max(delta.vec)), ylim = c(0,1), axes = F,
     xaxs = 'i', yaxs = 'i', xlab = '', ylab = '')
abline(h = seq(from = 0, to = 1, by = 0.2), v = seq(from = 0, to = 10000, by = 2000),
       col = col2alpha('gray', alpha = 0.5), lwd = 1.25, lty = 3)
polygon(c(delta.vec, rev(delta.vec)), c(efficacy.delta[1,], rev(efficacy.delta[3,])), 
        col = col2alpha(palette[2]), border = NA)
lines(delta.vec, efficacy.delta[1,], lwd = 1.5, col = palette[2])
lines(delta.vec, efficacy.delta[3,], lwd = 1.5, col = palette[2])
lines(delta.vec, efficacy.delta[2,], lwd = 1.5, col = palette[1])
print(approxfun(delta.vec,efficacy.delta[2,])(c(500,2000))) #efficacy observed when delta==(500m,2000m)
print((0.771 - approxfun(delta.vec,efficacy.delta[2,])(c(500,2000))) / 0.771) # proportion change
segments(x0 = 1000, y0 = 0, y1 = efficacy.delta.mean, lty = 3, col = '#222222', lwd = 1.5)
segments(x0 = 0, x1 = 1000, y0 = efficacy.delta.mean, lty = 3, col = '#222222', lwd = 1.5)
box()
axis(side = 1)
axis(side = 2, at = seq(from = 0, to = 1, by = 0.2), labels = seq(from = 0, to = 100, by = 20), las = 1)
mtext(side = 1, line = 2.3, expression('Width of cluster (m), ' * delta))
mtext(side = 3, line = 0, 'F', adj = 0, font = 2)
dev.off()
