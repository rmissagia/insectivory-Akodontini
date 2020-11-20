library(latticeExtra)
library(lattice)
library(viridis)
library(wesanderson)
d2 <- read.csv("alldata9.csv", row.names = 1)
as.data.frame(d2)->d2
str(d2)
PC1 <- d2$PC1
PC2 <- d2$PC2
BF <- d2$res_logBFxlogCS
MA1 <-d2$MADMinc
MA2 <-d2$res_MASMxlogCS
MA3 <-d2$MATMinc
MA4 <-d2$PC1ma
Diet <- d2$Diet

#plot data
xyplot(PC2 ~ PC1, d2, group = Diet, auto.key = TRUE, col=c("blue", "orange"), pch = list(8,19), xlim =c(-0.1,0.1))
xyplot(PC2 ~ PC1, d2, groups=factor(Diet,labels=c("insectivores","omnivores")),pch=c(8,19),auto.key=list(columns=2))

#contour plots
plot_DM <- levelplot(MA1 ~ PC1 * PC2, d2,pch = c(25,21)[Diet],panel = panel.levelplot.points,
                       cex = 1, col.regions = wes_palette("Zissou1", 100, type = "continuous"), col = "black", lwd = 2, ylim=c(-0.075,0.045), title = "Superficial masseter") +
  layer_(panel.2dsmoother(..., n = 200), col.regions=wes_palette("Zissou1", 100, type = "continuous"))
plot_DM

plot_SM <- levelplot(MA2 ~ PC1 * PC2, d2,pch = c(25,21)[Diet],panel = panel.levelplot.points,
                       cex = 1, col.regions = wes_palette("Zissou1", 100, type = "continuous"), col = "black", lwd = 2, ylim=c(-0.075,0.045), title = "Deep masseter") +
  layer_(panel.2dsmoother(..., n = 200), col.regions=wes_palette("Zissou1", 100, type = "continuous"))
plot_SM

plot_TM <- levelplot(MA3 ~ PC1 * PC2, d2,pch = c(25,21)[Diet],panel = panel.levelplot.points,
                       cex = 1, col.regions = wes_palette("Zissou1", 100, type = "continuous"), col = "black", lwd = 2, ylim=c(-0.075,0.045),title = "Temporal") +
  layer_(panel.2dsmoother(..., n = 200), col.regions=wes_palette("Zissou1", 100, type = "continuous"))
plot_TM

plot_BF <- levelplot(BF ~ PC1 * PC2, d2,pch = c(25,21)[Diet],panel = panel.levelplot.points,
                       cex = 1, col.regions = wes_palette("Zissou1", 100, type = "continuous"), col = "black", lwd = 10, ylim=c(-0.075,0.045),title = "Residuals of log_bite forcexlog_cs") +
  layer_(panel.2dsmoother(..., n = 200), col.regions=wes_palette("Zissou1", 100, type = "continuous"))
plot_BF

