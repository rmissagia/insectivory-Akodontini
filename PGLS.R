#directory and libraries
library(phytools)
library(caper)
library(geiger)
library(ape)

#Allometry analysis
#read_tree
tree <- read.tree("ako_tree.nwk")

#read_data
data <- read.csv("alldata6.csv", header = T)
data
#comparative_data_file
cd <- comparative.data(data=data, phy=tree, names.col=X, vcv=TRUE, vcv.dim=3, warn.dropped = TRUE)
?comparative.data

#Table 1
#models
m1<-pgls(BF_log~logCS, cd, lambda="ML", param.CI = 0.95)
summary(m1)
plot.pgls(m1)
#check for outliers
res<- residuals(m1, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m1$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m2<-pgls(MADMinc~logCS, cd, lambda="ML", param.CI = 0.95)
summary(m2)
plot.pgls(m2)
#check for outliers
res<- residuals(m2, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m2$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m3<-pgls(MASMinc~logCS, cd, lambda="ML", param.CI = 0.95)
summary(m3)
plot.pgls(m3)
#check for outliers
res<- residuals(m3, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m3$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m4<-pgls(MATMinc~logCS, cd, lambda="ML", param.CI = 0.95)
summary(m4)
plot.pgls(m4)
#check for outliers
res<- residuals(m4, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m4$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

#Table 2
#models
m1<-pgls(BF_log~Diet, cd, lambda="ML", param.CI = 0.95)
summary(m1)
plot.pgls(m1)
#check for outliers
res<- residuals(m1, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m1$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m2<-pgls(MADMinc~Diet, cd, lambda="ML", param.CI = 0.95)
summary(m2)
plot.pgls(m2)
#check for outliers
res<- residuals(m2, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m2$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m3<-pgls(MASMinc~Diet, cd, lambda="ML", param.CI = 0.95)
summary(m3)
plot.pgls(m3)
#check for outliers
res<- residuals(m3, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m3$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m4<-pgls(MATMinc~Diet, cd, lambda="ML", param.CI = 0.95)
summary(m4)
plot.pgls(m4)
#check for outliers
res<- residuals(m4, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m4$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m5<-pgls(logCS~Diet, cd, lambda="ML", param.CI = 0.95)
summary(m5)
plot.pgls(m5)
#check for outliers
res<- residuals(m5, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m5$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

#Table 3
#models
m1<-pgls(BF_log~logCS, cd, lambda="ML", param.CI = 0.95)
summary(m1)
plot.pgls(m1)
#check for outliers
res<- residuals(m1, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m1$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m2<-pgls(BF_log~MADMinc, cd, lambda="ML", param.CI = 0.95)
summary(m2)
plot.pgls(m2)
#check for outliers
res<- residuals(m2, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m2$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m3<-pgls(BF_log~MADMinc+logCS, cd, lambda="ML", param.CI = 0.95)
summary(m3)
plot.pgls(m3)
#check for outliers
res<- residuals(m3, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m3$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m4<-pgls(BF_log~MADMinc:logCS, cd, lambda="ML", param.CI = 0.95)
summary(m4)
plot.pgls(m4)
#check for outliers
res<- residuals(m4, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m4$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m5<-pgls(BF_log~MADMinc*logCS, cd, lambda="ML", param.CI = 0.95)
summary(m5)
plot.pgls(m5)
#check for outliers
res<- residuals(m5, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m5$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m6<-pgls(BF_log~MASMinc, cd, lambda="ML", param.CI = 0.95)
summary(m6)
plot.pgls(m6)
#check for outliers
res<- residuals(m6, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m6$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m7<-pgls(BF_log~MASMinc+logCS, cd, lambda="ML", param.CI = 0.95)
summary(m7)
plot.pgls(m7)
#check for outliers
res<- residuals(m7, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m7$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m8<-pgls(BF_log~MASMinc:logCS, cd, lambda="ML", param.CI = 0.95)
summary(m8)
plot.pgls(m8)
#check for outliers
res<- residuals(m8, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m8$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m9<-pgls(BF_log~MASMinc*logCS, cd, lambda="ML", param.CI = 0.95)
summary(m9)
plot.pgls(m9)
#check for outliers
res<- residuals(m9, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m9$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m10<-pgls(BF_log~MATMinc, cd, lambda="ML", param.CI = 0.95)
summary(m10)
plot.pgls(m10)
#check for outliers
res<- residuals(m10, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m10$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m11<-pgls(BF_log~MATMinc+logCS, cd, lambda="ML", param.CI = 0.95)
summary(m11)
plot.pgls(m11)
#check for outliers
res<- residuals(m11, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m11$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m12<-pgls(BF_log~MATMinc:logCS, cd, lambda="ML", param.CI = 0.95)
summary(m12)
plot.pgls(m12)
#check for outliers
res<- residuals(m11, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m12$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m13<-pgls(BF_log~MATMinc*logCS, cd, lambda="ML", param.CI = 0.95)
summary(m13)
plot.pgls(m13)
#check for outliers
res<- residuals(m13, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m13$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m14<-pgls(BF_log~MATMinc+MADMinc+MASMinc, cd, lambda="ML", param.CI = 0.95)
summary(m14)
plot.pgls(m14)
#check for outliers
res<- residuals(m14, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m14$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m15<-pgls(BF_log~MATMinc:MADMinc:MASMinc, cd, lambda="ML", param.CI = 0.95)
summary(m15)
plot.pgls(m15)
#check for outliers
res<- residuals(m15, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m15$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m16<-pgls(BF_log~MATMinc*MADMinc*MASMinc, cd, lambda="ML", param.CI = 0.95)
summary(m16)
plot.pgls(m16)
#check for outliers
res<- residuals(m16, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m16$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m17<-pgls(BF_log~MATMinc+MADMinc+MASMinc+logCS, cd, lambda="ML", param.CI = 0.95)
summary(m17)
plot.pgls(m17)
#check for outliers
res<- residuals(m17, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m17$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m18<-pgls(BF_log~MADMinc:MASMinc:MASMinc:logCS, cd, lambda="ML", param.CI = 0.95)
summary(m18)
plot.pgls(m18)
#check for outliers
res<- residuals(m18, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m18$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m19<-pgls(BF_log~MADMinc*MASMinc*MASMinc*logCS, cd, lambda="ML", param.CI = 0.95)
summary(m19)
plot.pgls(m19)
#check for outliers
res<- residuals(m19, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m19$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m20<- pgls(BF_log~Diet, cd, lambda="ML", param.CI = 0.95)
summary(m20)
plot.pgls(m20)
#check for outliers
res<- residuals(m20, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m20$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m21<- pgls(BF_log~logCS+Diet, cd, lambda="ML", param.CI = 0.95)
summary(m21)
plot.pgls(m21)
#check for outliers
res<- residuals(m21, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m21$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m22<- pgls(BF_log~logCS:Diet, cd, lambda="ML", param.CI = 0.95)
summary(m22)
plot.pgls(m22)
#check for outliers
res<- residuals(m22, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m22$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m23<- pgls(BF_log~logCS*Diet, cd, lambda="ML", param.CI = 0.95)
summary(m23)
plot.pgls(m23)
#check for outliers
res<- residuals(m23, phylo = TRUE)
res1<- res/sqrt(var(res))[1] #standardises residuals by sqrt of their variance
rownames(res1)<-rownames(m23$residuals) #matches the residuals up with the species names
rownames(res1)[(abs(res)>3)]#gives the names of the outliers

m24<- pgls(BF_log~1, cd, lambda="ML", param.CI = 0.95)
summary(m24)
plot.pgls(m24)

AICc_values<-c(m1$aicc, m2$aicc, m3$aicc, m4$aicc, m5$aicc, m6$aicc, m7$aicc, m8$aicc, m9$aicc, m10$aicc, m11$aicc, m12$aicc, m13$aicc, m14$aicc, m15$aicc, m16$aicc, m17$aicc, m18$aicc, m19$aicc, m20$aicc, m21$aicc, m22$aicc, m23$aicc, m24$aicc)
aicw(AICc_values)->AICc
models<-c("m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8","m9", "m10", "m11", "m12","m13","m14", "m15", "m16", "m17","m18","m19", "m20", "m21", "m22", "m23", "m24")
summary<-data.frame(models, AICc)
summary[order(summary$delta),]->summary
summary

#R code for PGLS regressions on procrustes coordinates (adapted from Navalon et al. 2018)
library(geomorph)

#Read Procrustes coordinates
ako_coords <- read.csv("Proc.coords.csv", sep=",",row.names=1)
ako_coords1 <-(arrayspecs(ako_coords[,1:ncol(ako_coords)], 17, 2))
ako_coords2 <- two.d.array(ako_coords1)
ako_mean <- mshape(ako_coords2)
v.ako_mean <- t(as.vector(ako_mean))
write.csv(v.ako_mean, "ako_mean.csv")

#Read logCS data
log.CS <- read.csv("logCS.csv", sep=",",row.names=1)
log.CS1 <- as.matrix(log.CS)

#Read diet group data
diet.data <- read.csv("UBF_groups.csv", sep=",",row.names=1)
diet.data1 <-as.matrix(diet.data)

diet.data1<-as.factor(diet.data1)
log.CS2<-as.numeric(log.CS1)

#Set specific geomorph environment

ako_data <- geomorph.data.frame(coords = ako_coords2, diet = diet.data1, CS = log.CS2, phy = tree)

#PGLS regressions of size and diet on procrustes data
PGLS.coords_CS <- procD.pgls(coords ~ CS, data = ako_data, phy = phy, iter = 999, print.progress = T, SS.type = "II")
summary(PGLS.coords_CS)

PGLS.coords_CS <- procD.pgls(coords ~ diet, data = ako_data, phy = phy, iter = 999, print.progress = T, SS.type = "II")
summary(PGLS.coords_CS)