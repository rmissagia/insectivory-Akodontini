library(geomorph)

landmarks <- readland.tps(file = "mandibleland_phyl_semcaenosus.tps", specID="imageID")##17 landmark dataset
Groups<-read.csv("groups_phyl.csv", header = T, sep = ",")#Specimen grouping data file
Spp<-as.factor(Groups[,1])#list of specific groupings

#GPA
GPA<-gpagen(landmarks)
plotAllSpecimens(GPA$coords)
lands2d<- two.d.array(GPA$coords)
Csize2d<- log(GPA$Csize)
GPA3<-arrayspecs(lands2d, 17, 2)
plotOutliers(GPA3, inspect.outliers = TRUE)

#centroid size
csize <- GPA$Csize
ln <- log(csize)
#saving centroid size and log of cs per species
Sppmeancsize<-aggregate(csize, by = list(Spp), FUN = mean) # By species
write.csv(Sppmeancsize, file = "csize_speciesjaw_phyl.csv")
Sppmeancsizelog<-aggregate(ln, by = list(Spp), FUN = mean) # By species
write.csv(Sppmeancsize, file = "logcsize_speciesjaw_phyl.csv")

gdf <- geomorph.data.frame(GPA, species = Spp)

#Take mean landmark positions
Sppmean<-aggregate(lands2d, by = list(Spp), FUN=mean)[,2:35] # By species
Sppmean3d<-arrayspecs(Sppmean, 17, 2, sep=".")

#PCA per species
SppmeanPCA <-  gm.prcomp(Sppmean3d, tol=0)
coords <- SppmeanPCA$x
coords
write.csv(coords, file = "pc_coords.csv")