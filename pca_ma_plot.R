library(FactoMineR)
library(factoextra)
library(phytools)
library(ggpubr)

ma <- read.csv("ma_facto.csv", header = TRUE)

#PCA
pca <- PCA(ma[,2:4], graph = TRUE)
coords<-pca$ind$coord
coords
write.csv(coords, file = "scores_MAinc.csv")
#The output includes the following components:
print(pca)
eig.val <- get_eigenvalue(pca)
eig.val
write.csv(eig.val, file = "eigenvaluePCAma.csv")

#plot
tree <- read.tree("ako_tree.nwk")
data <- read.csv("scores_MAinc.csv", header = TRUE, row.names = 1)
data
name.check(tree,data)
phylomorphospace(tree,data[,c(1,2)],xlab="PC1",ylab="PC2", label = "off")

ggscatter(data, x = "PC1", y = "PC2", color = "Diet", palette = c("#00AFBB", "#E7B800"))+
  stat_chull(aes(color = Diet, fill = Diet), alpha = 0.1, geom = "polygon")
