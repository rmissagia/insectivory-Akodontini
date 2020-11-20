library(convevol)
library(geiger)

#Read tree
tr<-read.nexus("tree.nexus")

#Read data
d<-as.matrix(read.csv("data.csv", row.names = 1))
name.check(tr, d)
#Convergent_tips
convtips<-c( "Oxymycterus_amazonicus", "Oxymycterus_dasytrichus", "Oxymycterus_delator", "Oxymycterus_hiska", "Oxymycterus_nasutus", "Oxymycterus_paramensis", "Oxymycterus_quaestor", "Oxymycterus_rufus", "Blarinomys_breviceps", "Brucepattersonius_soricinus", "Lenoxus_apicalis", "Scapteromys_aquaticus", "Scapteromys_tumidus", "Juscelinomys_huanchacae")

##Frequency of convergence
N<-convnum(tr, d, convtips, plot = TRUE, ellipse = NULL, plotellipse = NULL)
N
dev.off()

X<-convnumsig(tr, d, convtips, nsim = 1000, ellipse = NULL,
              plot = TRUE, plotellipse = NULL)
X
dev.off()

##Strength of convergence
S<-convratsig(tr, d, convtips, nsim = 1000)
S

##Summaries(100_sim)
#Frequency_of_convergence
as.character(N[[1]])->entries
as.character(X[[1]])->p_values
names<-c("entries", "p_value")
values<-c(entries, p_values)
data.frame(names, values)->summary_frequency
summary_frequency
write.csv(summary_frequency, file="summary_frequency.csv")

#Strength_of_convergence
as.vector(S$ObservedCs)->obs
as.vector(S$Cutoffs)->cutoff
as.vector(S$Pvals)->p_value
c("C1", "C2", "C3", "C4")->metric
summary_strength<- data.frame(metric, obs, cutoff, p_value)
summary_strength
write.csv(summary_strength, file="summary_strength.csv")