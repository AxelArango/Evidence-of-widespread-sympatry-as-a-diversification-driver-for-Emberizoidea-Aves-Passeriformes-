library(phytools)
library(caper)

setwd("~/Documents/Sympatry on Emberizoidea_Masterscript/Data/")
###load the single species phylogeny created with "masterscript_biogeobears-rates
embtree<-read.tree("family phy.tre")
###load the table with the transition rates, speciation rates, age and richness for Emberizoidea
speciesmatrix<-read.csv("pglsready_rage.csv",header=T)
#creat a comparative data object with the single tip treep and the all the rates
comprates<-comparative.data(embtree,speciesmatrix,names.col="family",vcv.dim = T)
comprates
###Local extinction rate explains speciation
extinction_pgls<-pgls(BAMM~extinction_rate,comprates,lambda = "ML")
###Dispersal rate explains speciation
dispertion_pgls<-pgls(BAMM~dispersal_rate,comprates,lambda="ML")
###Family age explains speciation
DivTime_pgls<-pgls(BAMM~age,comprates,lambda="ML")
###Family age explains richness
RichTime_pgls<-pgls(total_richness~age,comprates,lambda="ML")
###Richness explains speciation
RichDiv_pgls<-pgls(BAMM~total_richness,comprates,lambda="ML")
###Family age explains extinction rate
extinctionTime_pgls<-pgls(extinction_rate~age,comprates,lambda="ML")
###Speciation explains Richness
DivRich<-pgls(total_richness~BAMM, comprates,lambda="ML")

summary(extinction_pgls)
summary(dispertion_pgls)
summary(DivTime_pgls)
summary(RichTime_pgls)
summary(extinctionTime_pgls)
summary(RichDiv_pgls)

