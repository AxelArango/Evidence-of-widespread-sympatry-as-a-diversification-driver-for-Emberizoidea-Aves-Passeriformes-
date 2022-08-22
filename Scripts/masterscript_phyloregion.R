###masterscript phyloregion
library(phytools)
library(phyloregion)
library(letsR)
library(rgdal)
library(raster)
setwd("~/Documents/Sympatry on Emberizoidea_Masterscript/phyloregion")
##load maps and tree
embemaps<-load("~/Documents/Sympatry on Emberizoidea_Masterscript/Data/embemaps_v1.Rdata")####this object was generated extracting all the Emberizoidea species from the BirdLife interational distribution data set
etreex<-read.tree("~/Documents/Sympatry on Emberizoidea_Masterscript/Data/Breeding-phylogeny_MCC.tre")
#cleaning the maps so they only contain resident and breeding ranges with ceirtain native origin
mapsx<-embemaps
migrationa<-which(mapsx$SEASONAL==1)
migrationb<-which(mapsx$SEASONAL==2)
migration<-c(migrationa,migrationb)
mapsx<-mapsx[migration,]
mapsx<-mapsx[-which(mapsx$ORIGIN!=1),]
###transforming the shp files into a PAM
speciesmat1<-polys2comm(mapsx,res=1,species="species")
#matching the species in the phylogeny with the species in the PAM
phylomatch1<-match_phylo_comm(etreex,speciesmat1$comm_dat,delete_empty_rows=T)
#calculating the phylogenetic turn over
resphylo = "phylobeta.RData"
runslow=T
if (runslow)
{
  phylo1<-phylobeta(phylomatch1$comm,phylomatch1$phy)
  phylosim1<-phylo1$phylo.beta.sim
  nas<-which(is.na(phylosim1))
  phylosim1[nas]<-0
  save(phylosim1,file=resphylo)
} else {
  load(resphylo)
}
#calculating best clustering algorythm and optimal phyloregion
#seleccionar linkages
selink="linkage.RData"
runslow=T
if (runslow){
  selink1<-select_linkage(phylosim1)
  save(selink1,file=selink)
} else {
  load(selink)
}

optimal="optimal_bioreg.RData"
runslow=T
if(runslow){
  opbioreg<-optimal_phyloregion(phylosim1) 
  save(opbioreg,file=optimal)
} else {
  load(optimal)
}
###### calculating the phyloregions, in this case with a k of 8
resbioreg="phyloregion_k8.RData"
runslow=T
if(runslow){
  bioreg1<-phyloregion(phylosim1,k=8,shp=speciesmat1$poly_shp,method = "average")
  bioregx<-bioreg1
  s2x<-bioregx$membership$cluster[which(bioregx$membership$cluster==2)]
  s3x<-bioregx$membership$cluster[which(bioregx$membership$cluster==3)]
  bioregx[s2x]<-3
  bioregx[s3x]<-2
  s2x1<-bioregx$shp$cluster[bioregx$shp$cluster=="2"]
  s3x1<-bioregx$shp$cluster[bioregx$shp$cluster=="3"]
  bioregx$shp$cluster[s2x1]<-3
  bioregx$shp$cluster[s3x1]<-2
  ccolor<-bioregx$shp$COLOURS[bioregx$shp$cluster=="2"]
  bclolor<-bioregx$shp$COLOURS[bioregx$shp$cluster=="3"]
  bioregx$shp$COLOURS[bioregx$shp$cluster=="2"]<-bclolor
  bioregx$shp$COLOURS[bioregx$shp$cluster=="3"]<-ccolor
}else{
  load(resbioreg)
}

plot(bioregx)