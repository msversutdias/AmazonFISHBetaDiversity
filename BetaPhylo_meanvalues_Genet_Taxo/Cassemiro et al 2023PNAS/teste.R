library(ape)
tr=read.nexus("fcassemiro-South-American-freshwater-fish-2568adc/3167_Neotrop_Freshwater_Fish_spp_time_tree.tre")

plot(tr,cex=0.1)
tr$tip.label


tr$tip.label<-unlist(lapply(strsplit(tr$tip.label,split = "_"),function(x){
  paste0(x[2:3],collapse = "_")
}))

tr$tip.label
