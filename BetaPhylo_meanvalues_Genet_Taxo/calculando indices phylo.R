#mudando nome de acordo com a filogenia
muda_nome=function(x,sep=" "){
  unlist(lapply(strsplit(x = as.character(x),sep),function(x){
    paste(x[1],x[2],sep="_")
    }))
}

#carregando pacotes usados
library(letsR)
# library(maptools)
library(rgdal)
library(ape)

library(picante)


#Amazon: New data set 2019_Dec
fish_mat=read.csv(paste("/home/murilo/Insync/msdias@unb.br/OneDrive Biz/Amazon",
                        "/AmazonFish_analises/",
                        "StatistiquesSubBasinAmazon_02052016/",
                        "NewData_2019_10/Dataset1_SpeciesList_97bv_311019.csv",
                        sep=""),
                  # paste("/home/murilo/Insync/msdias@unb.br/OneDrive Biz/Amazon",
                  #       "/AmazonFish_analises/",
                  #       "StatistiquesSubBasinAmazon_02052016/",
                  #       "NewData_2017_04/MatrixSpecies_0417.csv",
                  #       sep="")
                  sep=";",header=T)

#leaving only valid species names
dim(fish_mat)
head(fish_mat)
table(fish_mat$Occurrence.Status)
fish_mat<-droplevels(fish_mat[fish_mat$Occurrence.Status %in% "valid",])

fish_mat$Referent.Species.Name<-muda_nome(fish_mat$Referent.Species.Name,sep="[.]")

table(grep(fish_mat$Referent.Species.Name,pattern = "[_]")==1:length(fish_mat$Referent.Species.Name))


# ##excluding additional marine fish families   2022_03_25
# #marine "Atherinopsidae", "Achiridae", "Belonidae", "Clupeidae","Engraulidae",
# #       "Hemiramphidae", "Potamotrygonidae", "Pristigasteridae", "Sciaenidae",
# #       "Tetraodontidae"
# #not marine according to Lovejoy et al.2006 "Atherinopsidae"
# fish_mat<-droplevels(fish_mat[!fish_mat$Family.Referent.Species %in%
#                                 c("Atherinopsidae", "Achiridae", "Belonidae",
#                                   "Clupeidae","Engraulidae", "Hemiramphidae",
#                                   "Potamotrygonidae", "Pristigasteridae", "Sciaenidae",
#                                   "Tetraodontidae"),])
# table(fish_mat$Family.Referent.Species)

####retirar as especies de origem marinha
head(fish_mat)
# marineIncSp<-c("Potamotrygon","Paratrygon","Heliotrygon","Plesiotrygon","Anchoviella","Anchovia",
#   "Lycengraulis","Pterengraulis","Amazonsprattus","Jurengraulis","Pellona","Ilisha","Pristigaster",
#   "Rhinosardina","Potamorhaphis","Belonion","Pseudotylusurus","Hyporhamphus",
#   "Plagioscion","Pachypops","Pachyurus","Petilipinnis","Colomesus")
# table(fish_mat$Genus.Referent.Species[fish_mat$Genus.Referent.Species %in% marineIncSp])
# fish_mat[fish_mat$Genus.Referent.Species%in% marineIncSp,]
# fish_mat<-fish_mat[!fish_mat$Genus.Referent.Species %in% marineIncSp,]


# ##eliminating Orestias species because they are from Titikaka lake
fish_mat<-droplevels(fish_mat[!fish_mat$Genus.Referent.Species %in% "Orestias",])

#unindo os codigos antigos e novos
fish2=read.csv(paste("/home/murilo/Insync/msdias@unb.br/OneDrive Biz/Amazon",
                     "/AmazonFish_analises/",
                     "StatistiquesSubBasinAmazon_02052016/",
                     "NewData_2019_10/stat97bv_311019.csv",
                     sep=""),
               sep=";",header=T)[,c("Sub_drainage","Basin2stat")]
head(fish2)

fish_mat<-merge(x=fish_mat,y=fish2,by.x="Sub_drainage",by.y="Sub_drainage",all=T)
head(fish_mat)
fish_mat$basin<-fish_mat$Basin2stat




fish_mat$occ<-1
library(reshape2)
fish_mat<-dcast(fish_mat,basin~Referent.Species.Name,fun.aggregate = sum, value.var = "occ")
dim(fish_mat)

fish_mat[1:10,1:5]
rownames(fish_mat)<- fish_mat$basin
fish_mat <- fish_mat[,-c(1)]
fish_mat[fish_mat>0]<-1




mami=fish_mat



psim<-list()
psor<-list()
psne<-list()
PDindex<-list()
PCDp<-list()
PCDc<-list()

#phylogeny
filos<-read.tree(file="actinopt_full.trees")#peixes #Rabosky et al. 2020 JBio & Chang et al. 2019 TREE
#filos<-read.tree(file="actinopt_12k_treePL.tre")#peixes


#loop over all posterior samplings
#library(foreach)
#library(doMC)
#registerDoMC(cores=4)

for(i in 1:100){
  #i=1
  #posterior samplings
  filo<-filos[[i]] #fazendo com uma arvore so
  
  #especies presentes tanto na filogenia quanto no arquivo mami
  sp<-intersect(filo$tip.label,colnames(mami));length(sp)
  
  #fazendo a matriz
  mat<-fish_mat[,sp]
  #mat[1:10,1:6]
  #dim(mat)
  
  #deixar na filo as sp da matrix 
  #library(picante)
  filo_gamb<-prune.sample(samp = mat, phylo = filo)
  #plot.phylo(filo_gamb,type="fan",edge.color="gray",cex=0.01)
  
  
  # #phyloBetaSorSim
   library(betapart)
   SIM<-phylo.beta.pair(x = mat,tree = filo_gamb,index.family = "sor")
   psor[[i]]<-as.matrix(SIM$phylo.beta.sor)
   psim[[i]]<-as.matrix(SIM$phylo.beta.sim)
   psne[[i]]<-as.matrix(SIM$phylo.beta.sne)
   print(i) 
}



#loop over all posterior samplings
for(i in 1:100){
  #i=1
  #posterior samplings
  filo<-filos[[i]] #fazendo com uma arvore so
  
  #especies presentes tanto na filogenia quanto no arquivo mami
  sp<-intersect(filo$tip.label,colnames(mami));length(sp)
  
  #fazendo a matriz
  mat<-fish_mat[,sp]
  #mat[1:10,1:6]
  #dim(mat)
  
  #deixar na filo as sp da matrix 
  #library(picante)
  filo_gamb<-prune.sample(samp = mat, phylo = filo)
  #plot.phylo(filo_gamb,type="fan",edge.color="gray",cex=0.5)
  
  
  ##Fast package for PhyloMetrics
  library(PhyloMeasures)
  # #Calculate mntd values for each community
  is.ultrametric(filo_gamb)

  # library(phytools)
  library(phangorn)
  if(!is.ultrametric(filo_gamb)){
    # ## compute the NNLS ultrametric tree
    filo_gamb<-nnls.tree(cophenetic(filo_gamb),filo_gamb,rooted=TRUE)
  }
  ## check
  is.ultrametric(filo_gamb)

  # table(filo_gamb$tip.label %in% colnames(mat))
  # mat <- mat[,filo_gamb$tip.label]


  # # compute Phylogenetic Diversity measures(MNTD,MPD,PD)
  pd_z <- pd.query(filo_gamb,mat,null.model = "uniform",
                 standardize = TRUE, reps = 1999)
  pd_p <- pd.pvalues(filo_gamb,mat,null.model = "uniform",
                 reps = 1999)
  mpd_z <- mpd.query(filo_gamb,mat,null.model = "uniform",
                          standardize = TRUE, reps = 1999)
  mpd_p <- mpd.pvalues(filo_gamb,mat,null.model = "uniform",
                          reps = 1999)
  mntd_z <- mntd.query(filo_gamb,mat,null.model = "uniform",
                            standardize = TRUE, reps = 1999)
  mntd_p <- mntd.pvalues(filo_gamb,mat,null.model = "uniform",
                       reps = 1999)
  SR<-rowSums(mat)
  PD <- data.frame(SR,
                   pd_z,pd_p,
                   mpd_z,mpd_p,
                   mntd_z,mntd_p,
                   basin=rownames(mat))
  PDindex[[i]] <- PD
  
  #Phylogenetic Community Dissimilarity
  # PCD<-pcd(comm = mat, tree = filo, reps = 100)
  # PCDp[[i]]<-as.matrix(PCD$PCDp)
  # PCDc[[i]]<-as.matrix(PCD$PCDc)
  # names(PCDp)[i]=paste0("i=",i)
  # names(PCDc)[i]=paste0("i=",i)
  # write.table(as.matrix(PCDp[[i]]),
  #           file = paste0("PCD/PCDp_",i,"_over100posteriorSampling_2020_05_06.csv"),
  #           sep=";")
  # write.table(as.matrix(PCDc[[i]]),
  #           file = paste0("PCD/PCDc_",i,"_over100posteriorSampling_2020_05_06.csv"),
  #           sep=";")
  # 
  print(i)
}



retorno<-function(x,val){
    if(is(x)=="dist"){ x=as.matrix(x) }
    return(x[,val])
}




#################
#mean and sd values of PD index over all posterior sampling
PDfinal<-sapply(c("SR","pd_z","pd_p","mpd_z","mpd_p","mntd_z","mntd_p"),FUN = function(x){
  rowMeans(sapply(PDindex,retorno,val=x))
})

# PDfinalSD<-sapply(c("SR","pd_z","pd_p","mpd_z","mpd_p","mntd_z","mntd_p"),FUN = function(x){
#   apply(sapply(PDindex,retorno,val=x),1,sd)
# })

rownames(PDfinal)<-rownames(mat)
#rownames(PDfinalSD)<-rownames(mat)
colnames(PDfinalSD)<-paste0(colnames(PDfinalSD),"_sd")
PDfinal<-data.frame(PDfinal,PDfinalSD)
write.csv(PDfinal,file = "PDfinal_index_over100posteriorSampling_1999permut_2021_04_21.csv")
#write.csv(PDfinalSD,file = "PDsdfinal_index_over100posteriorSampling_2020_04_15.csv")
#################



####calculating metrics with picante - 21/04/2021
#loop over all posterior samplings
for(i in 1:100){
  #i=1
  #posterior samplings
  filo<-filos[[i]] #fazendo com uma arvore so
  
  #especies presentes tanto na filogenia quanto no arquivo mami
  sp<-intersect(filo$tip.label,colnames(mami));length(sp)
  
  #fazendo a matriz
  mat<-fish_mat[,sp]
  #mat[1:10,1:6]
  #dim(mat)
  
  #deixar na filo as sp da matrix 
  #library(picante)
  filo_gamb<-prune.sample(samp = mat, phylo = filo)
  #plot.phylo(filo_gamb,type="fan",edge.color="gray",cex=0.5)
  
  
  ##Fast package for PhyloMetrics
  library(PhyloMeasures)
  # #Calculate mntd values for each community
  is.ultrametric(filo_gamb)
  
  # library(phytools)
  library(phangorn)
  if(!is.ultrametric(filo_gamb)){
    # ## compute the NNLS ultrametric tree
    filo_gamb<-nnls.tree(cophenetic(filo_gamb),filo_gamb,rooted=TRUE)
  }
  ## check
  is.ultrametric(filo_gamb)
  
  # table(filo_gamb$tip.label %in% colnames(mat))
  # mat <- mat[,filo_gamb$tip.label]
  
  
  # # compute Phylogenetic Diversity measures(MNTD,MPD,PD)
  pd.result <- ses.pd(tree = filo_gamb,samp = mat,
                        null.model = "richness",include.root = FALSE)
  mpd.result <- ses.mpd(dis = cophenetic(filo_gamb),samp = mat,
                        null.model = "richness")
  mntd.result <- ses.mntd(dis = cophenetic(filo_gamb),samp = mat,
                        null.model = "richness")
  
  
  SR<-rowSums(mat)
  PD <- data.frame(SR,
                   PD_z=pd.result$pd.obs.z,
                   PD_p=pd.result$pd.obs.p,
                   MPD_z=mpd.result$mpd.obs.z,
                   MPD_p=mpd.result$mpd.obs.p,
                   MNTD_z=mntd.result$mntd.obs.z,
                   MNTD_p=mntd.result$mntd.obs.p,
                   basin=rownames(mat))
  PDindex[[i]] <- PD
  
  #Phylogenetic Community Dissimilarity
  # PCD<-pcd(comm = mat, tree = filo, reps = 100)
  # PCDp[[i]]<-as.matrix(PCD$PCDp)
  # PCDc[[i]]<-as.matrix(PCD$PCDc)
  # names(PCDp)[i]=paste0("i=",i)
  # names(PCDc)[i]=paste0("i=",i)
  # write.table(as.matrix(PCDp[[i]]),
  #           file = paste0("PCD/PCDp_",i,"_over100posteriorSampling_2020_05_06.csv"),
  #           sep=";")
  # write.table(as.matrix(PCDc[[i]]),
  #           file = paste0("PCD/PCDc_",i,"_over100posteriorSampling_2020_05_06.csv"),
  #           sep=";")
  # 
  print(i)
}



retorno<-function(x,val){
  if(is(x)=="dist"){ x=as.matrix(x) }
  return(x[,val])
}




#################
#mean and sd values of PD index over all posterior sampling
PDfinal<-sapply(c("PD","SR","mpd.result","mntd.result"),FUN = function(x){
  rowMeans(sapply(PDindex,retorno,val=x))
})






#################
#Phylo Community Dissimilarity
arq<-dir(paste0(getwd(),"/PCD/"))
arq<-arq[grep(pattern = "PCDp_",x = arq)]
length(arq)
lapply(1:54,function(x){
  #x=1
    PCDp[[x]]<<-as.dist(read.table(paste0(getwd(),"/PCD/",arq[x]),sep=";"))
    print(x)
})
lapply(PCDp,function(x){as.matrix(x)[1:5,1:5]})

#PCDp mean values of PCDp over all posterior sampling
PCDpfinal<-sapply(rownames(as.matrix(PCDp[[1]])),FUN = function(x){
    #x="Abuna"
  rowMeans(sapply(PCDp,retorno,val=x))
})

PCDpfinalSD<-sapply(rownames(as.matrix(PCDp[[1]])),FUN = function(x){
  apply(sapply(PCDp,retorno,val=x),1,sd)
})
write.table(PCDpfinal,file = "PCDpMean_over54posteriorSampling_2020_05_08.csv",sep=";")
write.table(PCDpfinalSD,file = "PCDpSd_over54posteriorSampling_2020_05_08.csv",sep=";")


# #PCDc mean values of PCDc over all posterior sampling
# #As Composition doesn't change, PCDc values are all equal and 
# #PCDcSD is zero. So I do not calculate SD values
# PCDcfinal<-sapply(rownames(as.matrix(PCDc[[1]])),FUN = function(x){
#     #x="Abuna"
#   rowMeans(sapply(PCDc,retorno,val=x))
# })
# 
# write.csv(PCDcfinal,file = "PCDcMean_over14posteriorSampling_2020_05_08.csv")




#################
#PhyloBetaSim

#PSIM mean values of psim over all posterior sampling
PSIMfinal<-sapply(rownames(psim[[1]]),FUN = function(x){
  rowMeans(sapply(psim,retorno,val=x))
})

PSIMfinalSD<-sapply(rownames(psim[[1]]),FUN = function(x){
  apply(sapply(psim,retorno,val=x),1,sd)
})
write.csv(PSIMfinal,file = "PSIMmean_over100posteriorSampling_2020_05_12.csv")
write.csv(PSIMfinalSD,file = "PSIMsd_over100posteriorSampling_2020_05_12.csv")

#PSOR mean values of psim over all posterior sampling
PSORfinal<-sapply(rownames(psor[[1]]),FUN = function(x){
  rowMeans(sapply(psor,retorno,val=x))
})

PSORfinalSD<-sapply(rownames(psor[[1]]),FUN = function(x){
  apply(sapply(psor,retorno,val=x),1,sd)
})
write.csv(PSORfinal,file = "PSORmean_over100posteriorSampling_2020_05_12.csv")
write.csv(PSORfinalSD,file = "PSORsd_over100posteriorSampling_2020_05_12.csv")


#PSNE mean values of psim over all posterior sampling
PSNEfinal<-sapply(rownames(psne[[1]]),FUN = function(x){
  rowMeans(sapply(psne,retorno,val=x))
})

PSNEfinalSD<-sapply(rownames(psne[[1]]),FUN = function(x){
  apply(sapply(psne,retorno,val=x),1,sd)
})
write.csv(PSNEfinal,file = "PSNEmean_over100posteriorSampling_2020_05_12.csv")
write.csv(PSNEfinalSD,file = "PSNEsd_over100posteriorSampling_2020_05_12.csv")


#################


PD<-PDfinal
PD$basin<-rownames(PDfinal)

# #other: same results
# ses.mpd.result <- ses.mpd(mat, cophenetic(filo_gamb), 
#                           null.model = "taxa.labels",
#                           abundance.weighted = FALSE, runs = 99)
# ses.mntd.result <- ses.mntd(mat, cophenetic(filo_gamb), 
#                           null.model = "taxa.labels",
#                           abundance.weighted = FALSE, runs = 99)
# PD<-data.frame(ses.mpd.result,ses.mntd.result,basin=rownames(mat))


# # #merging fish data base
glob@data = data.frame(glob@data, PD[match(glob@data[,"BvNiv2"], PD[,"basin"]),])
glob@data[,20:25]
plot(glob$NbEspeceSi,glob$SR);(abline(a=0,b=1))
summary((glob$SR/glob$NbEspeceSi)[-30])
glob@data[30,]

#mudando nome dos slots dos poligonos
new_IDs = as.vector(glob$BvNiv2)
for (i in 1:length(slot(glob, "polygons"))){
  slot(slot(glob, "polygons")[[i]], "ID") = new_IDs[i]
}

library(ggplot2)
glob_df=fortify(glob)
glob_df<-merge(x=glob_df,y=PD,by.x="id",by.y="basin",all=T)
head(glob_df)


WorldData <- map_data('world')
WorldData<-fortify(WorldData)
wld<-ggplot()
wld<-wld+geom_map(data=WorldData,map=WorldData,
                  aes(x=long, y=lat, map_id=region),
                  color="white", fill="gray60", size=0.2)+xlim(-81,-50)+ylim(-21,6)

library(viridis)
ditch_the_axes <- theme(
  axis.text = element_blank(),
  #axis.line = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_rect(fill = FALSE),
  #panel.grid = element_blank(),
  axis.title = element_blank()
)
wld
wld+geom_polygon(data=glob_df,aes(long,lat,group=group,fill=mpd.result),size=0.2)+
  scale_fill_viridis(begin = 1, end = 0,name = "PD",guide = guide_colorbar(barwidth = 10, barheight = 1.5))+
  ditch_the_axes+
  theme(legend.position = c(0.27, 0.07),legend.direction = "horizontal")


plot(mpd.result~WaterColor,glob@data)
plot(mntd.result~WaterColor,glob@data)
