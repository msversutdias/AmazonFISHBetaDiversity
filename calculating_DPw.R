library(knitr)
library(knitcitations)
library(formatR)
#cleanbib()
options("citation_format" = "pandoc",digits=2)
opts_chunk$set(tidy=T,
               tidy.opts=list(width.cutoff=60))


panel.cor <- function(x, y, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- (cor(x, y))
  txt <- format(round(r,2))
  r <- abs(r)
  text(0.5, 0.5, txt, cex = 4*r)
}
var_expl=function(x){
  if(is.null(x$null)){
    print("This is a LM: R2.adjusted has been printed with coefficients!")
  }else     (paste("PseudoR2",round((x$null-x$deviance)/x$null,2)*100,"%"))
} 

# standardize function (centering and reducing)
Standard <- function(x){ (x-mean(x,na.rm=T))/(2*sd(x,na.rm=T))  }

#mudando nome de acordo com a filogenia
muda_nome=function(x,sep=" "){
  unlist(lapply(strsplit(x = as.character(x),sep),function(x){
    paste(x[1],x[2],sep="_")
  }))
}

#criar categoria de cores de acordo com quartis
catego_cores<-function(g){
  ecd_f<-ecdf(g)
  g_f<-g
  g_f[]  <-NA
  g_f[ecd_f(g)<=0.25]<-1
  g_f[ecd_f(g)>0.25 & ecd_f(g)<=0.5]<-2
  g_f[ecd_f(g)>0.5 & ecd_f(g)<=0.75]<-3
  g_f[ecd_f(g)>0.75]<-4
  return(g_f)
}


#plot function
# y=rnorm(10)
# d=data.frame(y,x=6*y+rnorm(10),z=-1.5*y+rnorm(10),
#              p=6*y+rnorm(10),r=6*y+rnorm(10))
# mod<-lm(y~x+z+p+r,d)
# groups=c(1,2,3,4,4)
plot_effect <- function(x,groups=NULL){
  #x=mod
  #g=c(1,3,3)
  f<-as.data.frame(summary(x)$coefficients)#[-1,]
  f$factor = rownames(f)
  #f<-f[order(f$factor),]
  P<-grep(pattern = "P",x = colnames(f))
  op<-par()$mar
  x_eixo = 1:nrow(f)
  par(mar=c(5,13,0.5,0.5))
  stripchart(f$Estimate~x_eixo,las=1,cex=3,
             xlim=c(min(f$Estimate-1.96*f$`Std. Error`),max(f$Estimate+1.96*f$`Std. Error`)),
             #col=ifelse(f[,P]<=0.05,"black","red"),
             pch=1,yaxt="n",
             ylab="",xlab="Coefficient estimates");abline(v=0,lty=1,lwd=2,col="gray")
  if(!is.null(groups)){
    g=groups
    cores<-gray.colors(n = length(unique(groups)),start = 0,end = 1,alpha = 0.4)
    cores<-cores[groups]
    for(i in seq_along(g)){
      rect(xleft = min(f$Estimate-1.96*f$`Std. Error`)-1, ybottom = i-0.5,
           xright = max(f$Estimate+1.96*f$`Std. Error`)+1, ytop = i+0.5,
           col=cores[i],border = "white")
    }
  }
  #abline(h=x_eixo,lty=2,col="lightgray")
  axis(side = 2,at = x_eixo, labels = f$factor,las=1)
  segments(x0 = f$Estimate-1*f$`Std. Error`, y0 = 1:length(f$factor),
           x1=f$Estimate+1*f$`Std. Error`,y1 = 1:length(f$factor),lwd=6)
  segments(x0 = f$Estimate-1.96*f$`Std. Error`, y0 = 1:length(f$factor),
           x1=f$Estimate+1.96*f$`Std. Error`,y1 = 1:length(f$factor),lwd=3)
  points(y=x_eixo,x=f$Estimate,las=1,cex=3,
         bg=ifelse(f[,P]<=0.05,"black","white"),
         col="black",pch=21)
  par(mar=op)
}
# plot_effect(mod,groups)




fish_mat=read.csv(paste("/home/murilo/Dropbox/Amazon",
                        "/AmazonFish_analises/",
                        "StatistiquesSubBasinAmazon_02052016/",
                        "NewData_2019_10/Dataset1_SpeciesList_97bv_311019.csv",
                        sep=""),
                  # paste("/home/murilo/Dropbox/Amazon",
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
fish_mat[10065,]

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


####checking if marine species are out
# head(fish_mat)
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
fish2=read.csv(paste("/home/murilo/Dropbox/Amazon",
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




#setting NA to 0
fish_mat[is.na(fish_mat)] <- 0




library(vegan)



library(ape)
library(GISTools)
library(RColorBrewer)
library(recluster)

D_fish<- betapart::beta.pair(fish_mat,index.family = "sorensen")
str(D_fish)


# #excluding species per site matrix
#rm(fish_mat)
fish_tree <- (hclust((D_fish$beta.sim),method = 'average'))
# cores<-cutree(fish_tree,h=0.5)
# cores<-rainbow(length(unique(cores)))[cores]
# names(cores)==as.phylo(fish_tree)$tip.label
# cores<-brewer.pal(n = 11, name = 'Paired')[cores]

#Definindo numero de cluster Kreft & Jetz 2010
d2 <- cophenetic(fish_tree)
#transform hclust to phylogenetic tree (ape)
tree=as.phylo(fish_tree)
#inspect explained diversity for different cuts of a tree
expl_div<-recluster.expl.diss(tree,D_fish$beta.sim,maxcl=30)
#expl_div
R2=expl_div$expl.div
clusters=expl_div$nclust#[-length(expl_div$nclust)]
library(segmented)
out.lm<-lm(R2~clusters)
o<-segmented(out.lm,seg.Z=~clusters,psi=list(clusters=c(15)),
             control=seg.control(display=FALSE))
slope(o)
o
k=round(o$psi[,2],0);k
plot(R2~clusters,main=paste("Clusters= ",k,sep=""));abline(v=k,col="red",lwd=3);
plot(o,add=T)
cores <- cutree(fish_tree, k=k)
coresDF<-data.frame(subdrainage=names(cores),cores=cores)


cores<-c(brewer.pal(n=k,name = "Paired"))[cores]
coresDF$cores<-cores

op<-par()
plot(as.phylo(fish_tree),type="f",las=1,cex=1.5,
     tip.col=cores,font=c(2));box();#grid()
mtext("a)", adj=0.01, line=-3, cex=3)

library(GISTools)
library(RColorBrewer)
amazon=rgdal::readOGR(paste("/home/murilo/Dropbox/Amazon",
                            "/AmazonFish_analises/",
                            "StatistiquesSubBasinAmazon_02052016/ShapefileBasin/",
                            sep=""),
                      layer = "BasinNivel2_v052016")


#g60<-cutree(tree = fish_tree,h = 0.5)+1
g60 <- data.frame(basin = as.phylo(fish_tree)$tip.label, Cores=cores)
amazon@data = data.frame(amazon@data, g60[match(amazon@data[,"BvNiv2"], g60[,"basin"]),])
#par(mar=c(1,1,1,1))
maps::map("world",ylim = slot(amazon,"bbox")[2,],
          xlim = slot(amazon,"bbox")[1,],fill=T,
          col=gray(0.9),
          resolution = 100);#box()
mtext("b)", adj=-0.05, line=-1.5, cex=3)
plot(amazon,col=as.vector(amazon$Cores),add=T)
# dev.off()

summary(as.vector(D_fish$beta.sim))
sd(as.vector(D_fish$beta.sim))
# 
library(vegan)
set.seed(55444)
NMDS <- metaMDS((D_fish$beta.sim),k = 2,#try = c(1,10), 
                trymax = 200,
                autotransform = F, noshare = F, 
                previous.best = T,parallel=4)


par(mar=c(1,1,1,1))
plot(scores(NMDS,display="sites"),las=1,col="black",
       pch=21,cex=4,bg=coresDF$cores,axes=F);box()
text(x = -0.5,y = -0.3,paste0("stress = ",round(NMDS$stress,2)),cex=1.5)
#title("c)", adj = 0.01, line = -1,cex=3)
mtext("c)", adj=0.01, line=-3, cex=3)
#dev.off()



par(op)

##cutting tree
library(cluster)
groups <- cutree(fish_tree,h = 0.4)
groups <- data.frame(basin = names(groups), groups)
prin_canal<-groups[groups$basin=="Amazon5",]$groups
prin_canal <- as.vector(groups[groups$groups%in%prin_canal,]$basin)
# plot(scores(NMDS),las=1,type="n")
# text(scores(NMDS),rownames(scores(NMDS)),cex=0.8,col=groups$groups)

#Mean distance from each group to the main channel
library(plyr)
groups <- ddply(.data = groups,.variables = 'groups',.fun = function(x){
  # bas<-c("Paru_Este","Tapajos2")
  # bas<-c("Jari")
  bas<-as.vector(x$basin)
  # bas<-as.vector(groups[groups$groups==6,]$basin)
  if( any(prin_canal%in%bas) ) {
    x$Bsim_m <- ( rep(0,length(bas)) )
    return(x)
    } else {
      Bsim<-as.matrix(betapart::beta.pair(fish_mat[c(prin_canal,bas),],
                                          index.family = "sorensen")$beta.sim)
      resu<-mean(Bsim[prin_canal,bas])
      x$Bsim_m <- ( rep(resu,length(bas)) )
      return(x)
    }
})



#merging fish data base
amazon@data = data.frame(amazon@data, groups[match(amazon@data[,"BvNiv2"], groups[,"basin"]),])
k <- 9

colors <- brewer.pal(k, "YlOrRd")
#op <- par()$mai

vari <- amazon@data$Bsim_m
#jpeg("Groups_species compositionBsim.jpeg",width = 4000,height = 4000,res = 400)
par(mai=c(0.1,0.1,1,0.1))
choropleth(amazon,vari,auto.shading(vari,n=k,cols = colors))
text(coordinates(amazon),labels = amazon@data$basin,cex=0.7)
choro.legend(-78,-15,auto.shading(vari,n=k,cols = colors));
# dev.off()
par(op)


##True Nestedness
library(vegan)
BACIAS <- as.vector(rownames(fish_mat))


rm(op,amazon,k,colors,groups,fish_tree,vari)





PBsim<-read.csv(paste0(getwd(),
  "/BetaPhylo_meanvalues_Genet_Taxo/PSIMmean_over100posteriorSampling_2020_05_12.csv"),
                sep=",",header=TRUE,row.names = 1)
PBsim[1:5,1:5]

PBsim<-as.dist(PBsim)

#plot(metaMDS(PBsim,k=2,autotransform = FALSE))
summary(as.vector(PBsim))
sd(as.vector(PBsim))
plot(D_fish$beta.sim~PBsim)
cor.test(D_fish$beta.sim,PBsim,method = "pearson")





###########################
# #phylogenetic diversity
# #phylogeny
filo<-read.tree(file="actinopt_12k_treePL.tre")#peixes

#especies presentes tanto na filogenia quanto no arquivo mami
sp<-intersect(filo$tip.label,colnames(fish_mat));length(sp)

#excluindo as especies do mapa para reter
#somente as especies presentes na filogenia
fish_mat_philo <- droplevels(fish_mat[,colnames(fish_mat) %in% sp])

#deixar na filo as sp da matrix
library(picante)
filo_gamb<-prune.sample(samp = fish_mat_philo, phylo = filo)
plot.phylo(filo_gamb,type="fan",edge.color="gray",cex=0.5)

library(betapart)
DsorPhylo<-phylo.beta.pair(x = fish_mat_philo,tree = filo_gamb,index.family = "sorensen")

plot(DsorPhylo$phylo.beta.sim,D_fish$beta.sim)
cor.test(DsorPhylo$phylo.beta.sim,D_fish$beta.sim)

###########################



### Calculating Dpw

#TREE GENETICS + POLYTOMIES
library(picante)

Dpw_mean<-file.exists("BetaPhylo_meanvalues_Genet_Taxo/DpwMean_30nullModelControlRichness_over100trees.csv")
# Dpw_mean<-file.exists("BetaPhylo_meanvalues_Genet_Taxo/DpwMean_over100trees_withoutMarineGroups.csv")#without marine
#Dpw_mean<-FALSE

if(Dpw_mean==FALSE){
  filos<-read.tree(file="BetaPhylo_meanvalues_Genet_Taxo/actinopt_full.trees")#Rabosky et al. 2020 JBio & Chang et al. 2019 TREE
  Dpw_poly<-list()

  #loop over all posterior samplings

  for(i in 1:100){
    #i=1
    #posterior samplings
    filo<-filos[[i]] #fazendo com uma arvore so
  
    #especies presentes tanto na filogenia quanto no arquivo mami
    sp<-intersect(filo$tip.label,colnames(fish_mat));length(sp)
  
    #fazendo a matriz
    mat<-fish_mat[,sp]
    #mat[1:10,1:6]
    #dim(mat)
  
    #deixar na filo as sp da matrix 
    #library(picante)
    filo_gamb<-prune.sample(samp = mat, phylo = filo)
    #plot.phylo(filo_gamb,type="fan",edge.color="gray",cex=0.5)
  
    # #Dpw
    # library(PhyloMeasures) 
    # Dpw_poly[[i]]<-cd.query(tree=filo_gamb,mat,standardize=FALSE)
    # rownames(Dpw_poly[[i]])<-colnames(Dpw_poly[[i]])<-rownames(mat)
    
    #Dpw with null model to control for richness 2022_08_24
    system.time({
      Dpw_poly[[i]]<-nullPhyloMetric(mat = mat, phylo = filo_gamb,nreps = 30)
    })
     
    print(i) 
    rm(mat2)
  }

  Dpw_poly

  #mean over all Dpw values
  retorno<-function(x,val){
    return(as.matrix(x)[,val])
  }

  Dpw_poly_mean<-sapply(rownames(as.matrix(Dpw_poly[[1]])),FUN = function(x){
    rowMeans(sapply(Dpw_poly,retorno,val=x))
  })

  write.table(x = Dpw_poly_mean,
              file = "BetaPhylo_meanvalues_Genet_Taxo/DpwMean_30nullModelControlRichness_over100trees.csv",
              sep = ";")

} else {
  Dpw_mean<-as.dist(read.table("BetaPhylo_meanvalues_Genet_Taxo/DpwMean_30nullModelControlRichness_over100trees.csv",
                               sep=";",header=TRUE))
  # Dpw_mean<-as.dist(read.table("BetaPhylo_meanvalues_Genet_Taxo/DpwMean_over100trees_withoutMarineGroups.csv",
  #                              sep=";",header=TRUE))
}



#TREE: GENETICS
dim(fish_mat_philo)
filo_gamb
Dpw_gen<-(cd.query(tree=filo_gamb,
                   matrix.a = fish_mat_philo,
                   standardize=FALSE))
rownames(Dpw_gen)<-rownames(fish_mat_philo)
colnames(Dpw_gen)<-rownames(fish_mat_philo)
Dpw_gen<-as.dist(Dpw_gen)

#Dpw genetic tree vs full tree
plot(Dpw_mean~Dpw_gen,asp=1);abline(a = 0,b = 1)
cor(Dpw_mean,Dpw_gen,method = "spearman")


