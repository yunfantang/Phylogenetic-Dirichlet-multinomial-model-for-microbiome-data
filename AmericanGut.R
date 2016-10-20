library(HMP)
source("Read_hdf5biom.R")
source("TwoSampleNode.R")

AG = read_hdf5_biom("AG/May 2016/otu_table.biom")
AGmeta = read.csv("AG/May 2016/ag-cleaned.csv", header=TRUE, as.is = TRUE)
AG_tree = read_tree_greengenes("AG/May 2016/97_otus.tree")

AG_OTU = AG$data
AG_Taxa = AG$rows

nOTU = 100

FecesID = AGmeta$BODY_SITE=="UBERON:feces" & AGmeta$BODY_PRODUCT=="UBERON:feces" & AGmeta$BODY_HABITAT=="UBERON:feces"
FecesID1 = unlist(lapply(names(AG_OTU[[1]]), function(x) x %in% AGmeta$X.SampleID[FecesID]))
FecesOTUSum = sapply(AG_OTU, function(x)sum(x[FecesID1] ))
OTUremFeces = FecesOTUSum > sort(FecesOTUSum,decreasing=TRUE)[nOTU+1]
stopifnot(sum(OTUremFeces)==nOTU)

L3G1_Freq = c("Daily", "Regularly (3-5 times/week)")
L3G2_Freq = c("Never", "Occasionally (1-2 times/week)", "Rarely (less than once/week)")

L7G1_Freq = c("Daily")
L7G2_Freq = c("Never", "Rarely (less than once/week)", "Regularly (3-5 times/week)", "Occasionally (1-2 times/week)")

L1G1_Freq = c("Daily", "Regularly (3-5 times/week)","Occasionally (1-2 times/week)")
L1G2_Freq = c("Never", "Rarely (less than once/week)")


Name_Freq_list = c("FERMENTED_PLANT_FREQUENCY", "FROZEN_DESSERT_FREQUENCY","MILK_CHEESE_FREQUENCY",
                   "POULTRY_FREQUENCY", "SALTED_SNACKS_FREQUENCY","SEAFOOD_FREQUENCY","SUGARY_SWEETS_FREQUENCY","VEGETABLE_FREQUENCY",
                   "FRUIT_FREQUENCY", "HIGH_FAT_RED_MEAT_FREQUENCY","ALCOHOL_FREQUENCY", "WHOLE_GRAIN_FREQUENCY")

 
Name_Freq = Name_Freq_list[2] #Choose a dietary habit
summary(as.factor(AGmeta[,Name_Freq] ))

DemographicID = FecesID
G1 = (AGmeta[,Name_Freq] %in% L3G1_Freq) & (AGmeta$BODY_SITE=="UBERON:feces") & DemographicID
G2 = (AGmeta[,Name_Freq] %in% L3G2_Freq) & (AGmeta$BODY_SITE=="UBERON:feces") & DemographicID
stopifnot(all(AGmeta$BODY_PRODUCT[G1+G2]=="UBERON:feces"), all(AGmeta$BODY_HABITAT[G1+G2]=="UBERON:feces"))
cat(Name_Freq, sum(G1), sum(G2), "\n")

G1Sample = AGmeta$X.SampleID[G1]
G2Sample = AGmeta$X.SampleID[G2]


#SampleDepth = sapply(1:length(AG_OTU[[1]]),function(i) sum(unlist(lapply(AG_OTU,function(x)x[i] ))))

stopifnot(all(unlist(lapply(1:(length(AG_OTU)-1), function(n)all(names(AG_OTU[[n]]==names(AG_OTU[n+1]))) ))))
SampleIdx1 = names(AG_OTU[[1]]) %in% G1Sample
SampleIdx2 = names(AG_OTU[[1]]) %in% G2Sample

OTU1 = lapply(AG_OTU[OTUremFeces], function(x) x[SampleIdx1]);  OTU2 = lapply(AG_OTU[OTUremFeces], function(x) x[SampleIdx2])

OTUname = sapply(AG_Taxa[OTUremFeces],function(x)x$id)
OTUtable1 = t(as.matrix(as.data.frame(OTU1))); rownames(OTUtable1) = OTUname
OTUtable2 = t(as.matrix(as.data.frame(OTU2))); rownames(OTUtable2) = OTUname
Taxa = t(as.data.frame( lapply(AG_Taxa[OTUremFeces],function(x) x$metadata$taxonomy))); rownames(Taxa) = OTUname
for(i in 1:nrow(Taxa))
  for(j in 1:ncol(Taxa)) Taxa[i,j] = substring(Taxa[i,j], 4,1000)

AG_tree1 = prune_taxa(AG_tree$tip.label %in%  OTUname, AG_tree)
RevSeq = sapply(1:ntaxa(AG_tree1), function(i) which(rownames(OTUtable1)==AG_tree1$tip.label[i]))
OTUtable1 = OTUtable1[RevSeq, ]; OTUtable2 = OTUtable2[RevSeq, ]; Taxa = Taxa[RevSeq,]

#rm(AG, AG_OTU); gc()
#rm(AG, AG_OTU, OTU1, OTU2, OTUtable1, OTUtable2, OTSample1, OTSample2); gc()



Or2 = (rowSums(OTUtable1)+rowSums(OTUtable2))>(sum(OTUtable1)+sum(OTUtable2))*0
PT = prune_taxa(Or2, AG_tree1)
OTSample1 = OTUtable1[Or2,]; OTSample2 = OTUtable2[Or2,]
TT = Taxa[Or2,]
stopifnot(all(rownames(OTSample1)==PT$tip.label),all(rownames(OTSample2)==PT$tip.label),all(rownames(TT)==PT$tip.label) )
PT$tip.label= paste(TT[,2],TT[,3],TT[,4], TT[,5],TT[,6])


PPS = GenPhyloStructure(PT)
PhyloChild = PPS$PhyloChild
PhyloParents = PPS$PhyloParents
Descendant= PPS$Descendant
NodeCounts = GenNodeCounts(PT,PhyloParents,PhyloChild,Descendant,OTSample1,OTSample2)
NodeCount1 = NodeCounts$NodeCount1
NodeCount2 = NodeCounts$NodeCount2



#Generate MoM p-values on each node
MoMPvalues = numeric(PT$Nnode)
MoM_p = list()
for(i in 1:PT$Nnode){
  NodeID = i+ntaxa(PT)
  x = cbind(NodeCount1[PhyloChild[NodeID,1],], NodeCount1[PhyloChild[NodeID,2],])
  y = cbind(NodeCount2[PhyloChild[NodeID,1],], NodeCount2[PhyloChild[NodeID,2],])
  x = x[rowSums(x)!=0,]; y = y[rowSums(y)!=0,]
  Mx = DM.MoM(x); My = DM.MoM(y)
  Mxy = DM.MoM(rbind(x,y))
  MoM_p[[i]] = list(H0_sol = Mxy$gamma, Ha_sol=c(Mx$gamma,My$gamma))
  MoMPvalues[i] = Xmcupo.sevsample(list(x=x,y=y))$`p value`
}



#Siginificant triplets are plotted in red
wthreshold = 16.579
C = 1
TripletCol = rep("lightblue", PT$Nnode)
TripletStat = numeric(PT$Nnode)
for(i in 2:PT$Nnode){
  NodeID = i+ntaxa(PT)
  Par = PhyloParents[NodeID]
  GPar = PhyloParents[Par]
  if(GPar==0) next
  TripletSep = qchisq(1-MoMPvalues[c(NodeID,Par,GPar)-ntaxa(PT)],1)
  TripletStat[i] = sum(TripletSep)
  if(TripletStat[i]>=wthreshold) {
    TripletCol[c(NodeID,Par,GPar)-ntaxa(PT)] = "red"
    cat(c(NodeID,Par,GPar), "\n" )
  }
  #cat(i, TripletSep, "\n")
}
NodeP = -log10(MoMPvalues); max(NodeP)
NodeP[NodeP>10] = 10
par(mar=c(1,1,1,1))
#plot(PT,type="phylogram", main=substring(Name_Freq,1,nchar(Name_Freq)-10),cex=0.5,show.tip.label=FALSE,use.edge.length=FALSE,show.node.label=FALSE)
plot(PT,type="phylogram",cex=0.5,show.tip.label=FALSE,use.edge.length=FALSE,show.node.label=FALSE)
nodelabels(pch=21, col="black", cex=NodeP/2, bg=TripletCol)

#plot(PT,type=treetype, cex=ncex,show.tip.label=tiplabel,use.edge.length=edgelen)
#CirclePlot(, "MoM",0.5,tiplabel=FALSE,offset=2,redblue=FALSE,edgelen=FALSE)
#####



#####
#Get the names of significant nodes to report in the paper
TaxaRank = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
NodeTaxa = rep("", PT$Nnode)
for(i in 1:PT$Nnode){
  NodeAllTaxa = TT[Descendant[i+ntaxa(PT),],]
  for(j in 1:6){
    ppt = NodeAllTaxa[,j]
    CommonTaxa = unique(ppt[ppt!=""])
    if(length(CommonTaxa) != 1) break
    NodeTaxa[i] = CommonTaxa
    names(NodeTaxa)[i] = TaxaRank[j]
  }
}
SigNodeTaxa = NodeTaxa[TripletCol=="red"]
r = lapply(TaxaRank, function(tt) unique(SigNodeTaxa[names(SigNodeTaxa)==tt]))
names(r) = TaxaRank
r
#####




#####
#TDM p-value and DM p-values at different taxonomic levels
cat(Name_Freq,max(TripletStat),"\n")
cat(sum(G1),"\t",sum(G2))
TDMpvalues(max(TripletStat))
DM_OTU = Xmcupo.sevsample(list(x = t(OTSample1), y = t(OTSample2)))
cat(DM_OTU$`p value`,"\t")
for(i in 6:3){
  OTSample1_merge = apply(OTSample1,2, function(x)tapply(x,TT[,i],sum))
  OTSample2_merge = apply(OTSample2,2, function(x)tapply(x,TT[,i],sum))
  DM_merge = Xmcupo.sevsample(list(x = t(OTSample1_merge), y = t(OTSample2_merge)))
  stopifnot(all(colSums(OTSample1_merge)==colSums(OTSample1)), all(colSums(OTSample2_merge)==colSums(OTSample2)))
  cat(DM_merge$`p value`, "\t")
}
cat(Name_Freq, sum(G1), sum(G2), "\n")
#####


