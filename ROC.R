library(HMP)
library(mosaic)
source("Read_hdf5biom.R")
source("TwoSampleNode.R")

AG = read_hdf5_biom("AG/May 2016/otu_table.biom")
AGmeta = read.csv("AG/May 2016/ag-cleaned.csv", header=TRUE, as.is = TRUE)
AG_tree = read_tree_greengenes("AG/May 2016/97_otus.tree")

AG_OTU = AG$data
AG_Taxa = AG$rows

FecesID = AGmeta$BODY_SITE=="UBERON:feces" & AGmeta$BODY_PRODUCT=="UBERON:feces" & AGmeta$BODY_HABITAT=="UBERON:feces"
FecesID1 = unlist(lapply(names(AG_OTU[[1]]), function(x) x %in% AGmeta$X.SampleID[FecesID]))
FecesOTUSum = sapply(AG_OTU, function(x)sum(x[FecesID1] ))

OTUremFeces = FecesOTUSum > sort(FecesOTUSum,decreasing=TRUE)[101]
stopifnot(sum(OTUremFeces)==100)
OTUname = sapply(AG_Taxa[OTUremFeces],function(x)x$id)
Taxa = t(as.data.frame( lapply(AG_Taxa[OTUremFeces],function(x) x$metadata$taxonomy))); rownames(Taxa) = OTUname
for(i in 1:nrow(Taxa))
  for(j in 1:ncol(Taxa)) Taxa[i,j] = substring(Taxa[i,j], 4,1000)
AG_tree1 = prune_taxa(AG_tree$tip.label %in%  OTUname, AG_tree)

DemographicID = AGmeta$RACE=="Caucasian" & AGmeta$SEX=="male" & AGmeta$ECONOMIC_REGION=="Far West"
SampleSubsetID = DemographicID & FecesID
SampleSubsetID1 = unlist(lapply(names(AG_OTU[[1]]), function(x) x %in% AGmeta$X.SampleID[SampleSubsetID]))

OTUBoth = t(as.matrix(as.data.frame(lapply(AG_OTU[OTUremFeces], function(x) x[SampleSubsetID1]))))
OTUBoth = OTUBoth[,1:(floor(ncol(OTUBoth)/2)*2)]
OTUname = sapply(AG_Taxa[OTUremFeces],function(x)x$id)
rownames(OTUBoth) = OTUname

RevSeq = sapply(1:ntaxa(AG_tree1), function(i) which(rownames(OTUBoth)==AG_tree1$tip.label[i]))
OTUBoth = OTUBoth[RevSeq,]; Taxa = Taxa[RevSeq,]
PT = AG_tree1


PPS = GenPhyloStructure(PT)
PhyloChild = PPS$PhyloChild
PhyloParents = PPS$PhyloParents
Descendant= PPS$Descendant
rm(AG, AG_OTU); gc()

GenTM_stat = function(OTSample1, OTSample2){
  NodeCounts = GenNodeCounts(PT,PhyloParents,PhyloChild,Descendant,OTSample1,OTSample2)
  NodeCount1 = NodeCounts$NodeCount1; NodeCount2 = NodeCounts$NodeCount2
  NodeStat = numeric(PT$Nnode)
  for(i in 1:PT$Nnode){
    NodeID = i+ntaxa(PT)
    x = cbind(NodeCount1[PhyloChild[NodeID,1],], NodeCount1[PhyloChild[NodeID,2],])
    y = cbind(NodeCount2[PhyloChild[NodeID,1],], NodeCount2[PhyloChild[NodeID,2],])
    x = x[rowSums(x)!=0,]; y = y[rowSums(y)!=0,]
    NodeStat[i] = qchisq(1-Xmcupo.sevsample(list(x=x,y=y))$`p value`,1)
  }
  TripletStat = numeric(PT$Nnode)
  DoubleStat = numeric(PT$Nnode)
  for(i in 2:PT$Nnode){
    NodeID = i+ntaxa(PT)
    Par = PhyloParents[NodeID]
    DoubleStat[i] = NodeStat[i] + NodeStat[Par-ntaxa(PT)]
    GPar = PhyloParents[Par]
    if(GPar==0) next
    #TripletSep = qchisq(1-MoMPvalues[c(NodeID,Par,GPar)-ntaxa(PT)],1)
    TripletStat[i] = sum(NodeStat[c(NodeID,Par,GPar)-ntaxa(PT)])
    #cat(i, TripletSep, "\n")
  }
  c(max(NodeStat), max(DoubleStat), max(TripletStat))
}
GenTM_Intp = function(OTSample1, OTSample2){
  NodeCounts = GenNodeCounts(PT,PhyloParents,PhyloChild,Descendant,OTSample1,OTSample2)
  NodeCount1 = NodeCounts$NodeCount1; NodeCount2 = NodeCounts$NodeCount2
  MoMPvalues = numeric(PT$Nnode)
  for(i in 1:PT$Nnode){
    NodeID = i+ntaxa(PT)
    x = cbind(NodeCount1[PhyloChild[NodeID,1],], NodeCount1[PhyloChild[NodeID,2],])
    y = cbind(NodeCount2[PhyloChild[NodeID,1],], NodeCount2[PhyloChild[NodeID,2],])
    x = x[rowSums(x)!=0,]; y = y[rowSums(y)!=0,]
    MoMPvalues[i] = Xmcupo.sevsample(list(x=x,y=y))$`p value`
  }
  MoMPvalues
}
ROC_TDM_Intp = function(xx,SubsetSize){
  G1 = sample(1:ncol(OTUBoth),ncol(OTUBoth)/2); G2 = setdiff(1:ncol(OTUBoth),G1)
  Group1 = OTUBoth[,G1]; Group2 = OTUBoth[,G2]
  OTU_shuffle1 = sample(1:nrow(OTUBoth),SubsetSize); OTU_shuffle2 = sample(OTU_shuffle1,SubsetSize)
  if(SubsetSize==2) {OTU_shuffle2 = OTU_shuffle1[c(2,1)]
  }else
    while(all(OTU_shuffle2==OTU_shuffle1) ) OTU_shuffle2 = sample(OTU_shuffle1,SubsetSize) #Force shuffle?
  Group2s = Group2; Group2s[OTU_shuffle1,] = Group2s[OTU_shuffle2,]
  
  DMp_H0 = Xmcupo.sevsample(list(x = t(Group1), y = t(Group2)))$`p value`
  DMp_Ha = Xmcupo.sevsample(list(x = t(Group1), y = t(Group2s)))$`p value`
  c(GenTM_Intp(Group1,Group2),GenTM_Intp(Group1,Group2s))
}
ROC = function(xx,SubsetSize,alpha){
  G1 = sample(1:ncol(OTUBoth),ncol(OTUBoth)/2); G2 = setdiff(1:ncol(OTUBoth),G1)
  Group1 = OTUBoth[,G1]; Group2 = OTUBoth[,G2]
  OTU_shuffle1 = sample(1:nrow(OTUBoth),SubsetSize); OTU_shuffle2 = sample(OTU_shuffle1,SubsetSize)
  if(SubsetSize==2) {OTU_shuffle2 = OTU_shuffle1[c(2,1)]
  }else
    while(all(OTU_shuffle2==OTU_shuffle1) ) OTU_shuffle2 = sample(OTU_shuffle1,SubsetSize) #Force shuffle?
  Group2s = Group2; Group2s[OTU_shuffle1,] = Group2s[OTU_shuffle1,]*(1-alpha)+Group2s[OTU_shuffle2,]*alpha
  
  DMp_H0 = Xmcupo.sevsample(list(x = t(Group1), y = t(Group2)))$`p value`
  DMp_Ha = Xmcupo.sevsample(list(x = t(Group1), y = t(Group2s)))$`p value`
  
  c(DMp_H0,DMp_Ha,GenTM_stat(Group1, Group2),GenTM_stat(Group1, Group2s))
}

ROC_CirclePlot = function(xx,SubsetSize){
  G1 = sample(1:ncol(OTUBoth),ncol(OTUBoth)/2); G2 = setdiff(1:ncol(OTUBoth),G1)
  Group1 = OTUBoth[,G1]; Group2 = OTUBoth[,G2]
  OTU_shuffle1 = sample(1:nrow(OTUBoth),SubsetSize); OTU_shuffle2 = sample(OTU_shuffle1,SubsetSize)
  if(SubsetSize==2) {OTU_shuffle2 = OTU_shuffle1[c(2,1)]
  }else
    while(all(OTU_shuffle2==OTU_shuffle1) ) OTU_shuffle2 = sample(OTU_shuffle1,SubsetSize) #Force shuffle?
  Group2s = Group2; Group2s[OTU_shuffle1,] = Group2s[OTU_shuffle2,]
  
  MoMpvalues = GenTM_Intp(Group1,Group2s)
  CirclePlot(MoMpvalues, "MoM",0.5,tiplabel=FALSE,offset=2,redblue=FALSE,edgelen=FALSE)
}


##AOAS: p-values for DM and TDM under H0
ROC_p = function(xx){
  G1 = sample(1:ncol(OTUBoth),ncol(OTUBoth)/2); G2 = setdiff(1:ncol(OTUBoth),G1)
  Group1 = OTUBoth[,G1]; Group2 = OTUBoth[,G2]
  DMp_H0 = Xmcupo.sevsample(list(x = t(Group1), y = t(Group2)))$`p value`
  DMmerge_p = numeric(3)
  for(i in 5:3){
    OTSample1_merge = apply(Group1,2, function(x)tapply(x,Taxa[,i],sum))
    OTSample2_merge = apply(Group2,2, function(x)tapply(x,Taxa[,i],sum))
    DM_merge = Xmcupo.sevsample(list(x = t(OTSample1_merge), y = t(OTSample2_merge)))
    stopifnot(all(colSums(OTSample1_merge)==colSums(Group1)), all(colSums(OTSample2_merge)==colSums(Group2)))
    DMmerge_p[6-i] = DM_merge$`p value`
  } 
  c(DMp_H0, DMmerge_p, GenTM_Intp(Group1,Group2))
  
}
m = 5000 #AOAS setting
QQpO = mclapply(1:m, ROC_p,mc.cores=4)
QQp = as.matrix(as.data.frame(QQpO))
histogram(QQp[1,], breaks=25, main="",xlab = "P-value",col="gray")
histogram(QQp[2,], breaks=25, main="",xlab = "P-value",col="gray")
histogram(QQp[3,], breaks=25, main="",xlab = "P-value",col="gray")
histogram(QQp[4,], breaks=25, main="",xlab = "P-value",col="gray")
histogram(as.vector(QQp[5:103,]), breaks=25, main="",xlab = "P-value",col="gray")



##AOAS: Sampling on one OTU
ROC_OTU = function(xx,increment){
  G1 = sample(1:ncol(OTUBoth),ncol(OTUBoth)/2); G2 = setdiff(1:ncol(OTUBoth),G1)
  Group1 = OTUBoth[,G1]; Group2 = OTUBoth[,G2]
  OTU_shuffle = sample(1:ntaxa(PT),1)
  Group2s = Group2
  Group2s[OTU_shuffle,] = Group2s[OTU_shuffle,] * (1+increment)
  DMp_H0 = Xmcupo.sevsample(list(x = t(Group1), y = t(Group2)))$`p value`
  DMp_Ha = Xmcupo.sevsample(list(x = t(Group1), y = t(Group2s)))$`p value`
  c(DMp_H0,DMp_Ha,GenTM_stat(Group1, Group2),GenTM_stat(Group1, Group2s))
}
m = 5000
OTUincrement = 1.5
QQ = mclapply(1:m, ROC_OTU,increment=OTUincrement,mc.cores=4)
QQ = as.matrix(as.data.frame(QQ))
colnames(QQ) = NULL
MT = ""
pplist = seq(0,1,length.out=100)
DM_TPR = sapply(pplist,function(pp) sum(QQ[2,]<pp)/m)
DM_FPR = sapply(pplist,function(pp) sum(QQ[1,]<pp)/m)
Tlist = seq(1,50,length.out=200)
TDM_TPR3 = sapply(Tlist,function(TT) sum(QQ[8,]>TT)/m)
TDM_FPR3 = sapply(Tlist,function(TT) sum(QQ[5,]>TT)/m)
TDM_TPR1 = sapply(Tlist,function(TT) sum(QQ[6,]>TT)/m)
TDM_FPR1 = sapply(Tlist,function(TT) sum(QQ[3,]>TT)/m)
plot(DM_FPR, DM_TPR,type="l", xlab="False positive rate", ylab="True positive rate",main=MT)
lines(TDM_FPR3, TDM_TPR3, col="red")
lines(TDM_FPR1, TDM_TPR1, col="blue")
legend(0.45,0.2, c("DM","PhyloDM 3 nodes", "PhyloDM 1 node"), lty=c(1,1),col=c("black","red","blue")) 




##AOAS: Sampling on the internal nodes
nChildren = sapply(1:PT$Nnode, function(x) sum(Descendant[x+ntaxa(PT),]))
ROC_IntNodes = function(xx,increment,minChildren){
  G1 = sample(1:ncol(OTUBoth),ncol(OTUBoth)/2); G2 = setdiff(1:ncol(OTUBoth),G1)
  Group1 = OTUBoth[,G1]; Group2 = OTUBoth[,G2]
  OTU_shuffle = which(Descendant[sample(which(nChildren>minChildren),1)+ntaxa(PT),])
  Group2s = Group2
  Group2s[OTU_shuffle,] = Group2s[OTU_shuffle,] * (1+increment)
  
  DMp_H0 = Xmcupo.sevsample(list(x = t(Group1), y = t(Group2)))$`p value`
  DMp_Ha = Xmcupo.sevsample(list(x = t(Group1), y = t(Group2s)))$`p value`
  
  c(DMp_H0,DMp_Ha,GenTM_stat(Group1, Group2),GenTM_stat(Group1, Group2s))
  
}
OTUincrement =0.75
minOTU = 5
m = 5000
QQ = mclapply(1:m, ROC_IntNodes,increment=OTUincrement,minChildren=minOTU,mc.cores=4)
QQ = as.matrix(as.data.frame(QQ))
MT = ""
pplist = seq(0,1,length.out=100)
DM_TPR = sapply(pplist,function(pp) sum(QQ[2,]<pp)/m)
DM_FPR = sapply(pplist,function(pp) sum(QQ[1,]<pp)/m)
Tlist = seq(1,50,length.out=200)
TDM_TPR3 = sapply(Tlist,function(TT) sum(QQ[8,]>TT)/m)
TDM_FPR3 = sapply(Tlist,function(TT) sum(QQ[5,]>TT)/m)
TDM_TPR1 = sapply(Tlist,function(TT) sum(QQ[6,]>TT)/m)
TDM_FPR1 = sapply(Tlist,function(TT) sum(QQ[3,]>TT)/m)
plot(DM_FPR, DM_TPR,type="l", xlab="False positive rate", ylab="True positive rate",main=MT)
lines(TDM_FPR3, TDM_TPR3, col="red")
lines(TDM_FPR1, TDM_TPR1, col="blue")
legend(0.45,0.2, c("DM","PhyloDM 3 nodes", "PhyloDM 1 node"), lty=c(1,1),col=c("black","red","blue")) 




##AOAS: OTU Power curve
m = 5000
FindTPR_OTU = function(OTUincrement){
  QQ = mclapply(1:m, ROC_OTU,increment=OTUincrement,mc.cores=4)
  QQ = as.matrix(as.data.frame(QQ))
  DM_FPR = quantile(QQ[1,],0.05)
  TDM1_FPR = quantile(QQ[3,],0.95)
  TDM3_FPR = quantile(QQ[5,],0.95)
  c(sum(QQ[2,]<DM_FPR),sum(QQ[6,]>TDM1_FPR), sum(QQ[8,]>TDM3_FPR))/m
}
IncrementList = seq(1,10,by=1)
PowerCurve_OTU = sapply(IncrementList,FindTPR_OTU)
plot(c(0,IncrementList),c(0.05,PowerCurve_OTU[1,]),type="o",pch=20,ylim=c(0,1),xlab="Increment x100%", ylab="Power")
lines(c(0,IncrementList),c(0.05,PowerCurve_OTU[2,]),type="o",pch=20,col="blue")
lines(c(0,IncrementList),c(0.05,PowerCurve_OTU[3,]),type="o",pch=20,col="red")
legend(4.5,0.2, c("DM","PhyloDM 3 nodes", "PhyloDM 1 node"), lty=c(1,1),col=c("black","red","blue")) 



##AOAS: IntNode Power curve
m = 5000
FindTPR_IntNodes = function(OTUincrement,minOTU){
  QQ = mclapply(1:m, ROC_IntNodes,increment=OTUincrement, minChildren=minOTU,mc.cores=4)
  QQ = as.matrix(as.data.frame(QQ))
  DM_FPR = quantile(QQ[1,],0.05)
  TDM1_FPR = quantile(QQ[3,],0.95)
  TDM3_FPR = quantile(QQ[5,],0.95)
  c(sum(QQ[2,]<DM_FPR),sum(QQ[6,]>TDM1_FPR), sum(QQ[8,]>TDM3_FPR))/m
}
IncrementList = seq(0.2,2,by=0.2)
PowerCurve_IntNode2 = sapply(IncrementList,FindTPR_IntNodes,minOTU=2)
PowerCurve_IntNode3 = sapply(IncrementList,FindTPR_IntNodes,minOTU=3)
PowerCurve_IntNode5 = sapply(IncrementList,FindTPR_IntNodes,minOTU=5)
plot(c(0,IncrementList),c(0.05,PowerCurve_IntNode2[1,]),type="o",pch=20,ylim=c(0,1),xlim=c(0,max(IncrementList)),xlab="Increment x100%", ylab="Power")
lines(c(0,IncrementList),c(0.05,PowerCurve_IntNode2[2,]),type="o",pch=20,col="blue")
lines(c(0,IncrementList),c(0.05,PowerCurve_IntNode2[3,]),type="o",pch=20,col="red")
legend(0.9,0.2, c("DM","PhyloDM 3 nodes", "PhyloDM 1 node"), lty=c(1,1),col=c("black","red","blue")) 

plot(c(0,IncrementList),c(0.05,PowerCurve_IntNode3[1,]),type="o",pch=20,ylim=c(0,1),xlim=c(0,max(IncrementList)),xlab="Increment x100%", ylab="Power")
lines(c(0,IncrementList),c(0.05,PowerCurve_IntNode3[2,]),type="o",pch=20,col="blue")
lines(c(0,IncrementList),c(0.05,PowerCurve_IntNode3[3,]),type="o",pch=20,col="red")
legend(0.9,0.2, c("DM","PhyloDM 3 nodes", "PhyloDM 1 node"), lty=c(1,1),col=c("black","red","blue")) 


plot(c(0,IncrementList),c(0.05,PowerCurve_IntNode5[1,]),type="o",pch=20,ylim=c(0,1),xlim=c(0,max(IncrementList)),xlab="Increment x100%", ylab="Power")
lines(c(0,IncrementList),c(0.05,PowerCurve_IntNode5[2,]),type="o",pch=20,col="blue")
lines(c(0,IncrementList),c(0.05,PowerCurve_IntNode5[3,]),type="o",pch=20,col="red")
legend(0.9,0.2, c("DM","PhyloDM 3 nodes", "PhyloDM 1 node"), lty=c(1,1),col=c("black","red","blue")) 

