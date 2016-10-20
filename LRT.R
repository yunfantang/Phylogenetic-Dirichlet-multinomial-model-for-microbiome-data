library(HMP)
library(mosaic)
source("Read_hdf5biom.R")
source("TwoSampleNode.R")
source("SingleNodeOptimize.R")

memory.limit(size=10000)
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


DemographicID = AGmeta$RACE=="Caucasian" & AGmeta$SEX=="female" & AGmeta$ECONOMIC_REGION=="Southwest" 
SampleSubsetID = DemographicID & FecesID
SampleSubsetID1 = unlist(lapply(names(AG_OTU[[1]]), function(x) x %in% AGmeta$X.SampleID[SampleSubsetID]))

OTUBoth = t(as.matrix(as.data.frame(lapply(AG_OTU[OTUremFeces], function(x) x[SampleSubsetID1]))))
OTUname = sapply(AG_Taxa[OTUremFeces],function(x)x$id)
rownames(OTUBoth) = OTUname

RevSeq = sapply(1:ntaxa(AG_tree1), function(i) which(rownames(OTUBoth)==AG_tree1$tip.label[i]))
OTUBoth = OTUBoth[RevSeq,]; Taxa = Taxa[RevSeq,]
PT = AG_tree1
rm(AG, AG_OTU); gc()


PPS = GenPhyloStructure(PT)
PhyloChild = PPS$PhyloChild
PhyloParents = PPS$PhyloParents
Descendant= PPS$Descendant


#Generate node counts
NodeCount = matrix(0,ntaxa(PT)+PT$Nnode,ncol(OTUBoth) )
NodeCount[1:ntaxa(PT),] = OTUBoth 
for(i in (ntaxa(PT)+1):(PT$Nnode+ntaxa(PT)))
  NodeCount[i,] = colSums(OTUBoth[Descendant[i,],])




############
#Begin LRT testing
LRT_all = dirmult(t(OTUBoth))
DM_loglik = LRT_all$loglik
TDM_loglik = numeric(PT$Nnode)
adj = c(-0.05,0.05); adj_mul = c(-1,1,-2,2,-3,3,-4,4,-5,5)
LRTfail_known = c(); LRTfail = c()
nlopt_algo = c("NLOPT_LD_LBFGS")

for(i in 1:PT$Nnode){
  print(i)
  NodeID = i+ntaxa(PT)
  x = cbind(NodeCount[PhyloChild[NodeID,1],], NodeCount[PhyloChild[NodeID,2],])
  DM_node = DM.MoM(x)
  DM_node_loglik= loglik(x,DM_node$gamma)
  if(i %in% LRTfail_known) TDM_loglik[i]=dirmult(x, init=GammaInit[which(LRTfail_known==i), ])$loglik 
  else  TDM_loglik[i]=dirmult(x)$loglik
  if(TDM_loglik[i]<DM_node_loglik) {
    adjflag = FALSE
    for(am in adj_mul){
      gamma_init = sum(DM_node$gamma)*(DM_node$pi+adj*am)
      if(all(gamma_init>0)) {
        temploglik = dirmult(x,init = gamma_init)$loglik
        if(temploglik>DM_node_loglik){
          TDM_loglik[i] = temploglik
          adjflag = TRUE
          break
        }
      }
    }
    if(!adjflag){
      adjflag1 = FALSE
      for(algo in nlopt_algo){
        NewLRT = OptimizeNode(x, algo=algo)
        if(NewLRT[1] > DM_node_loglik && NewLRT[2] < 1e-6){
          TDM_loglik[i] = NewLRT[1]
          adjflag1 = TRUE
          break
        }
      }
      if(!adjflag1){
        LRTfail = c(LRTfail,i)
        TDM_loglik[i] = loglik(x,DM_node$gamma)  
      }
    }
  }
}

cat(format(2*(sum(TDM_loglik) - DM_loglik), nsmall=4))


