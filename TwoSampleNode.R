#Two-sample comparison, node-by-node
#Input: OTSample1, OTSample2 (OTU on the row), PT, TT
#Output: A vector of p-values and cirle plot on the tree
library(geiger)
library(nloptr)
library(parallel)
FitBB = function(NodeID, NodeCount, method="MLE"){
  nonZeroIndex = NodeCount[NodeID,]!=0 
  NodeOTU = cbind(NodeCount[PhyloChild[NodeID,1],nonZeroIndex], NodeCount[PhyloChild[NodeID,2],nonZeroIndex])
  if(method=="MoM"){
    Est = DM.MoM(NodeOTU)
  }else Est = dirmult(NodeOTU)
  return (c((1-Est$theta)/Est$theta, Est$pi[1],Est$loglik))
}
V1 = function(a,b) {return(b-b*(a*b+1)/(a+1))}
V2 = function(a,b) {return (b*(1-b)/(a+1))}

GenPhyloStructure = function(PT){
  PhyloParents = numeric(PT$Nnode+ntaxa(PT))
  PhyloChild = matrix(0,PT$Nnode+ntaxa(PT),2)
  for(i in 1:nrow(PT$edge)){
    i1 = PT$edge[i,1]
    i2 = PT$edge[i,2]
    PhyloParents[i2] = i1
    if(PhyloChild[i1,1]==0) {PhyloChild[i1,1]=i2}
    else {PhyloChild[i1,2]=i2}
  }
  
  #Calculate descendant matrix:
  Descendant = matrix(FALSE, PT$Nnode+ntaxa(PT), ntaxa(PT))
  for(i in 1:ntaxa(PT)) Descendant[i,i] = TRUE
  processed = logical(PT$Nnode+ntaxa(PT))
  processed[1:ntaxa(PT)] = TRUE
  while(!all(processed)){
    for(i in (ntaxa(PT)+1):(ntaxa(PT)+PT$Nnode)){
      if(all(processed[PhyloChild[i,]])){
        Descendant[i,Descendant[PhyloChild[i,1],]] = TRUE
        Descendant[i,Descendant[PhyloChild[i,2],]] = TRUE
        processed[i] = TRUE
      }
    }
  }
  return(list(PhyloChild=PhyloChild,PhyloParents=PhyloParents,Descendant=Descendant))
}

GenNodeCounts = function(PT, PhyloParents, PhyloChild, Descendant,OTSample1, OTSample2){
  NodeCount1 = matrix(0,ntaxa(PT)+PT$Nnode,ncol(OTSample1) )
  NodeCount2 = matrix(0,ntaxa(PT)+PT$Nnode,ncol(OTSample2))
  NodeCount1[1:ntaxa(PT),] = OTSample1 
  NodeCount2[1:ntaxa(PT),] = OTSample2
  for(i in (ntaxa(PT)+1):(PT$Nnode+ntaxa(PT))){
    NodeCount1[i,] = colSums(OTSample1[Descendant[i,],])
    NodeCount2[i,] = colSums(OTSample2[Descendant[i,],])
  } 
  return (list(NodeCount1=NodeCount1,NodeCount2=NodeCount2))
}


#MoMPvalues = numeric(PT$Nnode)
#MoM_p = list()
#for(i in 1:PT$Nnode){
#  NodeID = i+ntaxa(PT)
#  x = cbind(NodeCount1[PhyloChild[NodeID,1],], NodeCount1[PhyloChild[NodeID,2],])
#  y = cbind(NodeCount2[PhyloChild[NodeID,1],], NodeCount2[PhyloChild[NodeID,2],])
#  x = x[rowSums(x)!=0,]; y = y[rowSums(y)!=0,]
#  Mx = DM.MoM(x); My = DM.MoM(y)
#  Mxy = DM.MoM(rbind(x,y))
#  MoM_p[[i]] = list(H0_sol = Mxy$gamma, Ha_sol=c(Mx$gamma,My$gamma))
#  MoMPvalues[i] = Xmcupo.sevsample(list(x=x,y=y))$`p value`
#}






CirclePlot = function(NodePvalues, method, ncex=1,tiplabel=FALSE,offset=1,redblue=TRUE,edgelen=TRUE,treetype="phylogram"){
  NodeP <<- -log10(NodePvalues)
  NodeP[NodeP>10] <<- 10
  #print(NodeP)
  nC <<- rep("blue",PT$Nnode)
  nC[NodePvalues<0.05] <<- "red"
  par(mar=c(1,1,1,1))
  plot(PT,type=treetype, main=paste(substring(Name_Freq,1,nchar(Name_Freq)-10), method),cex=ncex,show.tip.label=tiplabel,use.edge.length=edgelen)
  
  #plot(PT,type=treetype, cex=ncex,show.tip.label=tiplabel,use.edge.length=edgelen)
  if(redblue) nodelabels(pch=21, col="black", cex=NodeP/offset, bg=nC)
  else nodelabels(pch=21, col="black", cex=NodeP/offset)
}



