library(cubature)
library(R2Cuba)
#L_INT_B = 1e-5
dchisqM = function(x,df){
  if(df>1) return (dchisq(x,df))
  else {
    x[x<0] = -1
    x[(x>=0) &(x<1e-9)] = 1e-9
    dchisq(x,df)
  }
}

InitLabel = function(){
  label = numeric(PT$Nnode)
  Clabel = 0
  cS = -1
  for(i in 1:PT$Nnode){
    if(label[i]!=0) next
    cl = PhyloChild[i+ntaxa(PT),1]; cr = PhyloChild[i+ntaxa(PT),2]
    flag = FALSE
    if(cl>ntaxa(PT)){
      if(PhyloChild[cl,1]>ntaxa(PT)){
        Clabel = Clabel + 1
        label[c(i,cl-ntaxa(PT),PhyloChild[cl,1]-ntaxa(PT))] = Clabel
        next
      }else if (PhyloChild[cl,2]>ntaxa(PT)){
        Clabel = Clabel + 1
        label[c(i,cl-ntaxa(PT),PhyloChild[cl,2]-ntaxa(PT))] = Clabel
        next
      }
    }
    if(cr>ntaxa(PT)){
      if(PhyloChild[cr,1]>ntaxa(PT)){
        Clabel = Clabel + 1
        label[c(i,cr-ntaxa(PT),PhyloChild[cr,1]-ntaxa(PT))] = Clabel
        next
      }else if (PhyloChild[cr,2]>ntaxa(PT)){
        Clabel = Clabel + 1
        label[c(i,cr-ntaxa(PT),PhyloChild[cr,2]-ntaxa(PT))] = Clabel
        next
      }
    }
  }
  for(i in 1:PT$Nnode){
    if(label[i]!=0) next
    cl = PhyloChild[i+ntaxa(PT),1]; cr = PhyloChild[i+ntaxa(PT),2]
    flag = FALSE
    if(cl>ntaxa(PT)){
      Clabel = Clabel + 1
      label[c(i,cl-ntaxa(PT))] = Clabel
      next
    }
    if(cr>ntaxa(PT)){
      Clabel = Clabel + 1
      label[c(i,cr-ntaxa(PT))] = Clabel
    next
    }
    Clabel = Clabel + 1
    label[i] = Clabel
  }
  label
}
   
  
InitGroup3 = function(label){
  Group3 = matrix(0, PT$Nnode-3,3)
  Cntr = 0
  for(i in 3:PT$Nnode){
    if(i==PhyloChild[1+ntaxa(PT),2]-ntaxa(PT)) next
    
    p = PhyloParents[i+ntaxa(PT)]
    G3Seq = c(PhyloParents[p]-ntaxa(PT), p-ntaxa(PT), i)
    if(label[G3Seq[1]] == label[G3Seq[2]] && label[G3Seq[2]] == label[G3Seq[3]]) next
    Cntr = Cntr + 1
    Group3[Cntr,] = c(PhyloParents[p]-ntaxa(PT), p-ntaxa(PT), i)
  }  
  Group3[1:Cntr,]
}


InitInteraction2 = function(G3,label){
  Interaction2 = list()
  Cntr = 0
  for(i in 1:nrow(G3)){
    if(i>=1)
      for(j in 1:(i-1)){
        if(length(intersect(G3[i,], G3[j,]))==2){
          if(length(Interaction2)< i) Interaction2[[i]] = G3[j,]
          else Interaction2[[i]] = c(Interaction2[[i]], G3[j,])
        }
      }
    if(length(Interaction2)<i ) {
      Interaction2[[i]] = -1
      lseq1 = label[G3[i,]]
      if (length(unique(lseq1))==2) {
        attr(Interaction2[[i]],"type") = "112"
        dup = lseq1[duplicated(lseq1)]
        attr(Interaction2[[i]],"idxdist") = c(dup, setdiff(lseq1,dup))
      }
      else if (length(unique(lseq1))==3){
        attr(Interaction2[[i]],"type") = "123"
        attr(Interaction2[[i]],"idxdist") = lseq1
      }
      next
    }
    seq1 = G3[i,]; seq2 = Interaction2[[i]]; allseq = c(seq1, seq2)
    lseq1 = label[seq1]; lseq2 = label[seq2]; lallseq = label[allseq]
    if(length(seq2)==3){
      if(length(unique(lallseq))==4) {
        attr(Interaction2[[i]],"type") = "123124"
        dup = intersect(lseq1, lseq2)
        DistSeq = c(dup,setdiff(lseq1,dup),setdiff(lseq2,dup))
        attr(Interaction2[[i]],"idxdist") = DistSeq
      }else if(length(unique(lallseq))==2) {
        attr(Interaction2[[i]],"type") = "112122"
        dup = (lseq1)[duplicated(lseq1)]
        DistSeq = c(dup, setdiff(lseq1,dup))
        attr(Interaction2[[i]],"idxdist") = DistSeq
      }else if(length(unique(lallseq))==3){
        if(length(unique(lseq1))==2&&length(unique(lseq2))==2) {
          attr(Interaction2[[i]],"type") = "112113"
          dup = (lseq1)[duplicated(lseq1)]
          DistSeq = c(dup, setdiff(lseq1,dup), setdiff(lseq2,dup))
          attr(Interaction2[[i]],"idxdist") = DistSeq
        }
        if(length(unique(lseq1))==2&&length(unique(lseq2))==3) {
          attr(Interaction2[[i]],"type") = "112123"
          dup = (lseq1)[duplicated(lseq1)]
          DistSeq = c(dup, setdiff(lseq1,dup), setdiff(lseq2,lseq1))
          attr(Interaction2[[i]],"idxdist") = DistSeq
        }
        if(length(unique(lseq1))==3&&length(unique(lseq2))==2) {
          attr(Interaction2[[i]],"type") = "123112"
          dup = (lseq2)[duplicated(lseq2)]
          DistSeq = c(dup, setdiff(lseq2,dup),setdiff(lseq1, lseq2))
          attr(Interaction2[[i]],"idxdist") = DistSeq
        }
        
      }
    }else if (length(seq2)==6){
      if(length(unique(lseq1))==2&&length(unique(lseq2))==3) {
        attr(Interaction2[[i]],"type") = "112122123"
        dup = (lseq1)[duplicated(lseq1)]
        DistSeq = c(dup, setdiff(lseq1,dup), setdiff(lseq2,lseq1))
        attr(Interaction2[[i]],"idxdist") = DistSeq
      }else if(length(unique(lseq1))==2&&length(unique(lseq2))==4) {
        attr(Interaction2[[i]],"type") = "112123124"
        dup = (lseq1)[duplicated(lseq1)]
        DistSeq = c(dup, setdiff(lseq1,dup), setdiff(lseq2,lseq1))
        attr(Interaction2[[i]],"idxdist") = DistSeq
      }else if(length(unique(lseq1))==3&&length(unique(lseq2))==2) {
        attr(Interaction2[[i]],"type") = "123112122"
        dup = unique(lseq2)
        DistSeq = c(dup, setdiff(lseq1,dup))
        attr(Interaction2[[i]],"idxdist") = DistSeq
      }else if(length(unique(lseq1))==3&&length(unique(lseq2))==3) {
        attr(Interaction2[[i]],"type") = "123122124"
        tt = table(lallseq)
        d1 = as.numeric(names(tt)[which(tt==3)])
        d2 = as.numeric(names(tt)[which(tt==4)])
        DistSeq = c(d1, d2, setdiff(lseq1,c(d1,d2)),setdiff(lseq2,lseq1))
        attr(Interaction2[[i]],"idxdist") = DistSeq
      }
    } 
  }
  Interaction2
}

InitInteraction1 = function(G3, label){
  Interaction1 = list()
  Cntr = 0
  for(i in 1:nrow(G3)){
    iC = 0
    Interaction1[[i]] = list()
    for(j in 1:(i-1)){
      if(length(intersect(G3[i,], G3[j,]))==1){
        iC = iC + 1
        Interaction1[[i]][[iC]] = G3[j,]
      }
    }
    if(length(Interaction1[[i]])==0) next
    for(j in 1:length(Interaction1[[i]])){
      seq1 = G3[i,]; seq2 = Interaction1[[c(i,j)]]
      lseq1 = label[seq1]; lseq2 = label[seq2]
      if(length(unique(lseq1))==2 ){
        dup = (lseq1)[duplicated(lseq1)]
        ndup = setdiff(lseq1,dup)
        if(sum(lseq2==dup)==2){
          attr(Interaction1[[c(i,j)]],"type") = "e112113"
          DistSeq = c(dup, ndup ,setdiff(lseq2,dup))
          attr(Interaction1[[c(i,j)]],"idxdist") = DistSeq
        }else if(sum(lseq2==dup)==1 && length(unique(lseq2))==2){
          attr(Interaction1[[c(i,j)]],"type") = "e112133"
          DistSeq = c(dup, ndup,setdiff(lseq2,dup))
          attr(Interaction1[[c(i,j)]],"idxdist") = DistSeq
        }else if(sum(lseq2==dup)==1 && length(unique(lseq2))==3){
          attr(Interaction1[[c(i,j)]],"type") = "e112134"
          DistSeq = c(dup, ndup,setdiff(lseq2,dup))
          attr(Interaction1[[c(i,j)]],"idxdist") = DistSeq
        }else if(sum(lseq2==ndup)==2){
          attr(Interaction1[[c(i,j)]],"type") = "e112223"
          DistSeq = c(dup, ndup ,setdiff(lseq2,ndup))
          attr(Interaction1[[c(i,j)]],"idxdist") = DistSeq
        }else if(sum(lseq2==ndup)==1 && length(unique(lseq2))==2){
          attr(Interaction1[[c(i,j)]],"type") = "e112233"
          DistSeq = c(dup, ndup ,setdiff(lseq2,ndup))
          attr(Interaction1[[c(i,j)]],"idxdist") = DistSeq
        }else if(sum(lseq2==ndup)==1 && length(unique(lseq2))==3){
          attr(Interaction1[[c(i,j)]],"type") = "e112234"
          DistSeq = c(dup, ndup ,setdiff(lseq2,ndup))
          attr(Interaction1[[c(i,j)]],"idxdist") = DistSeq
        }
      }else if (length(unique(lseq1))==3){
        dup = intersect(lseq1, lseq2)
        if( sum(lseq2==dup)==2){
          attr(Interaction1[[c(i,j)]],"type") = "e123114"
          DistSeq = c(lseq1, setdiff(lseq2,dup))
          attr(Interaction1[[c(i,j)]],"idxdist") = DistSeq
        }else if(sum(lseq2==dup)==1 && length(setdiff(lseq2,dup))==1 ){
          attr(Interaction1[[c(i,j)]],"type") = "e123144"
          DistSeq = c(lseq1, setdiff(lseq2,dup))
          attr(Interaction1[[c(i,j)]],"idxdist") = DistSeq
        }else if(sum(lseq2==dup)==1 && length(setdiff(lseq2,dup))==2 ){
          attr(Interaction1[[c(i,j)]],"type") = "e123145"
          DistSeq = c(lseq1, setdiff(lseq2,dup))
          attr(Interaction1[[c(i,j)]],"idxdist") = DistSeq
        }
        
      }
      
    }
  }
  Interaction1 
}


InitInteraction0 = function(G3, label){
  Interaction0 = list()
  Cntr = 0
  for(i in 1:nrow(G3)){
    iC = 0
    Interaction0[[i]] = list()
    if(i==1) next
    for(j in 1:(i-1)){
      if(length(intersect(G3[i,], G3[j,]))==0 && length(intersect(label[G3[i,]],label[G3[j,]]))>0) {
        iC = iC + 1
        Interaction0[[i]][[iC]] = G3[j,]
      }
    }
  }
  Interaction0
}

InitIndep = function(G3,label){
  Indep = list()
  Cntr = 0
  for(i in 1:nrow(G3)){
    iC = 0
    Indep[[i]] = list()
    if(i==1) next
    for(j in 1:(i-1)){
      if(length(intersect(label[G3[i,]],label[G3[j,]]))==0) {
        iC = iC + 1
        Indep[[i]][[iC]] = G3[j,]
      }
    }
  }
  Indep
}



label = InitLabel()
labeldf = sapply(1:length(unique(label)), function(x)sum(label==x))

PT_G3 = PT
PT_G3$node.label = paste(1:PT$Nnode,label)
plot(PT_G3,show.node.label=TRUE,show.tip.label=FALSE)
#Color the nodes

library(geiger)
#PT_G3$node.label = paste(1:PT$Nnode)
PT_G3$node.label = label
par(mar=c(1,1,1,1))
plot(PT_G3,show.node.label=TRUE,show.tip.label=FALSE,use.edge.length = FALSE,cex=0.7)
RainColor = sample(rainbow(max(label)), replace=FALSE, size=max(label))
#nodelabels(pch=21, col="black", bg=RainColor[label], cex=1)



G3 = InitGroup3(label)
Interaction2 = InitInteraction2(G3, label)
Interaction1 = InitInteraction1(G3, label)
Interaction0 = InitInteraction0(G3, label)
Indep = InitIndep(G3, label)


temp = function(n){
  if(length(Interaction2[[n]])==1) return (0)
  else return (length(Interaction2[[n]])/3)
}
stopifnot(sapply(1:(nrow(G3)), function(n)length(Interaction0[[n]])+length(Interaction1[[n]])+temp(n)+length(Indep[[n]])) == 0:(nrow(G3)-1))


#For each C, need:
#f31(x), F31(x)
#f32(x), F32(x)
#f311(x), F311s(x) computed from f31
#f21(x), F21(x)
#f22(x), F22(x)
#f11(x)


f31 = function(x) pchisq(C-x,2)*dchisqM(x,1)/pchisqC3
f21 = function(x) pchisq(C-x,1)*dchisqM(x,1)/pchisqC2
f11 = function(x) (x<C)*dchisqM(x,1)/pchisqC1
f32 = function(x) pchisq(C-sum(x),1)*prod(dchisqM(x,1))/pchisqC3
f32s = function(x) pchisq(C-x,1)*dchisqM(x,2)/pchisqC3
f22 = function(x) (sum(x)<C)*prod(dchisqM(x,1))/pchisqC2
f22s = function(x) (x<C)*dchisqM(x,2)/pchisqC2
f33 = function(x) (sum(x)<C)*prod(dchisqM(x,1))/pchisqC3
f3132s = function(x) (sum(x)<C)*dchisqM(x[1],1)*dchisqM(x[2],2)/pchisqC3

f_1 = function(x,d) switch(d, f11(x),f21(x),f31(x))
f_2 = function(x,d) switch(d-1, f22(x), f32(x))
f_2s = function(x,d) switch(d-1, f22s(x), f32s(x))

F_1_1s = function(x,d1,d2){
  if(x>=2*C-2*1e-3) return (1)
  if(x<1e-3) return (0)
  dd = paste(d1,d2, sep="")
  Fg = switch(dd, "33"=F3131s_g, "32"=F3121s_g, "23"=F3121s_g, "22"=F2121s_g, "21"=F2111s_g, "12"=F2111s_g, "31"=F3111s_g, "13"=F3111s_g, "11"=F1111s_g)
  Fg[round(x/1e-3)]
}
F_2s = function(x,d){
  if(x>=C-1e-3) return (1)
  if(x<1e-3) return (0)
  switch(d-1, pchisq(x,2)/pchisq(C,2),F32s_g[round(x/1e-3)])
}
F_1 = function(x,d){
  if(x>=C-1e-3) return (1)
  if(x<1e-3) return (0)
  switch(d, pchisq(x,1)/pchisq(C,1), F21_g[round(x/1e-3)], F31_g[round(x/1e-3)])
}

f11_a = function(x,a) (x<a)*dchisqM(x,1)/pchisq(a,1)
f21_a = function(x,a) pchisq(a-x,1)*dchisqM(x,1)/pchisq(a,2)
f31_a = function(x,a) pchisq(a-x,2)*dchisqM(x,1)/pchisq(a,3)
f_1_a = function(x,a,df) switch(df, f11_a(x,a), f21_a(x,a), f31_a(x,a))
F_1_a = function(x,a,df){
  if(x<=0) return (0) else if (x>a) return (1)
  ff = switch(df, f11_a, f21_a, f31_a)
  integrate(function(y)ff(y,a),0,x)$value
}

                        


CalcCDF = function(){
  IntGrid = seq(1e-3,C,by=1e-3)
  IntGrids2 = seq(1e-3,2*C,by=1e-3)
  
  pchisqC3 <<- pchisq(C,3)
  pchisqC2 <<- pchisq(C,2)
  pchisqC1 <<- pchisq(C,1)
  
  F31_g <<- unlist(mclapply(IntGrid,function(upper)integrate(f31,0,upper,subdivisions=1000)$value, mc.cores=4))
  F21_g <<- unlist(mclapply(IntGrid,function(upper)integrate(f21,0,upper,subdivisions=1000)$value, mc.cores=4))
  F32s_g <<- unlist(mclapply(IntGrid,function(upper)integrate(f32s,0,upper,subdivisions=1000)$value, mc.cores=4))
  
   
  #Mixture of 2d and cuhre
  i = 0
  F3131s_g <<- numeric(length(IntGrids2))
  F3121s_g <<- numeric(length(IntGrids2))
  F3111s_g <<- numeric(length(IntGrids2))
  F2121s_g <<- numeric(length(IntGrids2))
  F2111s_g <<- numeric(length(IntGrids2))
  F1111s_g <<- numeric(length(IntGrids2))
  while(1){
    i = i + 1
    s1 = i*100-99; s2 = i*100 
    F3131s_g[s1:s2] <<- unlist(mclapply(IntGrids2[s1:s2],function(u)cuhre(2,1,function(x)f31(x[1])*f31(x[2])*(sum(x)<u),lower=rep(0,2),upper=rep(u,2),flags =list(verbose=0),max.eval=5e7)$value, mc.cores=4))
    F3121s_g[s1:s2] <<- unlist(mclapply(IntGrids2[s1:s2],function(u)cuhre(2,1,function(x)f31(x[1])*f21(x[2])*(sum(x)<u),lower=rep(0,2),upper=rep(u,2),flags =list(verbose=0),max.eval=5e7)$value, mc.cores=4))
    F3111s_g[s1:s2] <<- unlist(mclapply(IntGrids2[s1:s2],function(u)cuhre(2,1,function(x)f31(x[1])*f11(x[2])*(sum(x)<u),lower=rep(0,2),upper=rep(u,2),flags =list(verbose=0),max.eval=5e7)$value, mc.cores=4))
    F2121s_g[s1:s2] <<- unlist(mclapply(IntGrids2[s1:s2],function(u)cuhre(2,1,function(x)f21(x[1])*f21(x[2])*(sum(x)<u),lower=rep(0,2),upper=rep(u,2),flags =list(verbose=0),max.eval=5e7)$value, mc.cores=4))
    F2111s_g[s1:s2] <<- unlist(mclapply(IntGrids2[s1:s2],function(u)cuhre(2,1,function(x)f21(x[1])*f11(x[2])*(sum(x)<u),lower=rep(0,2),upper=rep(u,2),flags =list(verbose=0),max.eval=5e7)$value, mc.cores=4))
    F1111s_g[s1:s2] <<- unlist(mclapply(IntGrids2[s1:s2],function(u)cuhre(2,1,function(x)f11(x[1])*f11(x[2])*(sum(x)<u),lower=rep(0,2),upper=rep(u,2),flags =list(verbose=0),max.eval=5e7)$value, mc.cores=4))
    cat(i,s1,s2,"\n")
    if(min(F3131s_g[s2],F3121s_g[s2],F3111s_g[s2],F2121s_g[s2],F2111s_g[s2],F1111s_g[s2])>0.1) break
  }
  ns1 = s2 + 1; ns2 = length(IntGrids2)
  F3131s_g[ns1:ns2] <<- unlist(mclapply(IntGrids2[ns1:ns2],function(u)cuhre(1,1,function(x)f31(x)*F_1(u-x,3),lower=0,upper=u,flags =list(verbose=0), max.eval=5e7)$value, mc.cores=4))
  F3121s_g[ns1:ns2] <<- unlist(mclapply(IntGrids2[ns1:ns2],function(u)cuhre(1,1,function(x)f31(x)*F_1(u-x,2),lower=0,upper=u,flags =list(verbose=0), max.eval=5e7)$value, mc.cores=4))
  F3111s_g[ns1:ns2] <<- unlist(mclapply(IntGrids2[ns1:ns2],function(u)cuhre(1,1,function(x)f31(x)*F_1(u-x,1),lower=0,upper=u,flags =list(verbose=0), max.eval=5e7)$value, mc.cores=4))
  F2121s_g[ns1:ns2] <<- unlist(mclapply(IntGrids2[ns1:ns2],function(u)cuhre(1,1,function(x)f21(x)*F_1(u-x,2),lower=0,upper=u,flags =list(verbose=0), max.eval=5e7)$value, mc.cores=4))
  F2111s_g[ns1:ns2] <<- unlist(mclapply(IntGrids2[ns1:ns2],function(u)cuhre(1,1,function(x)f21(x)*F_1(u-x,1),lower=0,upper=u,flags =list(verbose=0), max.eval=5e7)$value, mc.cores=4))
  F1111s_g[ns1:ns2] <<- unlist(mclapply(IntGrids2[ns1:ns2],function(u)cuhre(1,1,function(x)f11(x)*F_1(u-x,1),lower=0,upper=u,flags =list(verbose=0), max.eval=5e7)$value, mc.cores=4))
}





P_112 = function(df){
  IntegrandF = function(x) f_1(x,df[2])*(1-F_2s(C-x,df[1]))
  cuhre(ndim=1,ncomp=1,integrand=IntegrandF,lower=0,upper=C,max.eval=5e7,flags=list(verbose=0))$value
}
P_123 = function(df){
  IntegrandF = function(x) f_1(x,df[1])*(1-F_1_1s(C-x,df[2],df[3]))
  cuhre(ndim=1,ncomp=1,integrand=IntegrandF,lower=0,upper=C,max.eval=5e7,flags=list(verbose=0))$value
}
P_112122 = function(df){
  IntegrandFa = function(x) f_1(x[1],df[1])*f_1(x[2],df[2])*(1-F_1_a(C-sum(x),C-x[1],df[1]-1))*F_1_a(C-sum(x),C-x[2],df[2]-1)
  cuhre(ndim=2,ncomp=1,integrand=IntegrandFa,lower=rep(0,2),upper=rep(C,2),max.eval=5e7, flags=list(verbose=0))$value
}
P_112113 = function(df){
  IntegrandF = function(x) f_2(x,df[1])*(1-F_1(C-sum(x),df[2]))*F_1(C-sum(x),df[3])
  cuhre(ndim=2,ncomp=1,integrand=IntegrandF,lower=rep(0,2),upper=rep(C,2),max.eval=5e7,flags=list(verbose=0))$value
}
P_112123 = function(df){
  IntegrandFa = function(x) f_1(x[1],df[1])*f_1(x[2],df[2])*(1-F_1_a(C-sum(x),C-x[1],df[1]-1))*F_1(C-sum(x),df[3])
  cuhre(ndim=2,ncomp=1,integrand=IntegrandFa,lower=rep(0,2),upper=rep(C,2),max.eval=5e7,flags=list(verbose=0))$value
}
P_123112 = function(df){
  IntegrandFa = function(x) f_1(x[1],df[1])*f_1(x[2],df[2])*F_1_a(C-sum(x),C-x[1],df[1]-1)*(1-F_1(C-sum(x),df[3]))
  cuhre(ndim=2,ncomp=1,integrand=IntegrandFa,lower=rep(0,2),upper=rep(C,2),max.eval=5e7,flags=list(verbose=0))$value
}
P_123124 = function(df){
  IntegrandF = function(x) f_1(x[1],df[1])*f_1(x[2],df[2])*(1-F_1(C-sum(x),df[3]))*F_1(C-sum(x),df[4])
  cuhre(ndim=2,ncomp=1,integrand=IntegrandF,lower=rep(0,2),upper=rep(C,2),max.eval=5e7,flags=list(verbose=0))$value
}
P_112122123 = function(df){
  IntegrandFa = function(x) f_1(x[1],df[1])*f_1(x[2],df[2])*(1-F_1_a(C-sum(x),C-x[1],df[1]-1))*F_1_a(C-sum(x),C-x[2],df[2]-1)*F_1(C-sum(x),df[3])
  cuhre(ndim=2,ncomp=1,integrand=IntegrandFa,lower=rep(0,2),upper=rep(C,2),max.eval=5e7,flags=list(verbose=0))$value
}
P_112123124 = function(df){
  IntegrandFa = function(x) f_1(x[1],df[1])*f_1(x[2],df[2])*(1-F_1_a(C-sum(x),C-x[1],df[1]-1))*F_1(C-sum(x),df[3])*F_1(C-sum(x),df[4])
  cuhre(ndim=2,ncomp=1,integrand=IntegrandFa,lower=rep(0,2),upper=rep(C,2),max.eval=5e7,flags=list(verbose=0))$value
}
P_123112122 = function(df){
  IntegrandFa = function(x) f_1(x[1],df[1])*f_1(x[2],df[2])*(1-F_1(C-sum(x),df[3]))*F_1_a(C-sum(x),C-x[1],df[1]-1)*F_1_a(C-sum(x),C-x[2],df[2]-1)
  cuhre(ndim=2,ncomp=1,integrand=IntegrandFa,lower=rep(0,2),upper=rep(C,2),max.eval=5e7,flags=list(verbose=0))$value
}
P_123122124 = function(df){
  IntegrandFa = function(x) f_1(x[1],df[1])*f_1(x[2],df[2])*(1-F_1(C-sum(x),df[3]))*F_1_a(C-sum(x),C-x[2],df[2]-1)*F_1(C-sum(x),df[4])
  cuhre(ndim=2,ncomp=1,integrand=IntegrandFa,lower=rep(0,2),upper=rep(C,2),max.eval=5e7,flags=list(verbose=0))$value
}


#Error Probabilities
Pe_112113 = function(df){
  IntegrandF = function(x) f33(x)*(1-F_1(C-x[1]-x[2],df[2]))*(1-F_1(C-x[2]-x[3],df[3]))
  suave(ndim=3,ncomp=1,integrand=IntegrandF,lower=rep(0,3),upper=rep(C,3),max.eval=5e7,flags=list(verbose=0,final=1,pseudo.random=0,mersenne.seed=NULL))$value
  #cuhre(ndim=3,ncomp=1,integrand=IntegrandF,lower=rep(0,3),upper=rep(C,3),max.eval=5e7)$value
}
Pe_112133 = function(df){
  IntegrandF = function(x) f_2(x,df[1])*(1-F_1(C-sum(x),df[2]))*(1-F_2s(C-x[2],df[3]))
  cuhre(ndim=2,ncomp=1,integrand=IntegrandF,lower=rep(0,2),upper=rep(C,2),max.eval=5e7,flags=list(verbose=0))$value
}
Pe_112134 = function(df){
  IntegrandF = function(x) f_2(x,df[1])*(1-F_1(C-sum(x),df[2]))*(1-F_1_1s(C-x[2],df[3],df[4]))
  cuhre(ndim=2,ncomp=1,integrand=IntegrandF,lower=rep(0,2),upper=rep(C,2),max.eval=5e7,flags=list(verbose=0))$value
}
Pe_112223 = function(df){
  IntegrandF = function(x) f_2(x,df[2])*(1-F_2s(C-x[1],df[1]))*(1-F_1(C-sum(x),df[3]))
  cuhre(ndim=2,ncomp=1,integrand=IntegrandF,lower=rep(0,2),upper=rep(C,2),max.eval=5e7,flags=list(verbose=0))$value
}
Pe_112233 = function(df){
  IntegrandF = function(x) f_1(x,df[2])*(1-F_2s(C-x,df[1]))*(1-F_2s(C-x,df[3]))
  cuhre(ndim=1,ncomp=1,integrand=IntegrandF,lower=0,upper=C,max.eval=5e7,flags=list(verbose=0))$value
}
Pe_112234 = function(df){
  IntegrandF = function(x) f_1(x,df[2])*(1-F_2s(C-x,df[1]))*(1-F_1_1s(C-x,df[3],df[4]))
  cuhre(ndim=1,ncomp=1,integrand=IntegrandF,lower=0,upper=C,max.eval=5e7,flags=list(verbose=0))$value
}
Pe_123114 = function(df) {
  IntegrandF = function(x) f_2(x,df[1])*(1-F_1_1s(C-x[1],df[2],df[3]))*(1-F_1(C-sum(x),df[4]))
  cuhre(ndim=2,ncomp=1,integrand=IntegrandF,lower=rep(0,2),upper=rep(C,2),max.eval=5e7,flags=list(verbose=0))$value
}
Pe_123144 = function(df){
  IntegrandF = function(x) f_1(x,df[1])*(1-F_1_1s(C-x,df[2],df[3]))*(1-F_2s(C-x,df[4]))
  cuhre(ndim=1,ncomp=1,integrand=IntegrandF,lower=0,upper=C,max.eval=5e7,flags=list(verbose=0))$value
}
Pe_123145 = function(df){
  IntegrandF = function(x) f_1(x,df[1])*(1-F_1_1s(C-x,df[2],df[3]))*(1-F_1_1s(C-x,df[4],df[5]))
  cuhre(ndim=1,ncomp=1,integrand=IntegrandF,lower=0,upper=C,max.eval=5e7,flags=list(verbose=0))$value
}

#Generic s Proabilities Function
Ps = function(l1, l2, labeldf){
  cl = intersect(l1, l2)
  ll1 = l1[-which(l1==cl)]; ll2 = l2[-which(l2==cl)]
  
  componentf = function(x, ll){
    if(length(ll)==1) 1-F_1(C-x,labeldf[ll])
    else if(length(ll)==2 && length(unique(ll))==2) 1-F_1_1s(C-x,labeldf[ll[1]],labeldf[ll[2]])
    else if(length(ll)==2 && length(unique(ll))==1) 1-F_2s(C-x,labeldf[ll[1]])
  }
  
  if( sum(l1==cl)+sum(l2==cl) == 3){
    IntegrandF = function(x) f3132s(x)*componentf(x[1],ll1)*componentf(x[2],ll2)
  }else if (sum(l1==cl)+sum(l2==cl) == 2){ 
    IntegrandF = function(x) f_2(x,labeldf[cl])*componentf(x[1],ll1)*componentf(x[2],ll2)}
  
  cuhre(ndim=2,ncomp=1,integrand=IntegrandF,lower=rep(0,2),upper=rep(C,2),max.eval=5e7,flags=list(verbose=0))$value
}

FindPfunction = function(type){
  switch(type, "112"=P_112,"123"=P_123, "112122"=P_112122,"112113"=P_112113,"112123"=P_112123,"123112"=P_123112,"123124"=P_123124,"112122123"=P_112122123,"112123124"=P_112123124,"123112122"=P_123112122,"123122124"=P_123122124)
}
FindPefunction = function(type){
  switch(type, "e112113"=Pe_112113, "e112133"=Pe_112133,"e112134"=Pe_112134,"e112223"=Pe_112223,"e112233"=Pe_112233,"e112234"=Pe_112234,"e123114"=Pe_123114,"e123144"=Pe_123144,"e123145"=Pe_123145)
}
IndepProb = function(l, labeldf){
  if(length(unique(l))==3) return(P_123(labeldf[l]))
  dup = l[duplicated(l)]
  cl = c(dup, setdiff(l,dup))
  P_112(labeldf[cl])
}

##Wrap main program in a function:
TDMpvalues = function(CC){
  C <<- CC
  CalcCDF()
  GridCProb = prod(unlist(lapply(labeldf,function(x)pchisq(C,x))))
  message("Stage I compelte")
  
  UpperPvalue = 1-GridCProb
  for(i in 1:length(Interaction2)){
    P_func = FindPfunction(attributes(Interaction2[[i]])$type)
    CurrentProb = P_func(labeldf[attributes(Interaction2[[i]])$idxdist] )
    UpperPvalue = UpperPvalue + GridCProb*CurrentProb
    #cat(i, attributes(Interaction2[[i]])$type, CurrentProb, UpperPvalue, "\n")
  }
  message("Stage II complete")
  PuBound_Int1 = 0
  for(i in 1:nrow(G3)){
    if(length(Interaction1[[i]])==0) next
    for(j in Interaction1[[i]]) {
      
      Pe_func = FindPefunction(attributes(j)$type)
      CurrentProb = Pe_func(labeldf[attributes(j)$idxdist])
      PuBound_Int1 = PuBound_Int1 + GridCProb*CurrentProb
    }
    #cat(i, PuBound_Int1, "\n")
  }
  message("Stage III complete")
  PuBound_int0 = 0
  for(i in 1:nrow(G3)){
    if(length(Interaction0[[i]])==0) next
    for(j in Interaction0[[i]]) PuBound_int0 = PuBound_int0 + GridCProb*Ps(label[G3[i,]], label[j], labeldf)
    #cat(i, PuBound_int0, length(Interaction0[[i]]), "\n")
    cat(i," ")
  }
  message("Stage IV complete")
  PuBound_indep = 0
  for(i in 1:nrow(G3)){
    if(length(Indep[[i]])==0) next
    for(j in Indep[[i]])  PuBound_indep = PuBound_indep + GridCProb*IndepProb(label[G3[i,]],labeldf)*IndepProb(label[j],labeldf)
    #cat(i, PuBound_indep, length(Indep[[i]]),"\n")
  }
  
  LowerPvalue = UpperPvalue - (PuBound_Int1 + PuBound_int0 + PuBound_indep)
  cat(UpperPvalue, PuBound_Int1 + PuBound_int0 + PuBound_indep,"\n")
}
##This is as far as all the functions go






##Begin main program
C = max(TripletStat)
CalcCDF()
GridCProb = prod(unlist(lapply(labeldf,function(x)pchisq(C,x))))

UpperPvalue = 1-GridCProb
for(i in 1:length(Interaction2)){
  P_func = FindPfunction(attributes(Interaction2[[i]])$type)
  CurrentProb = P_func(labeldf[attributes(Interaction2[[i]])$idxdist] )
  UpperPvalue = UpperPvalue + GridCProb*CurrentProb
  cat(i, attributes(Interaction2[[i]])$type, CurrentProb, UpperPvalue, "\n")
}

#cat(1-GridCProb, UpperPvalue - (1-GridCProb))

PuBound_Int1 = 0
for(i in 1:nrow(G3)){
  if(length(Interaction1[[i]])==0) next
  for(j in Interaction1[[i]]) {
    
    Pe_func = FindPefunction(attributes(j)$type)
    CurrentProb = Pe_func(labeldf[attributes(j)$idxdist])
    PuBound_Int1 = PuBound_Int1 + GridCProb*CurrentProb
  }
  cat(i, PuBound_Int1, "\n")
}

PuBound_int0 = 0
for(i in 1:nrow(G3)){
  if(length(Interaction0[[i]])==0) next
  for(j in Interaction0[[i]]) PuBound_int0 = PuBound_int0 + GridCProb*Ps(label[G3[i,]], label[j], labeldf)
  cat(i, PuBound_int0, length(Interaction0[[i]]), "\n")
}

PuBound_indep = 0
for(i in 1:nrow(G3)){
  if(length(Indep[[i]])==0) next
  for(j in Indep[[i]])  PuBound_indep = PuBound_indep + GridCProb*IndepProb(label[G3[i,]],labeldf)*IndepProb(label[j],labeldf)
  cat(i, PuBound_indep, length(Indep[[i]]),"\n")
}

LowerPvalue = UpperPvalue - (PuBound_Int1 + PuBound_int0 + PuBound_indep)
cat(UpperPvalue, PuBound_Int1 + PuBound_int0 + PuBound_indep)



##End of main program

















##Simiulation verification
rm(AG, AG_OTU, OTU1, OTU2, OTUtable1, OTUtable2, OTSample1, OTSample2); gc()

AllCluster = matrix(0,PT$Nnode-2,3)
Cntr = 0
for(i in 3:PT$Nnode){
  if(i==PhyloChild[1+ntaxa(PT),2]-ntaxa(PT)) next
  p = PhyloParents[i+ntaxa(PT)]
  Cntr = Cntr + 1
  AllCluster[Cntr,] = c(PhyloParents[p]-ntaxa(PT), p-ntaxa(PT), i)
}  
AllCluster = AllCluster[1:Cntr,]

FindMaxZ = function(x) {
  z = apply(AllCluster, 1, function(row) sum(x[row]))
  max(z)
}

MCn = 5e4
GenChiSqM = function(){matrix(rchisq(ntaxa(PT)*MCn,1),MCn)}



C = 25; EmpiricalPvalue_C25 = unlist( mclapply(1:500, function(nn) sum(apply(GenChiSqM(), 1, FindMaxZ)>C)/MCn, mc.cores = 4))
C = 20; EmpiricalPvalue_C20 = unlist( mclapply(1:500, function(nn) sum(apply(GenChiSqM(), 1, FindMaxZ)>C)/MCn, mc.cores = 4))
C = 15; EmpiricalPvalue_C15 = unlist( mclapply(1:500, function(nn) sum(apply(GenChiSqM(), 1, FindMaxZ)>C)/MCn, mc.cores = 4))

#hist(EmpiricalPvalue, xlim = c(LowerPvalue-0.0001, UpperPvalue+0.0001) ,breaks=10, main="Empirical p-value, n = 60, C = 17")
library(mosaic)
UpperPvalue = 0.1003442
LowerPvalue = UpperPvalue - 0.006534947
EP = EmpiricalPvalue_C15
nkOTU = as.character(ntaxa(PT))

histogram(EP, breaks=seq(min(EP),max(EP), length.out=21), main="",xlab = "P-value",col="gray")
ladd(panel.abline(v = UpperPvalue, col="red"))
ladd(panel.abline(v = LowerPvalue, col="red"))









#####Test part######
GenChiSq = function(df){ChiSqSample = matrix(rchisq(length(df)*3*5e5,1),5e5)}
ChiGridCheck = function(x,df) {
  flag = TRUE
  for(i in 1:length(df)) flag = flag * (sum(x[(i*3-2):(i*3-3+df[i])])<C)
  flag
}


Ps(c(5,18,19),c(5,5,3),labeldf)
df = c(3,2,3,1)
sum(apply(GenChiSq(df),1,function(x)ChiGridCheck(x,df)*(sum(x[c(1,4,10)])>C)*(sum(x[c(2,3,7)])>C)  ))/5e5/ prod(sapply(df,function(d)pchisq(C,d)))

df = c(1,2,1,3,1)
Pe_123145(df)
sum(apply(GenChiSq(df),1,function(x)ChiGridCheck(x,df)*(sum(x[c(1,4,7)])>C)*(sum(x[c(1,10,13)])>C)  ))/5e5/ prod(sapply(df,function(d)pchisq(C,d)))

df = c(3,3,2,3)
Pe_123144(df)
sum(apply(GenChiSq(df),1,function(x)ChiGridCheck(x,df)*(sum(x[c(1,4,7)])>C)*(sum(x[c(1,10,11)])>C)  ))/5e5/ prod(sapply(df,function(d)pchisq(C,d)))

df = c(3,3,2,3)
Pe_123114(df)
sum(apply(GenChiSq(df),1,function(x)ChiGridCheck(x,df)*(sum(x[c(1,4,7)])>C)*(sum(x[c(1,2,10)])>C)  ))/5e5/ prod(sapply(df,function(d)pchisq(C,d)))

df = c(3,3,2)
Pe_112113(df)
sum(apply(GenChiSq(df),1,function(x)ChiGridCheck(x,df)*(sum(x[c(1,2,4)])>C)*(sum(x[c(2,3,7)])>C)  ))/5e5/ prod(sapply(df,function(d)pchisq(C,d)))

df = c(3,3,2)
Pe_112133(df)
sum(apply(GenChiSq(df),1,function(x)ChiGridCheck(x,df)*(sum(x[c(1,2,4)])>C)*(sum(x[c(2,7,8)])>C)  ))/5e5/ prod(sapply(df,function(d)pchisq(C,d)))

df = c(3,1,2,3)
Pe_112134(df)
sum(apply(GenChiSq(df),1,function(x)ChiGridCheck(x,df)*(sum(x[c(1,2,4)])>C)*(sum(x[c(2,7,10)])>C)  ))/5e5/ prod(sapply(df,function(d)pchisq(C,d)))

df = c(2,3,2)
Pe_112223(df)
sum(apply(GenChiSq(df),1,function(x)ChiGridCheck(x,df)*(sum(x[c(1,2,4)])>C)*(sum(x[c(4,5,7)])>C)  ))/5e5/ prod(sapply(df,function(d)pchisq(C,d)))

df = c(3,3,3)
Pe_112233(df)
sum(apply(GenChiSq(df),1,function(x)ChiGridCheck(x,df)*(sum(x[c(1,2,4)])>C)*(sum(x[c(4,7,8)])>C)  ))/5e5/ prod(sapply(df,function(d)pchisq(C,d)))

df = c(3,3,2,1)
Pe_112234(df)
sum(apply(GenChiSq(df),1,function(x)ChiGridCheck(x,df)*(sum(x[c(1,2,4)])>C)*(sum(x[c(4,7,10)])>C)  ))/5e5/ prod(sapply(df,function(d)pchisq(C,d)))



ChiSqSample = matrix(rchisq(3*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[1:3])>C)*all(x<C)))/1e6/pchisq(C,1)^3

P_112(c(3,3))
1; ChiSqSample = matrix(rchisq(6*1e5,1),1e5)
sum(apply(ChiSqSample,1, function(x) (sum(x[1:3])<C)*(sum(x[4:6])<C)*(sum(x[2:4])>C) ))/1e5/pchisq(C,3)/pchisq(C,3)

P_112(c(3,2))
1; ChiSqSample = matrix(rchisq(6*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[1:3])<C)*(sum(x[4:5])<C)*(sum(x[2:4])>C) ))/1e6/pchisq(C,3)/pchisq(C,2)

P_123(c(1,3,2))
1; ChiSqSample = matrix(rchisq(9*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[1])<C)*(sum(x[4:6])<C)*(sum(x[7:8])<C)*(sum(x[c(1,4,7)])>C) ))/1e6/pchisq(C,3)/pchisq(C,2)/pchisq(C,1)


P_123122124(c(3,3,3,3))
1; ChiSqSample = matrix(rchisq(12*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[1:3])<C)*(sum(x[4:6])<C)*(sum(x[7:9])<C)*(sum(x[10:12])<C)*(sum(x[c(3:5)])<C)*(sum(x[c(3,4,7)])>C)*(sum(x[c(3,4,10)])<C)))/1e6/pchisq(C,3)/pchisq(C,3)/pchisq(C,3)/pchisq(C,3)

P_123122124(c(2,2,3,1))
1; ChiSqSample = matrix(rchisq(12*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[2:3])<C)*(sum(x[4:5])<C)*(sum(x[7:9])<C)*(sum(x[10:10])<C)*(sum(x[c(3:5)])<C)*(sum(x[c(3,4,7)])>C)*(sum(x[c(3,4,10)])<C)))/1e6/pchisq(C,3)/pchisq(C,2)/pchisq(C,2)/pchisq(C,1)


P_123112122(c(3,2,2))
1; ChiSqSample = matrix(rchisq(9*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[1:3])<C)*(sum(x[4:5])<C)*(sum(x[7:8])<C)*(sum(x[c(2:4)])<C)*(sum(x[c(3,4,7)])>C)*(sum(x[c(3,4,5)])<C)))/1e6/pchisq(C,3)/pchisq(C,2)/pchisq(C,2)


P_112123124(c(2,2,3,1))
1; ChiSqSample = matrix(rchisq(12*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[2:3])<C)*(sum(x[4:5])<C)*(sum(x[7:9])<C)*(sum(x[10:10])<C)*(sum(x[c(2:4)])>C)*(sum(x[c(3,4,7)])<C)*(sum(x[c(3,4,10)])<C)))/1e6/pchisq(C,3)/pchisq(C,2)/pchisq(C,2)/pchisq(C,1)

P_112122123(c(3,2,2))
1; ChiSqSample = matrix(rchisq(9*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[1:3])<C)*(sum(x[4:5])<C)*(sum(x[7:8])<C)*(sum(x[c(2:4)])>C)*(sum(x[c(3,4,7)])<C)*(sum(x[c(3,4,5)])<C)))/1e6/pchisq(C,3)/pchisq(C,2)/pchisq(C,2)


P_123124(c(3,3,2,3))
ChiSqSample = matrix(rchisq(12*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[1:3])<C)*(sum(x[4:6])<C)*(sum(x[7:8])<C)*(sum(x[10:12])<C)*(sum(x[c(2,4,7)])>C)*(sum(x[c(2,4,10)])<C)))/1e6/pchisq(C,3)/pchisq(C,3)/pchisq(C,3)/pchisq(C,2)


P_123112(c(3,2,3))
ChiSqSample = matrix(rchisq(9*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[1:3])<C)*(sum(x[4:5])<C)*(sum(x[7:9])<C)*(sum(x[2:4])<C)*(sum(x[c(2,4,7)])>C)))/1e6/pchisq(C,3)/pchisq(C,2)/pchisq(C,3)

P_112123(c(3,2,2))
ChiSqSample = matrix(rchisq(9*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[1:3])<C)*(sum(x[4:5])<C)*(sum(x[7:8])<C)*(sum(x[2:4])>C)*(sum(x[c(2,4,7)])<C)))/1e6/pchisq(C,3)/pchisq(C,2)/pchisq(C,2)


P_112122(c(3,2))
ChiSqSample = matrix(rchisq(6*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[1:3])<C)*(sum(x[4:5])<C)*(sum(x[2:4])>C)*(sum(x[3:5])<C)))/1e6/pchisq(C,3)/pchisq(C,2)

P_112122(c(2,3))
ChiSqSample = matrix(rchisq(6*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[1:3])<C)*(sum(x[4:5])<C)*(sum(x[2:4])<C)*(sum(x[3:5])>C)))/1e6/pchisq(C,3)/pchisq(C,2)

P_112122(c(2,2))
ChiSqSample = matrix(rchisq(4*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[1:2])<C)*(sum(x[3:4])<C)*(sum(x[1:3])<C)*(sum(x[2:4])>C)))/1e6/pchisq(C,2)/pchisq(C,2)

P_112122(c(3,3))
ChiSqSample = matrix(rchisq(6*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[1:3])<C)*(sum(x[4:6])<C)*(sum(x[2:4])>C)*(sum(x[3:5])<C)))/1e6/pchisq(C,3)/pchisq(C,3)


P_112113(c(3,3,3))
ChiSqSample = matrix(rchisq(9*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[1:3])<C)*(sum(x[4:6])<C)*(sum(x[7:9])<C)*(sum(x[2:4])>C)*(sum(x[c(2,3,7)])<C)))/1e6/pchisq(C,3)/pchisq(C,3)/pchisq(C,3)

P_112113(c(3,2,2))
ChiSqSample = matrix(rchisq(9*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[1:3])<C)*(sum(x[4:5])<C)*(sum(x[7:8])<C)*(sum(x[2:4])<C)*(sum(x[c(2,3,7)])>C)))/1e6/pchisq(C,3)/pchisq(C,2)/pchisq(C,2)

P_112113(c(3,1,1))
1; ChiSqSample = matrix(rchisq(9*1e6,1),1e6)
sum(apply(ChiSqSample,1, function(x) (sum(x[1:3])<C)*(sum(x[4:5])<C)*(sum(x[7])<C)*(sum(x[2:4])<C)*(sum(x[c(2,3,7)])>C)))/1e6/pchisq(C,3)/pchisq(C,2)/pchisq(C,2)

