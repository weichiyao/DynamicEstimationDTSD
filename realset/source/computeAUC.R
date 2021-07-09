computeAUC=function(dat, horizon3 = c(0, 1, 2))
{
  
  # library(AUC)
  library(pROC)
  
  # dat=read.table(file="d:/11000347/Desktop/landmarking_dynamic_prediction/discrete_time/bankruptcy_hoora/Allresultsforr.txt",header=TRUE)   # portable
  out=rbind(rep(0,10),rep(0,10))
  out=data.frame(out)
  # print(out)
  names(out)=c("year","t","u",
               "hbench","hsuppDTPO","hsupp0hlgr","hsephlgr","hSKHThlgr","hsupphlgr",
               "horizon3")
  
  if(horizon3 %in% c(1, 2)) {
    for(year in 1:3) {
      for(t in unique(dat[, "t"])) {
        dat1 = dat[dat[,"year"] == year & dat[,"t"]==t & dat[,"u"]==t + 3, ]
        
        # print(dat1)
        
        auc1=auc(dat1[,"y"],dat1[,"hbench"])
        auc2=auc(dat1[,"y"],dat1[,"hsuppDTPO"])
        auc3=auc(dat1[,"y"],dat1[,"hsupp0hlgr"])
        auc4=auc(dat1[,"y"],dat1[,"hsephlgr"])
        auc5=auc(dat1[,"y"],dat1[,"hSKHThlgr"])
        auc6=auc(dat1[,"y"],dat1[,"hsupphlgr"])
        
        
        out=rbind(out,c(year,t,t+3,auc1,auc2,auc3,auc4,auc5,auc6,horizon3))
        
        
      }}
  }
  
  
  if (horizon3 %in% c(0, 2)) {
    for (year in 1:3) {
      for (t in 0:3) {
        if (t == 0) {
          for(u in 4:6) {
            dat1=dat[dat[,"year"]==year & dat[,"t"]==t & dat[,"u"]==u,]
            
            auc1=auc(dat1[,"y"],dat1[,"hbench"])
            auc2=auc(dat1[,"y"],dat1[,"hsuppDTPO"])
            auc3=auc(dat1[,"y"],dat1[,"hsupp0hlgr"])
            auc4=auc(dat1[,"y"],dat1[,"hsephlgr"])
            auc5=auc(dat1[,"y"],dat1[,"hSKHThlgr"])
            auc6=auc(dat1[,"y"],dat1[,"hsupphlgr"])
            
            out=rbind(out,c(year,t,u,auc1,auc2,auc3,auc4,auc5,auc6,horizon3))
          }
        }
        if (t == 1) {
          for(u in 5:6) {
            dat1=dat[dat[,"year"]==year & dat[,"t"]==t & dat[,"u"]==u,]
            
            auc1=auc(dat1[,"y"],dat1[,"hbench"])
            auc2=auc(dat1[,"y"],dat1[,"hsuppDTPO"])
            auc3=auc(dat1[,"y"],dat1[,"hsupp0hlgr"])
            auc4=auc(dat1[,"y"],dat1[,"hsephlgr"])
            auc5=auc(dat1[,"y"],dat1[,"hSKHThlgr"])
            auc6=auc(dat1[,"y"],dat1[,"hsupphlgr"])
            
            
            out=rbind(out,c(year,t,u,auc1,auc2,auc3,auc4,auc5,auc6,horizon3))
          }
        }
        
        if (t == 2) {
          for(u in 6:6) {
            dat1=dat[dat[,"year"]==year & dat[,"t"]==t & dat[,"u"]==u,]
            
            auc1=auc(dat1[,"y"],dat1[,"hbench"])
            auc2=auc(dat1[,"y"],dat1[,"hsuppDTPO"])
            auc3=auc(dat1[,"y"],dat1[,"hsupp0hlgr"])
            auc4=auc(dat1[,"y"],dat1[,"hsephlgr"])
            auc5=auc(dat1[,"y"],dat1[,"hSKHThlgr"])
            auc6=auc(dat1[,"y"],dat1[,"hsupphlgr"])
            
            out=rbind(out,c(year,t,u,auc1,auc2,auc3,auc4,auc5,auc6,horizon3))
          }
        }
        
        
      }}
  }
  
  out <- out[-c(1, 2), ]
  
  if (horizon3 == 2) {
    out0=out[out[,"horizon3"]==0,]
    out1=out[out[,"horizon3"]==1,]
    
    print(list(apply(out0,2,mean),apply(out1,2,mean)))
    
    # par(mfrow=c(2,2))
    
    #b0=list(out0[,"superpp"],out0[,"schmid"],out0[,"frtrn"],out0[,"tu"],out0[,"average"])
    #boxplot(b0,names=c("superpp","schmid","frtrn","tu","average"),ylim=c(.43,1))
    
    #b1=list(out1[,"superpp"],out1[,"schmid"],out1[,"frtrn"],out1[,"tu"],out1[,"average"])
    #boxplot(b1,names=c("superpp","schmid","frtrn","tu","average"),ylim=c(.43,1))
    
    b=list(out1[,"hbench"],out0[,"hbench"],
           out1[,"hsuppDTPO"],out0[,"hsuppDTPO"],
           out1[,"hsupp0hlgr"],out0[,"hsupp0hlgr"],
           out1[,"hsephlgr"],out0[,"hsephlgr"],
           out1[,"hSKHThlgr"],out0[,"hSKHThlgr"],
           out1[,"hsupphlgr"],out0[,"hsupphlgr"])
    
    boxplot(b)
    abline(v=2.5)
    abline(v=4.5)
    abline(v=6.5)
    abline(v=8.5)
    abline(v=10.5)
  } else {
    out <- out[out[, "horizon3"] == horizon3, ]
    print(apply(out, 2, mean))
  }
  return(out)
  
}

