computeroc=function(pred,y,ngroups){
  
  
  n=length(y)
  # m=floor(n/ngroups)
  dat=data.frame(pred,y)
  
  names(dat)=c("pred","y")
  
  # ny=sum(dat[,"y"])
  
  # dat1=dat[order(-pred),]
  
  # print(dat1)
  
  out=matrix(0,ngroups+1,2)
  
  for(i in 0:ngroups)
  {
    cut=i/ngroups
    dat[,"predbin"]=(dat[,"pred"]>cut)+0
    a=NROW(dat[dat[,"y"]==0 & dat[,"predbin"]==0,])
    b=NROW(dat[dat[,"y"]==0 & dat[,"predbin"]==1,])
    cc=NROW(dat[dat[,"y"]==1 & dat[,"predbin"]==0,])
    d=NROW(dat[dat[,"y"]==1 & dat[,"predbin"]==1,])
    
    out[i+1,1]=b*100/(a+b)
    out[i+1,2]=d*100/(cc+d)
    
  }
  
  # plot(out)
  
  print(out)
  
  out
  
  
}