computecumgain=function(pred,y,ngroups){
  
  
  n=length(y)
  m=floor(n/ngroups)
  dat=data.frame(pred,y)
  
  names(dat)=c("pred","y")
  
  ny=sum(dat[,"y"])
  
  dat1=dat[order(-pred),]
  
  # print(dat1)
  
  out=matrix(0,ngroups,2)
  
  indmin=1
  
  for(i in 1:ngroups)
  {
    if(i==ngroups)
    {
      indmax=n
      #out[i,1]=indmax*100/(n)
      out[i,1]=100
      out[i,2]=100*sum(dat1[1:indmax,"y"])/ny
      #print(indmin)
      #print(indmax)
    }
    if(i!=ngroups)
    {
      #indmax=indmin+m-1
      indmax=ceiling(i/ngroups*n)
      #out[i,1]=indmax*100/(n)
      out[i,1]=i/ngroups*100
      out[i,2]=100*sum(dat1[1:indmax,"y"])/ny
      #print(indmin)
      #print(indmax)
      indmin=indmax+1
    }
    
  }
  
  # print(out)
  out
  
  
}