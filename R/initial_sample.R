initial_sample=function(p0,p1,p2,p3,diff=0,selection.prob=0.9){
  D=sqrt((p1*(1-p1)+p2*(1-p2)+p3*(1-p3))/3)
  quant=qmvnorm(selection.prob,mean=c(0,0),sigma=matrix(c(1,0.5,0.5,1),nrow=2),tail = "lower.tail")$quantile
  temp=ceiling(2*D^2*(quant+diff/sqrt(2))^2/(p1-p3)^2)
  U_1=sqrt(temp)*(p3-p1)/(sqrt(2)*D)-diff/sqrt(2)
  U_2=sqrt(temp)*(p3-p2)/(sqrt(2)*D)-diff/sqrt(2)

  while(pmvnorm(lower=c(-Inf,-Inf),upper=c(U_1,U_2),mean=c(0,0),sigma=matrix(c(1,0.5,0.5,1),nrow=2))[1]<selection.prob){
    temp=temp+1
    U_1=sqrt(temp)*(p3-p1)/(sqrt(2)*D)-diff/sqrt(2)
    U_2=sqrt(temp)*(p3-p2)/(sqrt(2)*D)-diff/sqrt(2)
  }
  return(temp)
}
