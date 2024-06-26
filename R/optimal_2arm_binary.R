optimal_2arm_binary=function(p0,p1,p2,alpha=0.1,beta=0.2,tot_sample){
  simon2stage=ph2simon(p0, p1, alpha, beta, nmax=100)   ##sometime, if we set nmax=tot_sample, it may return error
  ind=(simon2stage$out[,4]==tot_sample)
  para=simon2stage$out[ind,1:6]
  typeI=as.numeric(binom.power(r1=para[1:4][1],n1=para[1:4][2],r=para[1:4][3],n=para[1:4][4],p0))
  power=as.numeric(binom.power(r1=para[1:4][1],n1=para[1:4][2],r=para[1:4][3],n=para[1:4][4],p1))
  message("$Two_Stage ","\n")
  result=data.frame( numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0), numeric(0),numeric(0))
  colnames(result)=c("alpha","beta","r1","n1","r2","n","ESS","PS")
  result[1,]=c(typeI,1-power,para)
  return(round(result, 4))
}
