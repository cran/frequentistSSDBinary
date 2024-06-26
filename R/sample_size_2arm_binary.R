sample_size_2arm_binary=function(p0,p1,p2,diff=0,selection.prob=0.9,alpha=0.1,beta=0.2){
  nsim <- 1000
  D=sqrt(p1*(1-p1)/2+(p2*(1-p2))/2)
  initial_n=ceiling(2*D^2*(qnorm(selection.prob)+diff/sqrt(2))^2/(p1-p2)^2)
  simon2stage=ph2simon(p0, p1, alpha, beta, nmax=100)
  ind=which(simon2stage$out[,4]>=initial_n)[1]
  message("The sample size for each arm is ","\n")
  if(diff==0){
    i<-0
    para<-simon2stage$out[ind+i,]
    result<-SSD.2arms_notext(r1=para[1:4][1],n1=para[1:4][2],r=para[1:4][3],n=para[1:4][4],diff=diff,p0,p=c(p1,p2),nsim=nsim)
    while(result$SSD.Arm.B<selection.prob){
      para<-simon2stage$out[ind+i,]
      result<-SSD.2arms_notext(r1=para[1:4][1],n1=para[1:4][2],r=para[1:4][3],n=para[1:4][4],diff=diff,p0,p=c(p1,p2),nsim=nsim)
      i=i+1
      if((initial_n+i)>100){#previous code is nmax
        stop("Exceed the max sample size")

      }
    }
    return(result$n)
  }
  if(diff!=0){
    i<-0
    para<-simon2stage$out[ind+i,]
    result<-SSD.2arms_notext(r1=para[1:4][1],n1=para[1:4][2],r=para[1:4][3],n=para[1:4][4],diff=diff,p0,p=c(p1,p2),nsim=nsim)
    while(result$Modified.SSD.Arm.B<selection.prob){
      para<-simon2stage$out[ind+i,]
      result<-SSD.2arms_notext(r1=para[1:4][1],n1=para[1:4][2],r=para[1:4][3],n=para[1:4][4],diff=diff,p0,p=c(p1,p2),nsim=nsim)
      i=i+1
      if((initial_n+i)>100){#previous code is nmax
        stop("Exceed the max sample size")

      }
    }
    return(result$n)
  }
}
