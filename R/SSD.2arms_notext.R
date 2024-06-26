SSD.2arms_notext <- function(r1, r, n1, n, p0, p1=NULL, p, nsim, diff=0.05, seed=0802) {

  set.seed(seed)

  if (n > 1000)
    stop ("sample size, n cannot exceed 1000")

  if (nsim <= 1) {
    stop(" nsim less than 2! ")
  }

  if (r < r1  ) {
    stop("r must be >= r1")
  }

  if (n <= n1) {
    stop("condition for n > n1 must be satisfied for a 2 stage design")
  }

  n2 <- n- n1

  p1<-ifelse(is.null(p1), p[1], p1)

  outcome2<-No.success<-n.Subj<-matrix(999, ncol=2, nrow=nsim)

  for(i in 1:nsim)	{

    for (a in 1:2) {

      Stage1<-rbinom(1, n1, p[a])

      outcome1<-ifelse(Stage1>r1, rbinom(1, n2, p[a]), NA)

      No.success[i,a]<-sum(outcome1, Stage1, na.rm=TRUE)

      outcome2[i,a]<-ifelse(outcome1+Stage1>r, 1, 0)

      n.Subj[i,a] <- ifelse(Stage1>r1, n, n1)

    }

  }


  Outcome<-apply(outcome2, 1, sum, na.rm=T)


  Outcome[is.na(Outcome)]<-0

  Prob.neg    <-length(Outcome[Outcome==0])/nsim
  Prob.pos    <-length(Outcome[Outcome==2])/nsim
  Prob.negpos <-length(Outcome[Outcome==1])/nsim

  Prob.ArmA<-sum(outcome2[,1], na.rm=TRUE)/nsim
  Prob.ArmB<-sum(outcome2[,2], na.rm=TRUE)/nsim

  mean.Subj<-apply(n.Subj,2,mean, na.rm=T)

  ### After 1st segment
  Prob.select.ArmA<-sum(outcome2[,1][outcome2[,2]==0 | is.na(outcome2[,2])],na.rm=TRUE)/nsim
  Prob.select.ArmB<-sum(outcome2[,2][outcome2[,1]==0 | is.na(outcome2[,1])],na.rm=TRUE)/nsim
  Prob.NoArm<-Prob.neg

  ### Selecting an Arm when both arms are positive  (2nd Segment)

  No.success.BothArms.select<-No.success[Outcome==2,]

  SSD.SelectArm<-function(x, diff)
  {
    NoArm<-ArmA<-ArmB<-NA
    if(x[2]/n==1&x[1]/n!=1){
      ArmB=1
      ArmA=0
    }
    if(x[2]/n!=1&x[1]/n==1){
      ArmB=0
      ArmA=1
    }
    if(x[2]/n==1&x[1]/n==1){
      if(diff==0){
        ArmA<-ifelse(runif(1,0,1)<0.5,0,1)
        ArmB<-1-ArmA
      }
      if(diff!=0){
        NoArm=1
      }
    }
    if(x[2]/n!=1&x[1]/n!=1){
      test_a=sqrt(n)*(x[1]/n-p0)/sqrt((x[1]/n)*(1-x[1]/n))
      test_b=sqrt(n)*(x[2]/n-p0)/sqrt((x[2]/n)*(1-x[2]/n))
      ArmB<-ifelse(test_b-test_a > diff,1,0)
      ArmA<-ifelse(test_b-test_a < -diff,1,0)
      ## for ties
      if(diff==0 & test_b-test_a  == 0) {
        ArmA<-ifelse(runif(1,0,1)<0.5,0,1)
        ArmB<-1-ArmA
      }

      ## for SSD_mod with no selection for ties or if diff <= e.g.0.05
      if(diff!=0) {
        NoArm<-ifelse(test_b-test_a <= diff & test_a-test_b <= diff,1,0)
      }
    }
    return(list(ArmA, ArmB, NoArm))
  }


  ### Original SSD

  if( length(No.success[Outcome==2,]) > 1 ) {
    SSD.SelectArm.2ndSeg<-matrix(unlist(
      apply(No.success[Outcome==2,],1,SSD.SelectArm, diff=0)   ###for each Two-stage successfull row, apply the selection criteria
    ),ncol=3,byrow=T)}   else {SSD.SelectArm.2ndSeg <- matrix(NA,ncol=3)  }

  ProbArmA.2ndSeg<-sum(SSD.SelectArm.2ndSeg[,1],na.rm=TRUE)/nsim
  ProbArmB.2ndSeg<-sum(SSD.SelectArm.2ndSeg[,2],na.rm=TRUE)/nsim
  ProbNoArm.2ndSeg<-sum(SSD.SelectArm.2ndSeg[,3],na.rm=TRUE)/nsim

  Overall.ArmA<-Prob.select.ArmA + ProbArmA.2ndSeg
  Overall.ArmB<-Prob.select.ArmB + ProbArmB.2ndSeg
  Overall.NoArm<-Prob.NoArm + ProbNoArm.2ndSeg


  ### Modified SSD ####

  if( length(No.success[Outcome==2,]) > 1 ) {
    SSD.SelectArm.2ndSeg<-matrix(unlist(
      apply(No.success[Outcome==2,],1,SSD.SelectArm, diff)
    ),ncol=3,byrow=T)}   else {SSD.SelectArm.2ndSeg <- matrix(NA,ncol=3)  }

  ProbArmA.2ndSeg.MOD<-sum(SSD.SelectArm.2ndSeg[,1],na.rm=TRUE)/nsim
  ProbArmB.2ndSeg.MOD<-sum(SSD.SelectArm.2ndSeg[,2],na.rm=TRUE)/nsim
  ProbNoArm.2ndSeg.MOD<-sum(SSD.SelectArm.2ndSeg[,3],na.rm=TRUE)/nsim

  Overall.ArmA.MOD<-Prob.select.ArmA + ProbArmA.2ndSeg.MOD
  Overall.ArmB.MOD<-Prob.select.ArmB + ProbArmB.2ndSeg.MOD
  Overall.NoArm.MOD<-Prob.NoArm + ProbNoArm.2ndSeg.MOD

  if(diff==0){
    soln<-data.frame("n"=n,
                     "SSD Arm A"=Overall.ArmA, "SSD Arm B"= Overall.ArmB, "SSD No Arm"=Overall.NoArm,
                     "diff"=diff,
                     "Mean N Arm A"=mean.Subj[1],"Mean N Arm B"=mean.Subj[2])


  }
  if(diff!=0){
    soln<-data.frame("n"=n,
                     "Modified SSD Arm A"=Overall.ArmA.MOD, "Modified SSD Arm B"=Overall.ArmB.MOD,
                     "Modified SSD No Arm"=Overall.NoArm.MOD, "diff"=diff,
                     "Mean N Arm A"=mean.Subj[1],"Mean N Arm B"=mean.Subj[2])


  }



  soln

}

