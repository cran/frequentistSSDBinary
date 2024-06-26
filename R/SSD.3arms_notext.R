
SSD.3arms_notext <- function(r1, r, n1, n, p0, p1=NULL, p, nsim, diff=0.05, seed=0802) {
  ProbNoArm.2ndSeg.MOD <- NULL
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

  outcome2<-No.success<-n.Subj<-matrix(999, ncol=3, nrow=nsim)

  for(i in 1:nsim)	{

    for (a in 1:3) {

      Stage1<-rbinom(1, n1, p[a])

      outcome1<-ifelse(Stage1>r1, rbinom(1, n2, p[a]), NA)

      No.success[i,a]<-sum(outcome1, Stage1, na.rm=TRUE)

      outcome2[i,a]<-ifelse(outcome1+Stage1>r, 1, 0)        ###NA is for failure of stage I ,0 is for success for stage I but fail at II, 1 means success

      n.Subj[i,a] <- ifelse(Stage1>r1, n, n1)

    }

  }


  Outcome<-apply(outcome2, 1, sum, na.rm=T)


  Outcome[is.na(Outcome)]<-0

  Prob.neg    <-length(Outcome[Outcome==0])/nsim
  Prob.pos    <-length(Outcome[Outcome==3])/nsim
  Prob.negpos <-length(Outcome[Outcome==1|Outcome==2])/nsim

  Prob.ArmA<-sum(outcome2[,1], na.rm=TRUE)/nsim
  Prob.ArmB<-sum(outcome2[,2], na.rm=TRUE)/nsim
  Prob.ArmC<-sum(outcome2[,3], na.rm=TRUE)/nsim

  mean.Subj<-apply(n.Subj,2,mean, na.rm=T)

  ### After 1st segment
  Ind_A=(outcome2[,2]==0 | is.na(outcome2[,2]))&(outcome2[,3]==0 | is.na(outcome2[,3]))
  Ind_B=(outcome2[,1]==0 | is.na(outcome2[,1]))&(outcome2[,3]==0 | is.na(outcome2[,3]))
  Ind_C=(outcome2[,1]==0 | is.na(outcome2[,1]))&(outcome2[,2]==0 | is.na(outcome2[,2]))
  Prob.select.ArmA<-sum(outcome2[,1][Ind_A],na.rm=TRUE)/nsim
  Prob.select.ArmB<-sum(outcome2[,2][Ind_B],na.rm=TRUE)/nsim
  Prob.select.ArmC<-sum(outcome2[,3][Ind_C],na.rm=TRUE)/nsim
  Prob.NoArm<-Prob.neg

  ### Selecting an Arm when both arms or three arms are positive  (2nd Segment)
  AB_ind=which((outcome2[,1]==1&outcome2[,2]==1&outcome2[,3]==0)|(outcome2[,1]==1&outcome2[,2]==1&is.na(outcome2[,3])))
  BC_ind=which((outcome2[,3]==1&outcome2[,2]==1&outcome2[,1]==0)|(outcome2[,3]==1&outcome2[,2]==1&is.na(outcome2[,1])))
  AC_ind=which((outcome2[,1]==1&outcome2[,3]==1&outcome2[,2]==0)|(outcome2[,1]==1&outcome2[,3]==1&is.na(outcome2[,2])))


  ##  this is for two postive arms
  SSD.SelectArm_Two<-function(x, diff){
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

  ###this is for three postive arms

  SSD.SelectArm_Three<-function(x, diff){
    NoArm<-ArmA<-ArmB<-ArmC<-NA

    if(diff!=0){
      if(x[1]/n==1&x[2]/n==1&x[3]/n==1){
        NoArm=1
      }
      if(x[1]/n==1&x[2]/n==1&x[3]/n!=1){
        NoArm=1
      }
      if(x[1]/n==1&x[2]/n!=1&x[3]/n==1){
        NoArm=1
      }
      if(x[1]/n!=1&x[2]/n==1&x[3]/n==1){
        NoArm=1
      }
    }
    if(diff==0){
      if(x[1]/n==1&x[2]/n==1&x[3]/n==1){
        temp=rmultinom(1,1 , prob=c(1/3,1/3,1/3))
        ArmA<-temp[1]
        ArmB<-temp[2]
        ArmC<-temp[3]
      }
      if(x[1]/n==1&x[2]/n==1&x[3]/n!=1){
        ArmA<-ifelse(runif(1,0,1)<0.5,0,1)
        ArmB<-1-ArmA
      }
      if(x[1]/n==1&x[2]/n!=1&x[3]/n==1){
        ArmA<-ifelse(runif(1,0,1)<0.5,0,1)
        ArmC<-1-ArmA
      }
      if(x[1]/n!=1&x[2]/n==1&x[3]/n==1){
        ArmB<-ifelse(runif(1,0,1)<0.5,0,1)
        ArmC<-1-ArmB
      }
    }

    if(x[1]/n==1&x[2]/n!=1&x[3]/n!=1){
      ArmA=1}
    if(x[1]/n!=1&x[2]/n==1&x[3]/n!=1){
      ArmB=1}
    if(x[1]/n!=1&x[2]/n!=1&x[3]/n==1){
      ArmC=1}



    if(x[2]/n!=1&x[1]/n!=1&x[3]/n!=1){
      test_a=sqrt(n)*(x[1]/n-p0)/sqrt((x[1]/n)*(1-x[1]/n))
      test_b=sqrt(n)*(x[2]/n-p0)/sqrt((x[2]/n)*(1-x[2]/n))
      test_c=sqrt(n)*(x[3]/n-p0)/sqrt((x[3]/n)*(1-x[3]/n))

      selection_A=(test_a-test_b > diff)&(test_a-test_c > diff)
      selection_B=(test_b-test_a > diff)&(test_b-test_c > diff)
      selection_C=(test_c-test_a > diff)&(test_c-test_b> diff)

      ArmB<-ifelse(selection_B,1,0)
      ArmA<-ifelse(selection_A,1,0)
      ArmC<-ifelse(selection_C,1,0)
      ## for ties
      if(diff==0 & (test_b==test_a)&(test_b==test_c)) {
        temp=rmultinom(1,1 , prob=c(1/3,1/3,1/3))
        ArmA<-temp[1]
        ArmB<-temp[2]
        ArmC<-temp[3]
      }
      if(diff==0 & (test_b==test_a)&(test_b>test_c)) {
        ArmA<-ifelse(runif(1,0,1)<0.5,0,1)
        ArmB<-1-ArmA
        ArmC<-0
      }
      if(diff==0 & (test_c==test_a)&(test_c>test_b)) {
        ArmA<-ifelse(runif(1,0,1)<0.5,0,1)
        ArmB<-0
        ArmC<-1-ArmA
      }
      if(diff==0 & (test_b==test_c)&(test_b>test_a)) {
        ArmB<-ifelse(runif(1,0,1)<0.5,0,1)
        ArmC<-1-ArmB
        ArmA<-0
      }

      ## for SSD_mod with no selection for ties or if diff <= e.g.0.05
      if(diff!=0) {
        NoArm<-ifelse(sum(selection_A,selection_B,selection_C)==0,1,0)
      }
    }
    return(list(ArmA, ArmB,ArmC, NoArm))
  }


  ### Original SSD

  ##########################Two arms###########################

  ##for AB arm postive

  if( length(No.success[AB_ind,c(1,2)]) > 1 ) {
    SSD.SelectArm.2nd.AB.Seg<-matrix(unlist(
      apply(matrix(No.success[AB_ind,c(1,2)],nrow=length(No.success[AB_ind,c(1,2)])/2),1,SSD.SelectArm_Two, diff=0)   ###for each Two-stage successfull row, apply the selection criteria
    ),ncol=3,byrow=T)}   else {SSD.SelectArm.2nd.AB.Seg <- matrix(NA,ncol=3)  }

  ProbArmA.2nd.AB.Seg<-sum(SSD.SelectArm.2nd.AB.Seg[,1],na.rm=TRUE)/nsim
  ProbArmB.2nd.AB.Seg<-sum(SSD.SelectArm.2nd.AB.Seg[,2],na.rm=TRUE)/nsim
  ProbNoArm.2nd.AB.Seg<-sum(SSD.SelectArm.2nd.AB.Seg[,3],na.rm=TRUE)/nsim

  ##for BC arm postive

  if( length(No.success[BC_ind,c(2,3)]) > 1 ) {
    SSD.SelectArm.2nd.BC.Seg<-matrix(unlist(
      apply(matrix(No.success[BC_ind,c(1,2)],nrow=length(No.success[BC_ind,c(1,2)])/2),1,SSD.SelectArm_Two, diff=0)   ###for each Two-stage successfull row, apply the selection criteria
    ),ncol=3,byrow=T)}   else {SSD.SelectArm.2nd.BC.Seg <- matrix(NA,ncol=3)  }

  ProbArmB.2nd.BC.Seg<-sum(SSD.SelectArm.2nd.BC.Seg[,1],na.rm=TRUE)/nsim
  ProbArmC.2nd.BC.Seg<-sum(SSD.SelectArm.2nd.BC.Seg[,2],na.rm=TRUE)/nsim
  ProbNoArm.2nd.BC.Seg<-sum(SSD.SelectArm.2nd.BC.Seg[,3],na.rm=TRUE)/nsim


  ##for AC arm postive

  if( length(No.success[AC_ind,c(1,3)]) > 1 ) {
    SSD.SelectArm.2nd.AC.Seg<-matrix(unlist(
      apply(matrix(No.success[AC_ind,c(1,2)],nrow=length(No.success[AC_ind,c(1,2)])/2),1,SSD.SelectArm_Two, diff=0)   ###for each Two-stage successfull row, apply the selection criteria
    ),ncol=3,byrow=T)}   else {SSD.SelectArm.2nd.AC.Seg <- matrix(NA,ncol=3)  }

  ProbArmA.2nd.AC.Seg<-sum(SSD.SelectArm.2nd.AC.Seg[,1],na.rm=TRUE)/nsim
  ProbArmC.2nd.AC.Seg<-sum(SSD.SelectArm.2nd.AC.Seg[,2],na.rm=TRUE)/nsim
  ProbNoArm.2nd.AC.Seg<-sum(SSD.SelectArm.2nd.AC.Seg[,3],na.rm=TRUE)/nsim

  ########################Three arms
  ###for ABC postive

  if( length(No.success[Outcome==3,]) > 1 ) {
    SSD.SelectArm.2nd.ABC.Seg<-matrix(unlist(
      apply(matrix(No.success[Outcome==3,],nrow=length(No.success[Outcome==3,])/3),1,SSD.SelectArm_Three, diff=0)   ###for each Two-stage successfull row, apply the selection criteria
    ),ncol=4,byrow=T)}   else {SSD.SelectArm.2nd.ABC.Seg <- matrix(NA,ncol=4)  }

  ProbArmA.2nd.ABC.Seg<-sum(SSD.SelectArm.2nd.ABC.Seg[,1],na.rm=TRUE)/nsim
  ProbArmB.2nd.ABC.Seg<-sum(SSD.SelectArm.2nd.ABC.Seg[,2],na.rm=TRUE)/nsim
  ProbArmC.2nd.ABC.Seg<-sum(SSD.SelectArm.2nd.ABC.Seg[,3],na.rm=TRUE)/nsim
  ProbNoArm.2nd.ABC.Seg<-sum(SSD.SelectArm.2nd.ABC.Seg[,4],na.rm=TRUE)/nsim



  Overall.ArmA<-Prob.select.ArmA + ProbArmA.2nd.AB.Seg+ProbArmA.2nd.AC.Seg+ProbArmA.2nd.ABC.Seg
  Overall.ArmB<-Prob.select.ArmB + ProbArmB.2nd.AB.Seg+ProbArmB.2nd.BC.Seg+ProbArmB.2nd.ABC.Seg
  Overall.ArmC<-Prob.select.ArmC + ProbArmC.2nd.AC.Seg+ProbArmC.2nd.BC.Seg+ProbArmC.2nd.ABC.Seg
  Overall.NoArm<-Prob.NoArm + ProbNoArm.2nd.AB.Seg+ProbNoArm.2nd.BC.Seg+ProbNoArm.2nd.AC.Seg+ProbNoArm.2nd.ABC.Seg


  ### Modified SSD ####



  ##########################Two arms###########################

  ##for AB arm postive

  if( length(No.success[AB_ind,c(1,2)]) > 1 ) {
    SSD.SelectArm.2nd.AB.Seg<-matrix(unlist(
      apply(matrix(No.success[AB_ind,c(1,2)],nrow=length(No.success[AB_ind,c(1,2)])/2),1,SSD.SelectArm_Two, diff=diff)   ###for each Two-stage successfull row, apply the selection criteria
    ),ncol=3,byrow=T)}   else {SSD.SelectArm.2nd.AB.Seg <- matrix(NA,ncol=3)  }

  ProbArmA.2nd.AB.Seg.MOD<-sum(SSD.SelectArm.2nd.AB.Seg[,1],na.rm=TRUE)/nsim
  ProbArmB.2nd.AB.Seg.MOD<-sum(SSD.SelectArm.2nd.AB.Seg[,2],na.rm=TRUE)/nsim
  ProbNoArm.2nd.AB.Seg.MOD<-sum(SSD.SelectArm.2nd.AB.Seg[,3],na.rm=TRUE)/nsim

  ##for BC arm postive

  if( length(No.success[BC_ind,c(2,3)]) > 1 ) {
    SSD.SelectArm.2nd.BC.Seg<-matrix(unlist(
      apply(matrix(No.success[BC_ind,c(2,3)],nrow=length(No.success[BC_ind,c(2,3)])/2),1,SSD.SelectArm_Two, diff=diff)   ###for each Two-stage successfull row, apply the selection criteria
    ),ncol=3,byrow=T)}   else {SSD.SelectArm.2nd.BC.Seg <- matrix(NA,ncol=3)  }

  ProbArmB.2nd.BC.Seg.MOD<-sum(SSD.SelectArm.2nd.BC.Seg[,1],na.rm=TRUE)/nsim
  ProbArmC.2nd.BC.Seg.MOD<-sum(SSD.SelectArm.2nd.BC.Seg[,2],na.rm=TRUE)/nsim
  ProbNoArm.2nd.BC.Seg.MOD<-sum(SSD.SelectArm.2nd.BC.Seg[,3],na.rm=TRUE)/nsim


  ##for AC arm postive

  if( length(No.success[AC_ind,c(1,3)]) > 1 ) {
    SSD.SelectArm.2nd.AC.Seg<-matrix(unlist(
      apply(matrix(No.success[AC_ind,c(1,3)],nrow=length(No.success[AC_ind,c(1,3)])/2),1,SSD.SelectArm_Two, diff=diff)   ###for each Two-stage successfull row, apply the selection criteria
    ),ncol=3,byrow=T)}   else {SSD.SelectArm.2nd.AC.Seg <- matrix(NA,ncol=3)  }

  ProbArmA.2nd.AC.Seg.MOD<-sum(SSD.SelectArm.2nd.AC.Seg[,1],na.rm=TRUE)/nsim
  ProbArmC.2nd.AC.Seg.MOD<-sum(SSD.SelectArm.2nd.AC.Seg[,2],na.rm=TRUE)/nsim
  ProbNoArm.2nd.AC.Seg.MOD<-sum(SSD.SelectArm.2nd.AC.Seg[,3],na.rm=TRUE)/nsim

  ########################Three arms
  ###for ABC postive

  if( length(No.success[Outcome==3,]) > 1 ) {
    SSD.SelectArm.2nd.ABC.Seg<-matrix(unlist(
      apply(matrix(No.success[Outcome==3,],nrow=length(No.success[Outcome==3,])/3),1,SSD.SelectArm_Three, diff=diff)   ###for each Two-stage successfull row, apply the selection criteria
    ),ncol=4,byrow=T)}   else {SSD.SelectArm.2nd.ABC.Seg <- matrix(NA,ncol=4)  }

  ProbArmA.2nd.ABC.Seg.MOD<-sum(SSD.SelectArm.2nd.ABC.Seg[,1],na.rm=TRUE)/nsim
  ProbArmB.2nd.ABC.Seg.MOD<-sum(SSD.SelectArm.2nd.ABC.Seg[,2],na.rm=TRUE)/nsim
  ProbArmC.2nd.ABC.Seg.MOD<-sum(SSD.SelectArm.2nd.ABC.Seg[,3],na.rm=TRUE)/nsim
  ProbNoArm.2nd.ABC.Seg.MOD<-sum(SSD.SelectArm.2nd.ABC.Seg[,4],na.rm=TRUE)/nsim



  Overall.ArmA.MOD<-Prob.select.ArmA + ProbArmA.2nd.AB.Seg.MOD+ProbArmA.2nd.AC.Seg.MOD+ProbArmA.2nd.ABC.Seg.MOD
  Overall.ArmB.MOD<-Prob.select.ArmB + ProbArmB.2nd.AB.Seg.MOD+ProbArmB.2nd.BC.Seg.MOD+ProbArmB.2nd.ABC.Seg.MOD
  Overall.ArmC.MOD<-Prob.select.ArmC + ProbArmC.2nd.AC.Seg.MOD+ProbArmC.2nd.BC.Seg.MOD+ProbArmC.2nd.ABC.Seg.MOD
  Overall.NoArm.MOD<-Prob.NoArm + ProbNoArm.2nd.AB.Seg.MOD+ProbNoArm.2nd.BC.Seg.MOD+ProbNoArm.2nd.AC.Seg.MOD+ProbNoArm.2nd.ABC.Seg.MOD







  if(diff==0){
    soln<-data.frame("n"=n,
                     "SSD Arm A"=Overall.ArmA, "SSD Arm B"= Overall.ArmB,"SSD Arm C"= Overall.ArmC, "SSD No Arm"=Overall.NoArm,
                     "diff"=diff,
                     "Mean N Arm A"=mean.Subj[1],"Mean N Arm B"=mean.Subj[2],"Mean N Arm C"=mean.Subj[3])


  }
  if(diff!=0){
    soln<-data.frame("n"=n,
                     "Modified SSD Arm A"=Overall.ArmA.MOD, "Modified SSD Arm B"=Overall.ArmB.MOD, "Modified SSD Arm C"=Overall.ArmC.MOD,
                     "Modified SSD No Arm"=Overall.NoArm.MOD, "diff"=diff,
                     "Mean N Arm A"=mean.Subj[1],"Mean N Arm B"=mean.Subj[2],"Mean N Arm C"=mean.Subj[3])


  }



  soln

}
