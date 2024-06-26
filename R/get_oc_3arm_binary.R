get_oc_3arm_binary=function(r1, r2, n1, n, p0,p1,p2,p3,diff=0,nsim, seed = 2483){
  SSD.3arms_notext(r1, r=r2, n1, n, p0,p=c(p1,p2,p3),diff=diff,nsim=nsim, seed = seed)
}
