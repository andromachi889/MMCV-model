


SimulateVaccination <- function(params, B, mu, ContactMatrix, max_age, burnin, final_t){
  
  init <- c((1/12)*InitPop,(1/12)*InitPop,(1/12)*InitPop,(1/12)*InitPop,(1/12)*InitPop,(1/12)*InitPop,(1/12)*InitPop,(1/12)*InitPop,(1/12)*InitPop,
            (1/12)*InitPop,(1/12)*InitPop,(1/12)*InitPop,rep(0,200))
  x <- as.numeric(init)
  
  params["epsbeta"] = 0.6
  
  if(length(mass_t)==1){
    
    oot<-lsoda(x,seq(0,burnin*365,1.0),mmcv::PentaSt,hmax = 0.1,parms=params,Cm=ContactMatrix,mu=mu,max_age=max_age,B=B,ir_A=ir_A, ir_R=ir_R, stoc=stoc)
    
    init2 <- oot[dim(oot)[1],2:(dim(oot)[2]-2*max_age-3)]
    
    init3 <- c(init2,rep(0,(4*12+8)*max_age))
    oot2 <- lsoda(init3,seq(0,(mass_t-2000)*365+1,1),mmcv::pentavalent,hmax = 0.1,parms=params,Cm=ContactMatrix,mu=fin_mu,max_age=max_age,B=birth_rate,ir_A=ir_A, ir_R=ir_R, stoc=stoc, epi_A=epi_A,epi_R=epi_R)
    
    inits <- oot2[dim(oot2)[1],2:((5*12+8)*max_age+1)]
    init4 <- campaign129_V(inits,max_age,mcov)
    
    oot3 <- lsoda(init4,seq((mass_t-2000)*365+1,(routine_t-2000)*365+1,1),mmcv::pentavalent,hmax = 0.1,parms=params,Cm=ContactMatrix,mu=fin_mu,max_age=max_age,B=birth_rate,ir_A=ir_A, ir_R=ir_R, stoc=stoc, epi_A=epi_A,epi_R=epi_R)
    inits2 <- oot3[dim(oot3)[1],2:((5*12+8)*max_age+1)]
    
    init5 <- catchup_split(inits2,max_age,ccov, target = target)
    
    oot4 <- lsoda(init5,seq((routine_t-2000)*365+1,(penta_t-2000)*365+1,1),mmcv::pentavalent,hmax = 0.1,parms=params,Cm=ContactMatrix,mu=fin_mu,max_age=max_age,B=birth_rate,ir_A=ir_A, ir_R=ir_R, stoc=stoc, epi_A=epi_A,epi_R=epi_R)
    inits3 <- oot4[dim(oot4)[1],2:((5*12+8)*max_age+1)]
    
    init6 <- campaign_penta(inits3,max_age,ccov, target = target_penta)
    
    oot5 <- lsoda(init6,seq((penta_t-2000)*365+1,(final_t-2000+1)*365+1,1),mmcv::pentavalent,hmax = 0.1,parms=params,Cm=ContactMatrix,mu=fin_mu,max_age=max_age,B=birth_rate,ir_A=ir_A, ir_R=ir_R, stoc=stoc, epi_A=epi_A,epi_R=epi_R)
    
    
    outs <- rbind(oot2[-dim(oot2)[1],],oot3[-dim(oot3)[1],],oot4[-dim(oot4)[1],],oot5)
  }
  
  
  if(length(mass_t)==2){
    
    oot<-lsoda(x,seq(0,burnin*365,1.0),mmcv::PentaSt,hmax = 0.1,parms=params,Cm=ContactMatrix,mu=mu,max_age=max_age,B=B,ir_A=ir_A, ir_R=ir_R, stoc=stoc)
    
    init2 <- oot[dim(oot)[1],2:(dim(oot)[2]-2*max_age-3)]
    
    init3 <- c(init2,rep(0,(4*12+8)*max_age))
    oot2 <- lsoda(init3,seq(0,(mass_t[1]-2000)*365+1,1),mmcv::pentavalent,hmax = 0.1,parms=params,Cm=ContactMatrix,mu=fin_mu,max_age=max_age,B=birth_rate,ir_A=ir_A, ir_R=ir_R, stoc=stoc, epi_A=epi_A,epi_R=epi_R)
    
    inits <- oot2[dim(oot2)[1],2:((5*12+8)*max_age+1)]
    init4 <- campaign129_V(inits,max_age,mcov[1])
    
    
    oot3_1 <- lsoda(init4,seq((mass_t[1]-2000)*365+1,(mass_t[2]-2000)*365+1,1),mmcv::pentavalent,hmax = 0.1,parms=params,Cm=ContactMatrix,mu=fin_mu,max_age=max_age,B=birth_rate,ir_A=ir_A, ir_R=ir_R, stoc=stoc, epi_A=epi_A,epi_R=epi_R)
    inits2_1 <- oot3_1[dim(oot3_1)[1],2:((5*12+8)*max_age+1)]
    init4_1 <- catchup_split(inits2_1,max_age,mcov[2], target = 29)
    
    oot3 <- lsoda(init4_1,seq((mass_t[2]-2000)*365+1,(routine_t-2000)*365+1,1),mmcv::pentavalent,hmax = 0.1,parms=params,Cm=ContactMatrix,mu=fin_mu,max_age=max_age,B=birth_rate,ir_A=ir_A, ir_R=ir_R, stoc=stoc, epi_A=epi_A,epi_R=epi_R)
    inits2 <- oot3[dim(oot3)[1],2:((5*12+8)*max_age+1)]
    
    init5 <- catchup_split(inits2,max_age,ccov, target = target)
    
    oot4 <- lsoda(init5,seq((routine_t-2000)*365+1,(penta_t-2000)*365+1,1),mmcv::pentavalent,hmax = 0.1,parms=params,Cm=ContactMatrix,mu=fin_mu,max_age=max_age,B=birth_rate,ir_A=ir_A, ir_R=ir_R, stoc=stoc, epi_A=epi_A,epi_R=epi_R)
    inits3 <- oot4[dim(oot4)[1],2:((5*12+8)*max_age+1)]
    
    init6 <- campaign_penta(inits3,max_age,1, target = target_penta)
    
    oot5 <- lsoda(init6,seq((penta_t-2000)*365+1,(final_t-2000+1)*365+1,1),mmcv::pentavalent,hmax = 0.1,parms=params,Cm=ContactMatrix,mu=fin_mu,max_age=max_age,B=birth_rate,ir_A=ir_A, ir_R=ir_R, stoc=stoc, epi_A=epi_A,epi_R=epi_R)
    
    
    outs <- rbind(oot2[-dim(oot2)[1],],oot3_1[-dim(oot3_1)[1],],oot3[-dim(oot3)[1],],oot4[-dim(oot4)[1],],oot5)
  }
  

  
  out <- cases_routine(outs,max_age)


  return(out)
  
  
}