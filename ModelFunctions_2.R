require(rootSolve)
require(deSolve)
require(Rcpp)
library(doParallel)
library(parallel)
library(foreach)
library(doSNOW)
library(mmcv)



'%!in%' <- function(x,y)!('%in%'(x,y))

#helper function to pack 2d matrix into array

index2 <- function(i,a,y)
{
  return (y*(i-1) + a)
}



Expand <- function(Cm,ageGroupBreaks,max_age)
{
  
  # Function to expand coarse-grained age matrix Cm
  # to annual age-cohorts
  
  # Input: Cm, contact matrix with age groups defined by
  #        ageGroupBreaks (vector)
  #        max_age (numeric)
  
  # Output: Matrix with annual age groups ( max_age, max_age)
  
  CmEx <- matrix(0,max_age,max_age)
  
  # ageGroupBreaks <- c(1,ageGroupBreaks,max_age)
  
  for(i in 1:(length(ageGroupBreaks)-1) )
  {
    for(j in 1:(length(ageGroupBreaks)-1))
    {
      for(a1 in seq(ageGroupBreaks[i],ageGroupBreaks[i+1],1))
      {
        for(a2 in seq(ageGroupBreaks[j],ageGroupBreaks[j+1],1))
        {
          CmEx[a1,a2] = Cm[i,j];
        }
      }
    }
  }
  return(CmEx)
}


ExpandInvasion <- function(rate,AgB,max_age){
  
  inf = numeric(max_age)
  inf[1] = rate[1]
  
for(i in 2:length(AgB)){
  inf[(AgB[i-1]+1):(AgB[i])] = rate[i]
}
  return(inf)
}



mk_initial_conditions <- function(max_age,params,Cm,mu,B, stoc)
{
  
  a = 1/365.0
  
  x = numeric(12*max_age)
  
  for(k in c(1:12))
  {
    x[index2(k,1,max_age)] = (B/(mu[1] + a))/12.0;
  }
  
  for( a1 in c(2:(max_age-1)))
  {
    for(k in c(1:12))
    {
      x[index2(k,a1,max_age)] = ((a*x[index2(k,a1-1,max_age)])/(mu[a1]+a));
    }
  }
  
  for(k in c(1:12))
  {
    x[index2(k,max_age,max_age)] = (a*x[index2(k,max_age-1,max_age)]/(mu[max_age]));
  }
  
  x2 <- c(x,rep(0,200))
  
  FPs<-steady(x2,0,mmcv::PentaSt,parms=params,positive=TRUE,Cm=Cm,mu=mu,max_age=max_age,B=B,ir_A=ir_A,ir_R=ir_R,stoc=stoc)$y
  
  # dim(FPs) <- c(max_age,12)
  
  return(cbind(FPs,0))
}



cases <- function(x,max_age){
  
  
  
  datA <- x[,(12*max_age+2):((12*max_age+2)+max_age-1)]
  datR <- x[,(13*max_age+2):((13*max_age+2)+max_age-1)]
  
  x[,1] = floor(x[,1]/365)
  
  testA = matrix(0,length(unique(x[,1])), max_age)
  testR = matrix(0,length(unique(x[,1])), max_age)
  
  for (i in 1:ncol(datA)){
    
    datA[,i] = c(0,diff(datA[,i]))
    datR[,i] = c(0,diff(datR[,i]))
    
    yA <- data.frame(time = x[,1], inc = datA[,i])
    yR <- data.frame(time = x[,1], inc = datR[,i])
    
    testA[,i] = aggregate(inc~time, data = yA, FUN = sum)$inc
    testR[,i] = aggregate(inc~time, data = yR, FUN = sum)$inc
    
  }
  
  test <- cbind(unique(x[,1]), testA,testR)
  return(test)
  
  
}




cases_routine <- function(x,max_age){
  
  
  
  datA <- x[,(63*max_age+2):((63*max_age+2)+max_age-1)]
  datR <- x[,(67*max_age+2):((67*max_age+2)+max_age-1)]
  
  x[,1] = floor(x[,1]/365)
  
  testA = matrix(0,length(unique(x[,1])), max_age)
  testR = matrix(0,length(unique(x[,1])), max_age)
  
  for (i in 1:ncol(datA)){
    
    datA[,i] = c(0,diff(datA[,i]))
    datR[,i] = c(0,diff(datR[,i]))
    
    yA <- data.frame(time = x[,1], inc = datA[,i])
    yR <- data.frame(time = x[,1], inc = datR[,i])
    
    testA[,i] = aggregate(inc~time, data = yA, FUN = sum)$inc
    testR[,i] = aggregate(inc~time, data = yR, FUN = sum)$inc
    
  }
  
  test <- cbind(unique(x[,1]), testA,testR)
  return(test)
  
  
}






##vaccination coverage for routine EPI

epirate <- function(x){
  x <- x[x$country_code==country & x$year>=2000 & x$activity_type=='routine',]
  y <- order(x$year)[-dim(x)[1]]
  x2 <- x[y,]
  
dat <- x2[,'coverage']
  
  return(dat)
  
}


carriage <- function(x){
  carA <- cbind(x[,102:201],x[,702:801])
  carriersA <- apply(carA,1,sum)
  
  carR <- cbind(x[,202:301],x[,802:901])
  carriersR <- apply(carR,1,sum)
  
  df <- data.frame(carriersA = carriersA, carriersR = carriersR, N = x[,1503])
  
  return(df)
}





carriage_vaccinated <- function(x){
  carA <- cbind(x[,102:201],x[,702:801],x[,1302:1401],x[,1902:2001],x[,2502:2601],x[,3102:3201],x[,3702:3801],x[,4302:4401],x[,4902:5001],x[,5502:5601])
  carriersA <- apply(carA,1,sum)
  
  carR <- cbind(x[,202:301],x[,802:901],x[,1402:1501],x[,2002:2101],x[,2602:2701],x[,3202:3301],x[,3802:3901],x[,4402:4501],x[,5002:5101],x[,5602:5701])
  carriersR <- apply(carR,1,sum)
  
  df <- data.frame(carriersA = carriersA, carriersR = carriersR, N = x[,6804])
  
  return(df)
}





# campaign129_V <- function(x,max_age,cov){
#   y <- vector("numeric")
#   cov_short <- c(0,rep(cov,4), rep(0,(max_age-5)))
#   cov_long <- c(0,0,0,0,0,rep(cov,25),rep(0,(max_age-30)))
#   
#   y[1:max_age] = (1-cov_short-cov_long)*x[1:(max_age)] 
#   y[(max_age+1):(2*max_age)] = (1-cov_short-cov_long)*x[(max_age+1):(2*max_age)]
#   y[(2*max_age+1):(3*max_age)] = (1-cov_short-cov_long)*x[(2*max_age+1):(3*max_age)]
#   y[(3*max_age+1):(4*max_age)] = x[(3*max_age+1):(4*max_age)]
#   y[(4*max_age+1):(5*max_age)] = x[(4*max_age+1):(5*max_age)]
#   y[(5*max_age+1):(6*max_age)] = (1-cov_short-cov_long)*x[(5*max_age+1):(6*max_age)]
#   y[(6*max_age+1):(7*max_age)] = (1-cov_short-cov_long)*x[(6*max_age+1):(7*max_age)]
#   y[(7*max_age+1):(8*max_age)] = (1-cov_short-cov_long)*x[(7*max_age+1):(8*max_age)]
#   y[(8*max_age+1):(9*max_age)] = (1-cov_short-cov_long)*x[(8*max_age+1):(9*max_age)]
#   y[(9*max_age+1):(10*max_age)] = x[(9*max_age+1):(10*max_age)]
#   y[(10*max_age+1):(11*max_age)] = x[(10*max_age+1):(11*max_age)]
#   y[(11*max_age+1):(12*max_age)] = (1-cov_short-cov_long)*x[(11*max_age+1):(12*max_age)]
#   
#   
#   
#   y[(2*max_age+1):(3*max_age)] = x[(2*max_age+1):(3*max_age)]
#   y[(3*max_age+1):(4*max_age)] = (1-cov_short-cov_long)*x[(3*max_age+1):(4*max_age)]
#   
#   y[(4*max_age+1):(5*max_age)] = cov_long*x[1:(max_age)]  
#   y[(5*max_age+1):(6*max_age)] = cov_long*x[(max_age+1):(2*max_age)]
#   y[(6*max_age+1):(7*max_age)] = 0
#   y[(7*max_age+1):(8*max_age)] = cov_long*x[(3*max_age+1):(4*max_age)]
#   
#   
#   y[(8*max_age+1):(9*max_age)] = cov_short*x[1:(max_age)]
#   y[(9*max_age+1):(10*max_age)] = cov_short*x[(max_age+1):(2*max_age)]
#   y[(10*max_age+1):(11*max_age)] = 0
#   y[(11*max_age+1):(12*max_age)] = cov_short*x[(3*max_age+1):(4*max_age)]
#   
#   
#   
#   y[(12*max_age+1):(13*max_age)] = 0
#   
#   return(y)
# }



campaign129_V <- function(x,max_age,cov){
  y <- vector("numeric")
  cov_short <- c(0,rep(cov,4), rep(0,(max_age-5)))
  cov_long <- c(0,0,0,0,0,rep(cov,25),rep(0,(max_age-30)))
  
  #not vaccinated
  
  y[1:max_age] = (1-cov_short-cov_long)*x[1:(max_age)] 
  y[(max_age+1):(2*max_age)] = (1-cov_short-cov_long)*x[(max_age+1):(2*max_age)]
  y[(2*max_age+1):(3*max_age)] = (1-cov_short-cov_long)*x[(2*max_age+1):(3*max_age)]
  y[(3*max_age+1):(4*max_age)] = x[(3*max_age+1):(4*max_age)]
  y[(4*max_age+1):(5*max_age)] =  x[(4*max_age+1):(5*max_age)] 
  y[(5*max_age+1):(6*max_age)] = (1-cov_short-cov_long)*x[(5*max_age+1):(6*max_age)]
  y[(6*max_age+1):(7*max_age)] = (1-cov_short-cov_long)*x[(6*max_age+1):(7*max_age)]
  y[(7*max_age+1):(8*max_age)] = (1-cov_short-cov_long)*x[(7*max_age+1):(8*max_age)]
  y[(8*max_age+1):(9*max_age)] = (1-cov_short-cov_long)*x[(8*max_age+1):(9*max_age)]
  y[(9*max_age+1):(10*max_age)] = x[(9*max_age+1):(10*max_age)]
  y[(10*max_age+1):(11*max_age)] = x[(10*max_age+1):(11*max_age)] 
  y[(11*max_age+1):(12*max_age)] = (1-cov_short-cov_long)*x[(11*max_age+1):(12*max_age)]
  
  #MAV long  
  
  y[(12*max_age+1):(13*max_age)] = cov_long*x[(0*max_age+1):(1*max_age)] 
  y[(13*max_age+1):(14*max_age)] = cov_long*x[(1*max_age+1):(2*max_age)]
  y[(14*max_age+1):(15*max_age)] = cov_long*x[(2*max_age+1):(3*max_age)]
  y[(15*max_age+1):(16*max_age)] = 0
  y[(16*max_age+1):(17*max_age)] = 0
  y[(17*max_age+1):(18*max_age)] = cov_long*x[(5*max_age+1):(6*max_age)]
  y[(18*max_age+1):(19*max_age)] = cov_long*x[(6*max_age+1):(7*max_age)]
  y[(19*max_age+1):(20*max_age)] = cov_long*x[(7*max_age+1):(8*max_age)]
  y[(20*max_age+1):(21*max_age)] = cov_long*x[(8*max_age+1):(9*max_age)]
  y[(21*max_age+1):(22*max_age)] = 0
  y[(22*max_age+1):(23*max_age)] = 0
  y[(23*max_age+1):(24*max_age)] = cov_long*x[(11*max_age+1):(12*max_age)]
  
  #MAV short  
  
  y[(24*max_age+1):(25*max_age)] = cov_short*x[(0*max_age+1):(1*max_age)]   
  y[(25*max_age+1):(26*max_age)] = cov_short*x[(1*max_age+1):(2*max_age)]
  y[(26*max_age+1):(27*max_age)] = cov_short*x[(2*max_age+1):(3*max_age)]
  y[(27*max_age+1):(28*max_age)] = 0
  y[(28*max_age+1):(29*max_age)] = 0
  y[(29*max_age+1):(30*max_age)] = cov_short*x[(5*max_age+1):(6*max_age)]
  y[(30*max_age+1):(31*max_age)] = cov_short*x[(6*max_age+1):(7*max_age)]
  y[(31*max_age+1):(32*max_age)] = cov_short*x[(7*max_age+1):(8*max_age)]
  y[(32*max_age+1):(33*max_age)] = cov_short*x[(8*max_age+1):(9*max_age)]
  y[(33*max_age+1):(34*max_age)] = 0
  y[(34*max_age+1):(35*max_age)] = 0
  y[(35*max_age+1):(36*max_age)] = cov_short*x[(11*max_age+1):(12*max_age)]
  
  #Penta long  
  
  y[(36*max_age+1):(37*max_age)] = 0   
  y[(37*max_age+1):(38*max_age)] = 0
  y[(38*max_age+1):(39*max_age)] = 0
  y[(39*max_age+1):(40*max_age)] = 0
  y[(40*max_age+1):(41*max_age)] = 0
  y[(41*max_age+1):(42*max_age)] = 0
  y[(42*max_age+1):(43*max_age)] = 0
  y[(43*max_age+1):(44*max_age)] = 0
  y[(44*max_age+1):(45*max_age)] = 0
  y[(45*max_age+1):(46*max_age)] = 0
  y[(46*max_age+1):(47*max_age)] = 0
  y[(47*max_age+1):(48*max_age)] = 0
  
  #Penta short  
  
  y[(48*max_age+1):(49*max_age)] = 0   
  y[(49*max_age+1):(50*max_age)] = 0
  y[(50*max_age+1):(51*max_age)] = 0
  y[(51*max_age+1):(52*max_age)] = 0
  y[(52*max_age+1):(53*max_age)] = 0
  y[(53*max_age+1):(54*max_age)] = 0
  y[(54*max_age+1):(55*max_age)] = 0
  y[(55*max_age+1):(56*max_age)] = 0
  y[(56*max_age+1):(57*max_age)] = 0
  y[(57*max_age+1):(58*max_age)] = 0
  y[(58*max_age+1):(59*max_age)] = 0
  y[(59*max_age+1):(60*max_age)] = 0
  
  #Incidence
  
  y[(60*max_age+1):(61*max_age)] = x[(60*max_age+1):(61*max_age)]    
  y[(61*max_age+1):(62*max_age)] = x[(61*max_age+1):(62*max_age)]
  y[(62*max_age+1):(63*max_age)] = x[(62*max_age+1):(63*max_age)]
  y[(63*max_age+1):(64*max_age)] = x[(63*max_age+1):(64*max_age)]
  y[(64*max_age+1):(65*max_age)] = x[(64*max_age+1):(65*max_age)]
  y[(65*max_age+1):(66*max_age)] = x[(65*max_age+1):(66*max_age)]
  y[(66*max_age+1):(67*max_age)] = x[(66*max_age+1):(67*max_age)]
  y[(67*max_age+1):(68*max_age)] = x[(67*max_age+1):(68*max_age)]
  
  
  return(y)
}




catchup_split <- function(x,max_age,cov,target){
  y <- vector("numeric")
  if(target<=4){
    cov_short <- c(0,rep(cov,target), rep(0,(max_age-(target+1))))
    cov_long <- rep(0,max_age)}
  if(target>4){
    cov_short <- c(0,rep(cov,4),rep(0,max_age-5))
    cov_long <- c(0,0,0,0,0,rep(cov,target-4),rep(0,(max_age-5-target+4)))}
  
  #not vaccinated
  
  y[1:max_age] = (1-cov_short-cov_long)*x[1:(max_age)] 
  y[(max_age+1):(2*max_age)] = (1-cov_short-cov_long)*x[(max_age+1):(2*max_age)]
  y[(2*max_age+1):(3*max_age)] = (1-cov_short-cov_long)*x[(2*max_age+1):(3*max_age)]
  y[(3*max_age+1):(4*max_age)] = x[(3*max_age+1):(4*max_age)]
  y[(4*max_age+1):(5*max_age)] =  x[(4*max_age+1):(5*max_age)] 
  y[(5*max_age+1):(6*max_age)] = (1-cov_short-cov_long)*x[(5*max_age+1):(6*max_age)]
  y[(6*max_age+1):(7*max_age)] = (1-cov_short-cov_long)*x[(6*max_age+1):(7*max_age)]
  y[(7*max_age+1):(8*max_age)] = (1-cov_short-cov_long)*x[(7*max_age+1):(8*max_age)]
  y[(8*max_age+1):(9*max_age)] = (1-cov_short-cov_long)*x[(8*max_age+1):(9*max_age)]
  y[(9*max_age+1):(10*max_age)] = x[(9*max_age+1):(10*max_age)]
  y[(10*max_age+1):(11*max_age)] = x[(10*max_age+1):(11*max_age)] 
  y[(11*max_age+1):(12*max_age)] = (1-cov_short-cov_long)*x[(11*max_age+1):(12*max_age)]
  
  #MAV long  
  
  y[(12*max_age+1):(13*max_age)] = cov_long*x[(0*max_age+1):(1*max_age)] + x[(12*max_age+1):(13*max_age)]
  y[(13*max_age+1):(14*max_age)] = cov_long*x[(1*max_age+1):(2*max_age)] + x[(13*max_age+1):(14*max_age)]
  y[(14*max_age+1):(15*max_age)] = cov_long*x[(2*max_age+1):(3*max_age)] + x[(14*max_age+1):(15*max_age)]
  y[(15*max_age+1):(16*max_age)] = x[(15*max_age+1):(16*max_age)]
  y[(16*max_age+1):(17*max_age)] = x[(16*max_age+1):(17*max_age)]
  y[(17*max_age+1):(18*max_age)] = cov_long*x[(5*max_age+1):(6*max_age)] + x[(17*max_age+1):(18*max_age)]
  y[(18*max_age+1):(19*max_age)] = cov_long*x[(6*max_age+1):(7*max_age)] + x[(18*max_age+1):(19*max_age)]
  y[(19*max_age+1):(20*max_age)] = cov_long*x[(7*max_age+1):(8*max_age)] + x[(19*max_age+1):(20*max_age)]
  y[(20*max_age+1):(21*max_age)] = cov_long*x[(8*max_age+1):(9*max_age)] + x[(20*max_age+1):(21*max_age)]
  y[(21*max_age+1):(22*max_age)] = x[(21*max_age+1):(22*max_age)]
  y[(22*max_age+1):(23*max_age)] = x[(22*max_age+1):(23*max_age)]
  y[(23*max_age+1):(24*max_age)] = cov_long*x[(11*max_age+1):(12*max_age)] + x[(23*max_age+1):(24*max_age)]
  
  #MAV short  
  
  y[(24*max_age+1):(25*max_age)] = cov_short*x[(0*max_age+1):(1*max_age)] + x[(24*max_age+1):(25*max_age)]  
  y[(25*max_age+1):(26*max_age)] = cov_short*x[(1*max_age+1):(2*max_age)] + x[(25*max_age+1):(26*max_age)]
  y[(26*max_age+1):(27*max_age)] = cov_short*x[(2*max_age+1):(3*max_age)] + x[(26*max_age+1):(27*max_age)]
  y[(27*max_age+1):(28*max_age)] = x[(27*max_age+1):(28*max_age)]
  y[(28*max_age+1):(29*max_age)] = x[(28*max_age+1):(29*max_age)]
  y[(29*max_age+1):(30*max_age)] = cov_short*x[(5*max_age+1):(6*max_age)] + x[(29*max_age+1):(30*max_age)]
  y[(30*max_age+1):(31*max_age)] = cov_short*x[(6*max_age+1):(7*max_age)] + x[(30*max_age+1):(31*max_age)]
  y[(31*max_age+1):(32*max_age)] = cov_short*x[(7*max_age+1):(8*max_age)] + x[(31*max_age+1):(32*max_age)]
  y[(32*max_age+1):(33*max_age)] = cov_short*x[(8*max_age+1):(9*max_age)] + x[(32*max_age+1):(33*max_age)]
  y[(33*max_age+1):(34*max_age)] = x[(33*max_age+1):(34*max_age)]
  y[(34*max_age+1):(35*max_age)] = x[(34*max_age+1):(35*max_age)]
  y[(35*max_age+1):(36*max_age)] = cov_short*x[(11*max_age+1):(12*max_age)] + x[(35*max_age+1):(36*max_age)]
  
  #Penta long  
  
  y[(36*max_age+1):(37*max_age)] = 0   
  y[(37*max_age+1):(38*max_age)] = 0
  y[(38*max_age+1):(39*max_age)] = 0
  y[(39*max_age+1):(40*max_age)] = 0
  y[(40*max_age+1):(41*max_age)] = 0
  y[(41*max_age+1):(42*max_age)] = 0
  y[(42*max_age+1):(43*max_age)] = 0
  y[(43*max_age+1):(44*max_age)] = 0
  y[(44*max_age+1):(45*max_age)] = 0
  y[(45*max_age+1):(46*max_age)] = 0
  y[(46*max_age+1):(47*max_age)] = 0
  y[(47*max_age+1):(48*max_age)] = 0
  
  #Penta short  
  
  y[(48*max_age+1):(49*max_age)] = 0   
  y[(49*max_age+1):(50*max_age)] = 0
  y[(50*max_age+1):(51*max_age)] = 0
  y[(51*max_age+1):(52*max_age)] = 0
  y[(52*max_age+1):(53*max_age)] = 0
  y[(53*max_age+1):(54*max_age)] = 0
  y[(54*max_age+1):(55*max_age)] = 0
  y[(55*max_age+1):(56*max_age)] = 0
  y[(56*max_age+1):(57*max_age)] = 0
  y[(57*max_age+1):(58*max_age)] = 0
  y[(58*max_age+1):(59*max_age)] = 0
  y[(59*max_age+1):(60*max_age)] = 0
  
  #Incidence
  
  y[(60*max_age+1):(61*max_age)] = x[(60*max_age+1):(61*max_age)]    
  y[(61*max_age+1):(62*max_age)] = x[(61*max_age+1):(62*max_age)]
  y[(62*max_age+1):(63*max_age)] = x[(62*max_age+1):(63*max_age)]
  y[(63*max_age+1):(64*max_age)] = x[(63*max_age+1):(64*max_age)]
  y[(64*max_age+1):(65*max_age)] = x[(64*max_age+1):(65*max_age)]
  y[(65*max_age+1):(66*max_age)] = x[(65*max_age+1):(66*max_age)]
  y[(66*max_age+1):(67*max_age)] = x[(66*max_age+1):(67*max_age)]
  y[(67*max_age+1):(68*max_age)] = x[(67*max_age+1):(68*max_age)]
  
  
  return(y)
}




catchup_split_NGA <- function(x,max_age,cov,target){
  y <- vector("numeric")
  if(target[1]>=2){
    cov_short <- rep(0,max_age)
    cov_long <- c(rep(0,target[1]),rep(cov,target[2]-target[1]+1),rep(0,max_age-target[2]-1))}
  if(target[1]==1 & target[2]>4){
    cov_short <- c(0,rep(cov,4),rep(0,max_age-5))
    cov_long <- c(0,0,0,0,0,rep(cov,target[2]-4),rep(0,(max_age-5-target[2]+4)))}
  
  #not vaccinated
  
  y[1:max_age] = (1-cov_short-cov_long)*x[1:(max_age)] 
  y[(max_age+1):(2*max_age)] = (1-cov_short-cov_long)*x[(max_age+1):(2*max_age)]
  y[(2*max_age+1):(3*max_age)] = (1-cov_short-cov_long)*x[(2*max_age+1):(3*max_age)]
  y[(3*max_age+1):(4*max_age)] = x[(3*max_age+1):(4*max_age)]
  y[(4*max_age+1):(5*max_age)] =  x[(4*max_age+1):(5*max_age)] 
  y[(5*max_age+1):(6*max_age)] = (1-cov_short-cov_long)*x[(5*max_age+1):(6*max_age)]
  y[(6*max_age+1):(7*max_age)] = (1-cov_short-cov_long)*x[(6*max_age+1):(7*max_age)]
  y[(7*max_age+1):(8*max_age)] = (1-cov_short-cov_long)*x[(7*max_age+1):(8*max_age)]
  y[(8*max_age+1):(9*max_age)] = (1-cov_short-cov_long)*x[(8*max_age+1):(9*max_age)]
  y[(9*max_age+1):(10*max_age)] = x[(9*max_age+1):(10*max_age)]
  y[(10*max_age+1):(11*max_age)] = x[(10*max_age+1):(11*max_age)] 
  y[(11*max_age+1):(12*max_age)] = (1-cov_short-cov_long)*x[(11*max_age+1):(12*max_age)]
  
  #MAV long  
  
  y[(12*max_age+1):(13*max_age)] = cov_long*x[(0*max_age+1):(1*max_age)] + x[(12*max_age+1):(13*max_age)]
  y[(13*max_age+1):(14*max_age)] = cov_long*x[(1*max_age+1):(2*max_age)] + x[(13*max_age+1):(14*max_age)]
  y[(14*max_age+1):(15*max_age)] = cov_long*x[(2*max_age+1):(3*max_age)] + x[(14*max_age+1):(15*max_age)]
  y[(15*max_age+1):(16*max_age)] = x[(15*max_age+1):(16*max_age)]
  y[(16*max_age+1):(17*max_age)] = x[(16*max_age+1):(17*max_age)]
  y[(17*max_age+1):(18*max_age)] = cov_long*x[(5*max_age+1):(6*max_age)] + x[(17*max_age+1):(18*max_age)]
  y[(18*max_age+1):(19*max_age)] = cov_long*x[(6*max_age+1):(7*max_age)] + x[(18*max_age+1):(19*max_age)]
  y[(19*max_age+1):(20*max_age)] = cov_long*x[(7*max_age+1):(8*max_age)] + x[(19*max_age+1):(20*max_age)]
  y[(20*max_age+1):(21*max_age)] = cov_long*x[(8*max_age+1):(9*max_age)] + x[(20*max_age+1):(21*max_age)]
  y[(21*max_age+1):(22*max_age)] = x[(21*max_age+1):(22*max_age)]
  y[(22*max_age+1):(23*max_age)] = x[(22*max_age+1):(23*max_age)]
  y[(23*max_age+1):(24*max_age)] = cov_long*x[(11*max_age+1):(12*max_age)] + x[(23*max_age+1):(24*max_age)]
  
  #MAV short  
  
  y[(24*max_age+1):(25*max_age)] = cov_short*x[(0*max_age+1):(1*max_age)] + x[(24*max_age+1):(25*max_age)]  
  y[(25*max_age+1):(26*max_age)] = cov_short*x[(1*max_age+1):(2*max_age)] + x[(25*max_age+1):(26*max_age)]
  y[(26*max_age+1):(27*max_age)] = cov_short*x[(2*max_age+1):(3*max_age)] + x[(26*max_age+1):(27*max_age)]
  y[(27*max_age+1):(28*max_age)] = x[(27*max_age+1):(28*max_age)]
  y[(28*max_age+1):(29*max_age)] = x[(28*max_age+1):(29*max_age)]
  y[(29*max_age+1):(30*max_age)] = cov_short*x[(5*max_age+1):(6*max_age)] + x[(29*max_age+1):(30*max_age)]
  y[(30*max_age+1):(31*max_age)] = cov_short*x[(6*max_age+1):(7*max_age)] + x[(30*max_age+1):(31*max_age)]
  y[(31*max_age+1):(32*max_age)] = cov_short*x[(7*max_age+1):(8*max_age)] + x[(31*max_age+1):(32*max_age)]
  y[(32*max_age+1):(33*max_age)] = cov_short*x[(8*max_age+1):(9*max_age)] + x[(32*max_age+1):(33*max_age)]
  y[(33*max_age+1):(34*max_age)] = x[(33*max_age+1):(34*max_age)]
  y[(34*max_age+1):(35*max_age)] = x[(34*max_age+1):(35*max_age)]
  y[(35*max_age+1):(36*max_age)] = cov_short*x[(11*max_age+1):(12*max_age)] + x[(35*max_age+1):(36*max_age)]
  
  #Penta long  
  
  y[(36*max_age+1):(37*max_age)] = 0   
  y[(37*max_age+1):(38*max_age)] = 0
  y[(38*max_age+1):(39*max_age)] = 0
  y[(39*max_age+1):(40*max_age)] = 0
  y[(40*max_age+1):(41*max_age)] = 0
  y[(41*max_age+1):(42*max_age)] = 0
  y[(42*max_age+1):(43*max_age)] = 0
  y[(43*max_age+1):(44*max_age)] = 0
  y[(44*max_age+1):(45*max_age)] = 0
  y[(45*max_age+1):(46*max_age)] = 0
  y[(46*max_age+1):(47*max_age)] = 0
  y[(47*max_age+1):(48*max_age)] = 0
  
  #Penta short  
  
  y[(48*max_age+1):(49*max_age)] = 0   
  y[(49*max_age+1):(50*max_age)] = 0
  y[(50*max_age+1):(51*max_age)] = 0
  y[(51*max_age+1):(52*max_age)] = 0
  y[(52*max_age+1):(53*max_age)] = 0
  y[(53*max_age+1):(54*max_age)] = 0
  y[(54*max_age+1):(55*max_age)] = 0
  y[(55*max_age+1):(56*max_age)] = 0
  y[(56*max_age+1):(57*max_age)] = 0
  y[(57*max_age+1):(58*max_age)] = 0
  y[(58*max_age+1):(59*max_age)] = 0
  y[(59*max_age+1):(60*max_age)] = 0
  
  #Incidence
  
  y[(60*max_age+1):(61*max_age)] = x[(60*max_age+1):(61*max_age)]    
  y[(61*max_age+1):(62*max_age)] = x[(61*max_age+1):(62*max_age)]
  y[(62*max_age+1):(63*max_age)] = x[(62*max_age+1):(63*max_age)]
  y[(63*max_age+1):(64*max_age)] = x[(63*max_age+1):(64*max_age)]
  y[(64*max_age+1):(65*max_age)] = x[(64*max_age+1):(65*max_age)]
  y[(65*max_age+1):(66*max_age)] = x[(65*max_age+1):(66*max_age)]
  y[(66*max_age+1):(67*max_age)] = x[(66*max_age+1):(67*max_age)]
  y[(67*max_age+1):(68*max_age)] = x[(67*max_age+1):(68*max_age)]
  
  
  return(y)
}






campaign_penta <- function(x,max_age,cov,target_penta){
  y <- vector("numeric")
  if(target_penta[1]==1){
    cov_short <- c(0,rep(cov,4),rep(0,max_age-5))
    cov_long <- c(0,0,0,0,0,rep(cov,target_penta[2]-4),rep(0,(max_age-5-target_penta[2]+4)))
  }
  if(target_penta[1]>1){
    cov_short <- c(0,0,0,0,0,rep(0,max_age-5))
    cov_long <- c(0,0,0,0,0,rep(cov,target_penta[2]-4),rep(0,(max_age-5-target_penta[2]+4)))
  }
  
  if(target_penta[1]==0){
    cov_short <- c(0,0,0,0,0,rep(0,max_age-5))
    cov_long <- c(0,0,0,0,0,rep(0,(max_age-5)))
  }
  
  #not vaccinated
  
  y[1:max_age] = (1-cov_short-cov_long)*x[1:(max_age)] 
  y[(max_age+1):(2*max_age)] = (1-cov_short-cov_long)*x[(max_age+1):(2*max_age)]
  y[(2*max_age+1):(3*max_age)] = (1-cov_short-cov_long)*x[(2*max_age+1):(3*max_age)]
  y[(3*max_age+1):(4*max_age)] = x[(3*max_age+1):(4*max_age)]
  y[(4*max_age+1):(5*max_age)] =  x[(4*max_age+1):(5*max_age)] 
  y[(5*max_age+1):(6*max_age)] = (1-cov_short-cov_long)*x[(5*max_age+1):(6*max_age)]
  y[(6*max_age+1):(7*max_age)] = (1-cov_short-cov_long)*x[(6*max_age+1):(7*max_age)]
  y[(7*max_age+1):(8*max_age)] = (1-cov_short-cov_long)*x[(7*max_age+1):(8*max_age)]
  y[(8*max_age+1):(9*max_age)] = (1-cov_short-cov_long)*x[(8*max_age+1):(9*max_age)]
  y[(9*max_age+1):(10*max_age)] = x[(9*max_age+1):(10*max_age)]
  y[(10*max_age+1):(11*max_age)] = x[(10*max_age+1):(11*max_age)] 
  y[(11*max_age+1):(12*max_age)] = (1-cov_short-cov_long)*x[(11*max_age+1):(12*max_age)]
  
  #MAV long  
  
  y[(12*max_age+1):(13*max_age)] = (1-cov_long)*x[(12*max_age+1):(13*max_age)]
  y[(13*max_age+1):(14*max_age)] = (1-cov_long)*x[(13*max_age+1):(14*max_age)]
  y[(14*max_age+1):(15*max_age)] = (1-cov_long)*x[(14*max_age+1):(15*max_age)]
  y[(15*max_age+1):(16*max_age)] = x[(15*max_age+1):(16*max_age)]
  y[(16*max_age+1):(17*max_age)] = x[(16*max_age+1):(17*max_age)]
  y[(17*max_age+1):(18*max_age)] = (1-cov_long)*x[(17*max_age+1):(18*max_age)]
  y[(18*max_age+1):(19*max_age)] = (1-cov_long)*x[(18*max_age+1):(19*max_age)]
  y[(19*max_age+1):(20*max_age)] = (1-cov_long)*x[(19*max_age+1):(20*max_age)]
  y[(20*max_age+1):(21*max_age)] = (1-cov_long)*x[(20*max_age+1):(21*max_age)]
  y[(21*max_age+1):(22*max_age)] = x[(21*max_age+1):(22*max_age)]
  y[(22*max_age+1):(23*max_age)] = x[(22*max_age+1):(23*max_age)]
  y[(23*max_age+1):(24*max_age)] = (1-cov_long)*x[(23*max_age+1):(24*max_age)]
  
  #MAV short  
  
  y[(24*max_age+1):(25*max_age)] = (1-cov_short)*x[(24*max_age+1):(25*max_age)]   
  y[(25*max_age+1):(26*max_age)] = (1-cov_short)*x[(25*max_age+1):(26*max_age)]
  y[(26*max_age+1):(27*max_age)] = (1-cov_short)*x[(26*max_age+1):(27*max_age)]
  y[(27*max_age+1):(28*max_age)] = x[(27*max_age+1):(28*max_age)]
  y[(28*max_age+1):(29*max_age)] = x[(28*max_age+1):(29*max_age)]
  y[(29*max_age+1):(30*max_age)] = (1-cov_short)*x[(29*max_age+1):(30*max_age)]
  y[(30*max_age+1):(31*max_age)] = (1-cov_short)*x[(30*max_age+1):(31*max_age)]
  y[(31*max_age+1):(32*max_age)] = (1-cov_short)*x[(31*max_age+1):(32*max_age)]
  y[(32*max_age+1):(33*max_age)] = (1-cov_short)*x[(32*max_age+1):(33*max_age)]
  y[(33*max_age+1):(34*max_age)] = x[(33*max_age+1):(34*max_age)]
  y[(34*max_age+1):(35*max_age)] = x[(34*max_age+1):(35*max_age)]
  y[(35*max_age+1):(36*max_age)] = (1-cov_short)*x[(35*max_age+1):(36*max_age)]
  
  #Penta long  
  
  y[(36*max_age+1):(37*max_age)] = cov_long*(x[1:(max_age)]+x[(12*max_age+1):(13*max_age)]) + x[(36*max_age+1):(37*max_age)]   
  y[(37*max_age+1):(38*max_age)] = cov_long*(x[(max_age+1):(2*max_age)]+x[(13*max_age+1):(14*max_age)]) + x[(37*max_age+1):(38*max_age)]
  y[(38*max_age+1):(39*max_age)] = cov_long*(x[(2*max_age+1):(3*max_age)]+x[(14*max_age+1):(15*max_age)]) + x[(38*max_age+1):(39*max_age)]
  y[(39*max_age+1):(40*max_age)] = x[(39*max_age+1):(40*max_age)]
  y[(40*max_age+1):(41*max_age)] = x[(40*max_age+1):(41*max_age)]
  y[(41*max_age+1):(42*max_age)] = cov_long*(x[(5*max_age+1):(6*max_age)]+x[(17*max_age+1):(18*max_age)]) + x[(41*max_age+1):(42*max_age)]
  y[(42*max_age+1):(43*max_age)] = cov_long*(x[(6*max_age+1):(7*max_age)]+x[(18*max_age+1):(19*max_age)]) + x[(42*max_age+1):(43*max_age)]
  y[(43*max_age+1):(44*max_age)] = cov_long*(x[(7*max_age+1):(8*max_age)]+x[(19*max_age+1):(20*max_age)]) + x[(43*max_age+1):(44*max_age)]
  y[(44*max_age+1):(45*max_age)] = cov_long*(x[(8*max_age+1):(9*max_age)]+x[(20*max_age+1):(21*max_age)]) + x[(44*max_age+1):(45*max_age)]
  y[(45*max_age+1):(46*max_age)] = x[(45*max_age+1):(46*max_age)]
  y[(46*max_age+1):(47*max_age)] = x[(46*max_age+1):(47*max_age)]
  y[(47*max_age+1):(48*max_age)] = cov_long*(x[(11*max_age+1):(12*max_age)]+x[(23*max_age+1):(24*max_age)]) + x[(47*max_age+1):(48*max_age)]
  
  #Penta short  
  
  y[(48*max_age+1):(49*max_age)] = cov_short*(x[1:(max_age)]+x[(24*max_age+1):(25*max_age)]) + x[(48*max_age+1):(49*max_age)]   
  y[(49*max_age+1):(50*max_age)] = cov_short*(x[(1*max_age+1):(2*max_age)]+x[(25*max_age+1):(26*max_age)]) + x[(49*max_age+1):(50*max_age)]
  y[(50*max_age+1):(51*max_age)] = cov_short*(x[(2*max_age+1):(3*max_age)]+x[(26*max_age+1):(27*max_age)]) + x[(50*max_age+1):(51*max_age)]
  y[(51*max_age+1):(52*max_age)] = x[(51*max_age+1):(52*max_age)]
  y[(52*max_age+1):(53*max_age)] = x[(52*max_age+1):(53*max_age)]
  y[(53*max_age+1):(54*max_age)] = cov_short*(x[(5*max_age+1):(6*max_age)]+x[(29*max_age+1):(30*max_age)]) + x[(53*max_age+1):(54*max_age)]
  y[(54*max_age+1):(55*max_age)] = cov_short*(x[(6*max_age+1):(7*max_age)]+x[(30*max_age+1):(31*max_age)]) + x[(54*max_age+1):(55*max_age)]
  y[(55*max_age+1):(56*max_age)] = cov_short*(x[(7*max_age+1):(8*max_age)]+x[(31*max_age+1):(32*max_age)]) + x[(55*max_age+1):(56*max_age)]
  y[(56*max_age+1):(57*max_age)] = cov_short*(x[(8*max_age+1):(9*max_age)]+x[(32*max_age+1):(33*max_age)]) + x[(56*max_age+1):(57*max_age)]
  y[(57*max_age+1):(58*max_age)] = x[(57*max_age+1):(58*max_age)]
  y[(58*max_age+1):(59*max_age)] = x[(58*max_age+1):(59*max_age)]
  y[(59*max_age+1):(60*max_age)] = cov_short*(x[(11*max_age+1):(12*max_age)]+x[(35*max_age+1):(36*max_age)]) + x[(59*max_age+1):(60*max_age)]
  
  #Incidence
  
  y[(60*max_age+1):(61*max_age)] = x[(60*max_age+1):(61*max_age)]    
  y[(61*max_age+1):(62*max_age)] = x[(61*max_age+1):(62*max_age)]
  y[(62*max_age+1):(63*max_age)] = x[(62*max_age+1):(63*max_age)]
  y[(63*max_age+1):(64*max_age)] = x[(63*max_age+1):(64*max_age)]
  y[(64*max_age+1):(65*max_age)] = x[(64*max_age+1):(65*max_age)]
  y[(65*max_age+1):(66*max_age)] = x[(65*max_age+1):(66*max_age)]
  y[(66*max_age+1):(67*max_age)] = x[(66*max_age+1):(67*max_age)]
  y[(67*max_age+1):(68*max_age)] = x[(67*max_age+1):(68*max_age)]
  
  
  return(y)
}





age_dstr <- function(x){
  year <- floor(x[,1]/365)
  tmp <- numeric(12)
  dat2 <- x[,-1]
  
  N <- matrix(0,length(year),max_age)
  
  for(ak in 1:max_age){
    # tmp[1] = index2(i = 1,a = ak, y = max_age)
    for(i in 1:12){
      tmp[i] <- index2(i,ak,max_age)
    }
    N[,ak] <- apply(dat2[,tmp],1,sum)
    
  }
  
  N = cbind(year,N)
  N1 = matrix(0,length(unique(year))-1,max_age+1)
  # N1 = matrix(0,length(unique(year)),max_age+1)
  
  for (i in 1:(length(unique(year))-1)){
    #   for (i in 1:(length(unique(year)))){
    
    N1[i,] = N[year==unique(year)[i],][1,]
  }
  
  # N1 = cbind(unique(year),N1)
  return(N1)
}



create_template <- function(x){
  # 
  
  # 
  mA <- cbind(x[,1]+2000,x[,2:101]*sum(N)*risk,x[,202])
  mCWYX <- cbind(x[,1]+2000,x[,102:201]*sum(N)*risk,x[,202])
  
  colnames(mA) <- c("time",paste(0:99, sep = ""),"sim")
  colnames(mCWYX) <- c("time",paste(0:99, sep = ""),"sim")
  
  mA <- data.frame(mA)
  mCWYX <- data.frame(mCWYX)
  
  mA2 <- melt(mA,id.vars = c("time","sim"))
  colnames(mA2) <- c("year","sim","age","casesA")
  mA2$age<-gsub("X","",as.character(mA2$age))
  
  mCWYX2 <- melt(mCWYX,id.vars = c("time","sim"))
  colnames(mCWYX2) <- c("year","sim","age","casesCWYX")
  mCWYX2$age<-gsub("X","",as.character(mA2$age))
  
  allcases <- data.frame(country= rep(country,length(mA2$year)),
                         year = mA2$year,
                         sim=mA2$sim,
                         age=mA2$age, 
                         casesA= mA2$casesA,
                         casesCWYX = mCWYX2$casesCWYX,
                         deathsA = 0.1*mA2$casesA,
                         deathsCWYX = 0.1*mCWYX2$casesCWYX)
  
  allcases <- allcases[allcases$year<=2100,]
  
  allcases2 <- allcases[order(allcases$year,allcases$sim),]
  
  
  LifeExp <-read.csv("C:/Users/ak889/Dropbox (Personal)/WORK/Penta_model/VIMC_runs_2024/Data/202310gavi-6_dds-202208_life_ex_both.csv")
  
  LifeExp = LifeExp[LifeExp$year>=year,]
  
  agb = unique(LifeExp$age_from)
  agb_life <- agb[-1]
  fin_exp = vector()
  life <- matrix(0,nrow = max_age, length(unique(LifeExp$year)))
  for(i in 1:length(unique(LifeExp$year))){
    expect = LifeExp[LifeExp$country_code==country & LifeExp$year==unique(LifeExp$year)[i],"value"]
    life[,i] = ExpandInvasion(expect,agb_life,max_age)
    fin_exp = c(fin_exp,rep(life[,i],5))
  }
  fin_exp = c(fin_exp,fin_exp[9901:10000])
  
  
  nsims <- length(unique(allcases2$sim))
  
  fin_exp_all <- rep(fin_exp,nsims)
  
  tmp <- list()
  for(i in 1:length(unique(allcases$sim))){
    tmp[[i]] <- allcases2[allcases2$sim==unique(allcases2$sim)[i],]
  }
  
  allcases2 <- do.call("rbind",tmp)
  
  
  allcases2$expect <- fin_exp_all
  #  allcases2$deathsA <- allcases2$casesA*0.1
  
  
  allcases2$yldA <- allcases2$deathsA*allcases2$expect   
  allcases2$yllA <- (allcases2$casesA-allcases2$deathsA)*0.072*0.26*allcases2$expect
  allcases2$dalysA <- allcases2$yldA+allcases2$yllA
  
  # allcases2$expect <- fin_exp_all
  # allcases2$deathsCWYX <- allcases2$casesCWYX*0.1
  
  
  allcases2$yldCWYX <- allcases2$deathsCWYX*allcases2$expect   
  allcases2$yllCWYX <- (allcases2$casesCWYX-allcases2$deathsCWYX)*0.072*0.26*allcases2$expect
  allcases2$dalysCWYX <- allcases2$yldCWYX+allcases2$yllCWYX
  
  
  
  
  newtest2 <- allcases2[order(allcases2$sim),]
  
  
  #############################################################################################################
  
  
  foo <- data.frame(disease = rep("MenA", nrow(newtest2)),
                    run_id = newtest2$sim,
                    year = newtest2$year,
                    age = newtest2$age,
                    country = rep(country,nrow(newtest2)),
                    country_name = rep(country_name,nrow(newtest2)),
                    cohort_size = rep(0,nrow(newtest2)),
                    yll_A = newtest2$yllA, 
                    yll_CWYX = newtest2$yllCWYX,
                    cases = newtest2$casesA,
                    cases_CWYX = newtest2$casesCWYX,
                    dalys = newtest2$dalysA,
                    dalys_CWYX = newtest2$dalysCWYX,
                    deaths = newtest2$deathsA,
                    deaths_CWYX = newtest2$deathsCWYX)
  
  return(foo)
  
}
