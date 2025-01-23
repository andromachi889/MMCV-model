
source("ModelFunctions_2.R")

  # country = "Burkina Faso"
country = country

birth <- read.csv("202110gavi-2_dds-201910_cbr_both.csv")
mortality <-read.csv("202110gavi-2_dds-201910_2_mort_rate_both.csv")
population <- read.csv("202110gavi-2_dds-201910_2_int_pop_both.csv")





# tmpdat <- data[data$Country==country,]
N <- population[population$country_code==country,"X2000"]
N <- N[1:100]
InitPop = N/sum(N)
# N <- N*risk
B <- birth[birth$country_code==country & birth$year==year, "value"]/365


max_age=100
a = 1/365

mu <- numeric(max_age)
mu[1] <- B/InitPop[1]-a
for(i in 2:(max_age-1)){
  mu[i] = a*(InitPop[i-1]/InitPop[i]-1)
}
mu[max_age] = a*(InitPop[max_age-1]/InitPop[max_age])

ageGroupBreaks = c(1,seq(5,100,5))

L= matrix(1,length(ageGroupBreaks),length(ageGroupBreaks))
L[3,3:4] = 4
L[4,3:5] = 4
L[5,4] = 4
L[5,5:6] = 6
L[6,5:7] = 6
L[7,6:7] = 6
L[7,8] = 2
for(i in 8:(length(ageGroupBreaks)-1)){
  L[i,(i-1):(i+1)] = 2
}
L[length(ageGroupBreaks),(length(ageGroupBreaks)-1):length(ageGroupBreaks)]=2
L=L[2:21,2:21]

# L=L-0.25*L
L = L+classify*L

ir2 = c(3.25, 3, 2, 1.75, 1.5, 1,0.75,0.5,0.25,0.05)
AgB = c(1,5,10,15,20,25,30,40,60,100)
ir_A <- ExpandInvasion(ir2,AgB,max_age)*0.0478/365
ir_R <- ExpandInvasion(ir2,AgB,max_age)*0.0478/365
 ir_R <- 0.2*ir_R

#mortality

mortality = mortality[mortality$year>=year,]

agb = unique(mortality$age_from)
agb_mort <- agb[-1]

fin_mu = vector()
mort <- matrix(0,nrow = max_age, length(unique(mortality$year)))
for(i in 1:length(unique(mortality$year))){
  ratemort = mortality[mortality$country_code==country & mortality$year==unique(mortality$year)[i],"value"]
  mort[,i] = ExpandInvasion(ratemort,agb_mort,max_age)
  mort[100,i] = 1
  fin_mu = c(fin_mu,rep(mort[,i],5))
}
fin_mu = c(fin_mu,fin_mu[8901:9000])

fin_mu <- fin_mu/365


##Birth rate
birth <- birth[birth$country_code ==country &birth$year>=2000,]
birth_rate <- birth[,'value']/365


#parameters
beta_A = 12.5/365
beta_R = (12.5/365)*0.4
phi = 0.0839/365
rho = 52/365
alpha_A = 12/365
alpha_R = 6/365
epsbeta = 0
delta = 0.9
ksi = 0.9


ws = 0.1/365    # 10 years
wl = 0.1/365

delta_R <- 0.6
ksi_R <- 0.9


params = c(phi=phi, rho=rho, alpha_A=alpha_A, alpha_R=alpha_R, epsbeta=epsbeta, beta_A=beta_A, beta_R=beta_R, ksi=ksi, delta=delta,ws = ws, wl=wl, delta_R=delta_R, ksi_R=ksi_R)


Cm = Expand(L,ageGroupBreaks,max_age)
addto <- as.numeric(rep(L[1,1],100))
Cm2 <- cbind(addto,Cm)
addto2 <- as.numeric(rep(L[1,1],101))
Cm2 <- rbind(addto2,Cm2)
Cm2 <- unname(Cm2, force = FALSE)


Cm <- Cm2[1:100,1:100]


epi_A <- rep(0,100)
epi_R <- rep(0,100)
