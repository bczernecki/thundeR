#'
#' ECAPE calculation procedure implemented with support of Dr John M Peters
#' 

#' CONSTANTS
Rd=287.04 #dry gas constant
Rv=461.5 #water vapor gas constant
epsilon=Rd/Rv
cp=1005 #specific heat of dry air at constant pressure
g=9.81 #gravitational acceleration
xlv=2501000 #reference latent heat of vaporization at the triple point temperature
xls=2834000 #reference latent heat of sublimation at the triple point temperature
cpv=1870 #specific heat of water vapor at constant pressure
cpl=4190 #specific heat of liquid water
cpi=2106 #specific heat of ice
ttrip=273.15 #triple point temperature
eref=611.2 #reference pressure at the triple point temperature
gamma=Rd/cp #POTENTIAL TEMPERATURE EXPONENT
Gamma_d=g/cp #DRY ADIABATIC LAPSE RATE
pref=611.65 #REFERENCE VAPOR PRESSURE OF WATER VAPOR AT TRIPLE POINT TEMPERATURE
T1=273.15
T2=253.15
L=250 #MIXING DEPTH
sigma=1.1
alpha=0.8
Pr=1/3 #PRANDTL NUMBER
ksq=0.18 #VON KARMAN CONSTANT

#'
#' FUNCTIONS
#'

heaviside <- function(x, a = 0){   
  result = (sign(x-a) + 1)/2
  return(result)
}

dpt2w <- function(p,t,dpt){
  RH = aiRthermo::dewpointdepression2rh(p*100, t+273.15, t-dpt, consts = export_constants())
  RH = ifelse(RH>100,100,RH)
  RH = ifelse(RH<0,0,RH)
  W = aiRthermo::rh2w(p*100, t+273.15, RH, consts = export_constants())
  return(W)
}


compute_rsat <- function(T,p,iceflag){
  omeg = ((T - T1)/(T2-T1))*heaviside((T - T1)/(T2-T1))*heaviside((1 - (T - T1)/(T2-T1))) + heaviside(-(1 - (T - T1)/(T2-T1)))
  if(iceflag==0){
    term1=(cpv-cpl)/Rv
    term2=(xlv-ttrip*(cpv-cpl))/Rv
    esl=exp((T-ttrip)*term2/(T*ttrip))*eref*(T/ttrip)^(term1)
    qsat=epsilon*esl/(p-esl)         
  }
  if(iceflag==1){
    term1=(cpv-cpl)/Rv
    term2=(xlv-ttrip*(cpv-cpl))/Rv
    esl_l=exp((T-ttrip)*term2/(T*ttrip))*eref*(T/ttrip)^(term1)
    qsat_l=epsilon*esl_l/(p-esl_l)
    term1=(cpv-cpi)/Rv
    term2=(xls-ttrip*(cpv-cpi))/Rv
    esl_i=exp((T-ttrip)*term2/(T*ttrip))*eref*(T/ttrip)^(term1)
    qsat_i=epsilon*esl_i/(p-esl_i)
    qsat=(1-omeg)*qsat_l + (omeg)*qsat_i
  }
  if(iceflag==2){
    term1=(cpv-cpi)/Rv
    term2=( xls-ttrip*(cpv-cpi))/Rv
    esl=exp((T-ttrip)*term2/(T*ttrip))*eref*(T/ttrip)^(term1)
    qsat=epsilon*esl/(p-esl)
  }
  return(qsat)
}

MSE0_F <- function(t0,q0,z0){
  MSE0 = cp * t0 + xlv * q0 + g * z0
  return(MSE0)
}

MSE0_star_F <- function(t0,p0,z0){
  rsat = compute_rsat(t0,p0,0)
  qsat = (1 - rsat) * rsat
  MSE0_star = cp * t0 + xlv * qsat + g * z0
  return(MSE0_star)
}

MSE0_bar_F <- function(MSE0){
  MSE0_bar=rep(0,length(MSE0))
  for(iz in 1:length(MSE0_bar)){
    MSE0_bar[iz]=mean(MSE0[1:iz])
  }
  return(MSE0_bar)
}

compute_NCAPE <- function(p0,z0,t0,q0,LFC,EL){
  MSE0 = MSE0_F(t0,q0,z0)
  MSE0_star = MSE0_star_F(t0,p0,z0)
  MSE0_bar = MSE0_bar_F(MSE0)
  int_arg = - ( g/(cp*t0) )*( MSE0_bar - MSE0_star)
  if((LFC + EL) > 0){
  ind_LFC = which(abs(z0 - LFC) == min(abs(z0 - LFC)))
  ind_EL = which(abs(z0 - EL) == min(abs(z0 - EL)))
  NCAPE = mean(int_arg[ind_LFC:ind_EL]) * (EL-LFC)
  NCAPE = ifelse(NCAPE<0,0,NCAPE)
  } else {
    NCAPE = 0
  }
  return (NCAPE)
}

compute_ETILDE <- function(p0,z0,q0,t0,LFC,EL,CAPE,V_SR){
  l=L/EL
  pitchfork=ksq*(alpha**2)*(3.14159265359**2)*L/(4*Pr*(sigma**2)*EL)
  vsr_tilde = V_SR/sqrt(2*CAPE)
  NCAPE = compute_NCAPE(p0,z0,t0,q0,LFC,EL)
  N_tilde = NCAPE/CAPE
  E_tilde = vsr_tilde**2 + ( -1 - pitchfork - (pitchfork/(vsr_tilde**2 ))*N_tilde + sqrt((1 + pitchfork + (pitchfork/(vsr_tilde**2 ))*N_tilde)**2 + (4*(pitchfork/(vsr_tilde**2 ))*(1 - pitchfork*N_tilde) ) ) )/( 2*pitchfork/(vsr_tilde**2) )
  eps = 2*ksq*L/(EL*Pr)
  varepsilon = 0.65*eps*(alpha**2)*(3.14159265359**2)*E_tilde/(4*(sigma**2)*EL*(vsr_tilde**2 ) ) #THIS IS THE FRACTIONAL ENTRAINMENT RATE
  Radius = sqrt(2*ksq*L/(Pr*varepsilon)) # UPDRAFT RADIUS
  return(c(E_tilde,varepsilon,Radius))
}

ECAPE_parcel <- function(p0,z0,t0,q0,LFC,EL,CAPE,V_SR,Tv,Tp){
  MSE0 = MSE0_F(t0,q0,z0)
  MSE0_star <- MSE0_star_F(t0,p0,z0)
  MSE0_bar <- MSE0_bar_F(MSE0)
  NCAPE <- compute_NCAPE(p0,z0,t0,q0,LFC,EL)
  E_tilde <- compute_ETILDE(p0,z0,q0,t0,LFC,EL,CAPE,V_SR)[1]
  Buoy_UD <- -g*(Tv-Tp)/(Tp)
  eps = 2*ksq*L/(EL*Pr)
  vsr_tilde = V_SR/sqrt(2*CAPE)
  varepsilon = 0.65*eps*(alpha^2)*(pi^2)*E_tilde/(4*(sigma^2)*EL*(vsr_tilde^2 ))
  B_ent = Buoy_UD*exp(-varepsilon*z0) + (g/(cp*t0))*(1 - exp(-varepsilon*z0) )*(MSE0_bar-MSE0_star) 
  Te = -(g*Tv)/(B_ent-g)
  return(Te-273.15)
}

ECAPE_value <- function(p0,z0,t0,q0,LFC,EL,CAPE,V_SR){
  ECAPE = compute_ETILDE(p0,z0,q0,t0,LFC,EL,CAPE,V_SR)[1]*CAPE
  return(ECAPE)
}

ENTRAINMENT_rate <- function(p0,z0,t0,q0,LFC,EL,CAPE,V_SR){
  ENTRAINMENT = compute_ETILDE(p0,z0,q0,t0,LFC,EL,CAPE,V_SR)[2]
  return(ENTRAINMENT)
}

UPDRAFT_width <- function(p0,z0,t0,q0,LFC,EL,CAPE,V_SR){
  UPDRAFT = compute_ETILDE(p0,z0,q0,t0,LFC,EL,CAPE,V_SR)[3]
  return(UPDRAFT)
}

#'
#' 
#'
