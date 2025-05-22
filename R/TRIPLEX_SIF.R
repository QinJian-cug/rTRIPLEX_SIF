library(pracma)
RC_PSII <- read.csv('./data/RC_PSII.csv') # Read the standard database of chlorophyll fluorescence by relative path
PPFD_Chlf <- function(wavelength, radiance)
{
  # Define the constant of convert-function
  h <- 6.62e-34 #
  c <- 3e8 #
  Na <- 6.02e23 #

  #
  rad_umol <- (pi * radiance * wavelength * 1e6 ) / (h * c * Na * 1e3 * 1e9)
  # PAR_ppfd <- sum(rad_umol, na.rm = TRUE)

  rad_umol_temp <- rep(0, length = length(radiance))
  for (i in 1: (length(radiance) - 1))
    {
    delta_wl <- wavelength[i + 1] - wavelength[i]
    rad_umol_temp[i] <- (rad_umol[i] + rad_umol[i + 1]) / 2 * delta_wl
    }

  PAR_wl <- which(wavelength > 640 & wavelength < 850)
  # 在选定的波长范围内对光子通量密度求和 (忽略 NA 值)
  PAR_ppfd <- sum(rad_umol_temp[PAR_wl], na.rm = TRUE) # µmol/m²/s

  return(PAR_ppfd)
}

#' rTRIPLEX_SIF
#'
#' @param Input_variable A table as described in \code{\link{Inputpara}} containing the information about input variables.
#' @param Input_parameter A table as described in \code{\link{Inputvariable}} containing the information about input parameters.
#' @returns A list with class "result". More details on the output is \code{\link{result}}
#' @export
#'
#' @examples
#' library(rTRIPLEXSIF)
#' data("onemonth_exam")
#' data("Inputpara")
#' out<- TRIPLEX_SIF (Input_variable = onemonth_exam,
#'                    Input_parameter = Inputpara)

TRIPLEX_SIF<- function(Input_variable,Input_parameter)
{

  ## To maintain user's original options
  #oldpar <- par(no.readonly = TRUE)
  #on.exit(par(oldpar))

  # Measured net ecosystem production and evapotranspiration
  Input_variable$GPP<-Input_variable$GPPob*1800/1000/44*12
  Input_variable$OETS<-0.43*Input_variable$LE/(597-0.564*Input_variable$Ta)

  # Input variable assignment
  Cappm<-Input_variable$Cof # CO2 concentration in the atmosphere, ppm
  Tem<-Input_variable$Ta # Air temperature, degrees Celsius
  Day<-Input_variable$DOY # Day of year, day
  time<-Input_variable$time # Local time, h
  RH<-Input_variable$RH # Relative humidity, %
  PPFD<-Input_variable$PPFD # Photosynthetic photon flux density, μmol m-2 s-1
  SWC<-Input_variable$SVWC30cm # Soil volumetric moisture at depth of 30 cm, %
  Rn<-Input_variable$Rn # Net radiation at the canopy surface, W m-2
  G<-Input_variable$G # Soil heat flux, W m-2
  VPD1<-Input_variable$VPDhpa/10 # Vapor pressure deficit, kPa

  # The new variable appended of SIF
  NIRv<- Input_variable$NIRv # Normalized difference vegetation index
  fPAR<- Input_variable$fPAR # The ratio of absorbed PAR by canopy
  SIFob<- Input_variable$SIF # Solar induced fluorescence, mW m-2 nm-1 sr-1

  # Procedure variable assignment for loop computation，初始化归零
  x<-nrow(Input_variable)
  Cippm<-rep(0,x);Ca<-rep(0,x);Ci<-rep(0,x);Error1<-rep(0,x);COp<-rep(0,x);
  Kc<-rep(0,x);Ko<-rep(0,x);K<-rep(0,x);fN<-rep(0,x);fT<-rep(0,x);Vm<-rep(0,x);
  Vc<-rep(0,x);Jmax<-rep(0,x); Vj1<-rep(0,x);Rd<-rep(0,x);
  Vjsh1<-rep(0,x);Ld<-rep(0,x);Et<-rep(0,x);to<-rep(0,x);h<-rep(0,x);d<-rep(0,x);
  sinB<-rep(0,x);kb1<-rep(0,x);kb2<-rep(0,x);pcb<-rep(0,x);X<-rep(0,x);Y<-rep(0,x);
  Z<-rep(0,x);Icsh<-rep(0,x);bsh<-rep(0,x);csh<-rep(0,x);Jsh1<-rep(0,x);
  Jsh2<-rep(0,x);bsunsh<-rep(0,x);csunsh<-rep(0,x);Jsunsh1<-rep(0,x);
  Jsunsh2<-rep(0,x);Vjsh2<-rep(0,x);Vjsh<-rep(0,x);Vjsun<-rep(0,x);VcVjsun<-rep(0,x);
  VcVjsh<-rep(0,x);Asunsh<-rep(0,x);
  Acsun<-rep(0,x);Acsh<-rep(0,x);Acanopy<-rep(0,x);AV<-rep(0,x);Lsh<-rep(0,x)
  Vj<-rep(0,x);kb<-rep(0,x);Lsun<-rep(0,x);ZZ<-rep(0,x);
  Icsun<-rep(0,x);Ic<-rep(0,x);Ratio1<-rep(0,x)
  Asun<-rep(0,x);Ash<-rep(0,x);Ratio<-rep(0,x);CJ4<-rep(0,x);fsun<-rep(0,x);
  fsh<-rep(0,x);GPPsec<-rep(0,x);Re30min<-rep(0,x);Rh30min<-rep(0,x);
  NEP30min<-rep(0,x);gs<-rep(0,x);Rmf<-rep(0,x);Rms<-rep(0,x);Rmr<-rep(0,x);
  Rm<-rep(0,x);Rgf<-rep(0,x);Rgs<-rep(0,x);Rgr<-rep(0,x);Rg<-rep(0,x);
  Rm30min<-rep(0,x);Rg30min<-rep(0,x);NPP<-rep(0,x);
  NPP30min<-rep(0,x);GPP30min<-rep(0,x);Z0M<-rep(0,x);Z0V<-rep(0,x);
  d<-rep(0,x);svt<-rep(0,x);ra<-rep(0,x);LES<-rep(0,x);gsms<-rep(0,x);ETS<-rep(0,x);
  fesc <- rep(0,x); SIFps2 <- rep(0,x); SIFps <- rep(0,x); SIFfc <- rep(0,x);
  qL <- rep(0,x); PSmax <- rep(0, x);
  SIFfc_vector <- matrix(0, nrow = 211, ncol = x)
  SIFfull <<- rep(0, x); J <<- rep(0, x);

  # Calculation process of TRIPLEX-CW-Flux model
  for(i in 1:x){
    # Rubisco-limited gross photosynthesis rates. (Assum: Rubisco-limited gross photosynthesis rates in shaded leaves and sunlit leaves are equal)
    j<-0
    # Setting initial value of intercellular CO2 concentration, assume: Cippm=Cappm
    Cippm[i] <- Cappm[i] - j

    Ca[i]<-Cappm[i]/1000000*1.013*100000 # Unit conversion: from ppm to Pa
    Ci[i]<-Cippm[i]/1000000*1.013*100000 # Unit conversion: from ppm to Pa

    Vm25<-Input_parameter$Vm25 # Maximum carboxylation rate at 25℃, μmol m-2 s-1
    Rgas<-Input_parameter$Rgas # Molar gas constant, m3 Pa mol-1 K-1
    O2<-Input_parameter$O2 # Oxygen concentration in the atmosphere, Pa
    N<-Input_parameter$N # Leaf nitrogen content, %
    Nm<-Input_parameter$Nm # Maximum nitrogen content, %
    fN<-N/Nm # Nitrogen limitation term
    Error1[i]<-Ci[i]/Ca[i] # Just for check

    COp[i]<-1.92*10^-4*O2*1.75^((Tem[i]-25)/10) # CO2 concentration point without dark respiration, Pa
    Kc[i]<-30*2.1^((Tem[i]-25)/10) # Michaelis–Menten constants for CO2, Pa
    Ko[i]<-30000*1.2^((Tem[i]-25)/10) # Michaelis–Menten constants for O2, Pa
    K[i]<-Kc[i]*(1+O2/Ko[i]) # Function of enzyme kinetics, Pa
    fT[i]<-(1+exp((-220000+710*(Tem[i]+273))/(Rgas*(Tem[i]+273))))^-1 # Temperature limitation term
    Vm[i]<-Vm25*2.4^((Tem[i]-25)/10)*fT[i]*fN # Maximum carboxylation rate, μmol m-2 s-1
    Vc[i]<-Vm[i]*(Ci[i]-COp[i])/(Ci[i]+K[i]) # Rubisco-limited gross photosynthesis rates, μmol m-2 s-1

    # SIFob convert to SIFfull according to Liu(2022)
    fPSII <- Input_parameter$fPSII # The contribution of PSII fluorescence to SIF
    # fPAR[i] <- 0.79 * (0.516 * WDRVI[i] + 0.726) # The radio of absorbed photosynthesis active radiation of canopy
    fesc[i] <- NIRv[i] / fPAR[i] # The escaped probability of SIF photon from leaf level to canopy level
    svd_result <- svd(RC_PSII) # Using SVD method to decompose the standard database of chlorophyll fluorescence
    fc <- svd_result$v[,1] # Function can convert SIF from PSII wavelength to total wavelength by using SVD

    SIFps2[i] <- SIFob[i] * fPSII # Convert observed SIF to SIF in PSII
    SIFps[i] <- SIFps2[i] / (0.9 * fesc[i]) # Downscale SIF from canopy scale to photosystem scale
    SIFfc_vector[,i] <- SIFps[i] / fc[760-640+1] * fc # Upscale the single point SIF to full wavelength SIF by using vector calculate
    SIFfull[i] <- PPFD_Chlf(640 : 850, SIFfc_vector[,i]) # Integrate the full wavelength SIF and convert unit from mW/m2/nm/sr to umol/m2/s

    # Light-limited gross photosynthesis rates for big leaf
    Kdf <- Input_parameter$Kdf # The radio of kd and kf, kd is how fast the heat loss happens, and kf is the rate of fluorscence emission
    qL[i] <- exp(-0.001 * PPFD[i]) * 1 #The redox status of PSII reaction centers
    PSmax[i] <- Tem[i]^2 * (-0.0011) + Tem[i] * 0.036 + 0.44 #The maximum photochemical efficiency of PSII
    J[i] <- SIFfull[i] * qL[i] * ((1 + Kdf) * PSmax[i]) / (1 - PSmax[i])

    #Jmax[i]<-29.1+1.64*Vm[i] # Light-saturated rate of electron transport in the photosynthetic carbon reduction cycle in leaf cells, μmol m-2 s-1
    #J[i]<-Jmax[i]*PPFD[i]/(PPFD[i]+2.1*Jmax[i]) # Electron transport rate, μmol m-2 s-1

    Vj1[i]<-J[i]*(Ci[i]-COp[i])/(4.5*Ci[i]+10.5*COp[i])
    Vj[i]<-ifelse(Vj1[i]>0,Vj1[i],0) # Light-limited gross photosynthesis rates, μmol m-2 s-1
    # Vj[i]<-max(J[i]*(Ci[i]-COp[i])/(4.5*Ci[i]+10.5*COp[i]),0) # Light-limited gross photosynthesis rates, μmol m-2 s-1

    # Solar geometry
    Ls<-Input_parameter$Ls # Standard longitude of time zone
    Le<-Input_parameter$Le # Local longitude, degree
    latitude<-Input_parameter$latitude # Local latitude, degree
    lat<-latitude*3.14/180 # Convert angle to radian, radian
    LAI<-Input_parameter$LAI # Leaf area index, m2 m-2
    Ld[i]<-2*3.14*(Day[i]-1)/365 # Day angle
    Et[i]<-0.017+0.4281*cos(Ld[i])-7.351*sin(Ld[i])-3.349*cos(2*Ld[i])-9.731*sin(Ld[i]) # Equation of time, min
    to[i]<-12+(4*(Ls-Le)-Et[i])/60 # Solar noon, h
    h[i]<-3.14*(time[i]-to[i])/12 # Hour angle of sun, radians
    d[i]<--23.4*3.14/180*cos(2*3.14*(Day[i]+10)/365) # Solar declination angle, radians
    sinB[i]<-sin(lat)*sin(d[i])+cos(lat)*cos(d[i])*cos(h[i]) # Solar elevation angle, radians

    # Sun/Shade model
    # Leaf area index for sunlit and shaded of the canopy
    kb[i]<-ifelse(sinB[i]>0,0.5/sinB[i],0) # Beam radiation extinction coefficient of canopy
    LAIc<-LAI # Leaf area index of canopy, m2 m-2
    Lsun[i]<-ifelse(kb[i]==0,0,(1-exp(-kb[i]*LAIc))/kb[i]) # Sunlit leaf area index of canopy, m2 m-2
    Lsh[i]<-LAIc-Lsun[i] # Shaded leaf area index of canopy, m2 m-2

    # Canopy reflection coefficients
    lsc<-0.15 # Leaf scattering coefficient of PAR
    ph<-(1-sqrt(1-lsc))/(1+sqrt(1-lsc)) # Reflection coefficient of a canopy with horizontal leaves
    pcd<-0.036 # Canopy reflection coefficient for diffuse PAR
    I<-2083 # Total incident PAR, μmol m-2 ground s-1
    fd<-0.159 # Fraction of diffuse irradiance
    Ib0<-I*(1-fd) # Beam irradiance
    Id0<-I*fd # Diffuse irradiance
    kd<-0.719 # Diffuse and scattered diffuse PAR extinction coefficient

    kb1[i]<-kb[i]*sqrt(1-lsc)
    kb2[i]<-0.46/sinB[i] # Beam and scattered beam PAR extinction coefficient
    pcb[i]<-1-exp(-2*ph*kb[i]/(1+kb[i])) # Canopy reflection coefficient for beam PAR


    # Absorbed irradiance of sunlit and shaded fractions of canopies
    X[i]<-Ib0*(1-lsc)*(1-exp(-kb[i]*LAIc))*sinB[i] # Direct-beam irradiance absorbed by sunlit leaves
    Y[i]<-Id0*(1-pcd)*(1-exp(-(kb[i]+kb1[i])*Lsun[i]))*(kd/(kd+kb[i]))*sinB[i] # Diffuse irradiance absorbed by sunlit leaves
    ZZ[i]<-Ib0*((1-pcb[i])*(1-exp(-(kb1[i]+kb[i])*Lsun[i]))*kb1[i]/(kb1[i]+kb[i])-(1-lsc)*(1-exp(-2*kb[i]*Lsun[i]))/2)*sinB[i]
    Z[i]<-ifelse(is.na(ZZ[i]),0,ZZ[i]) # Scattered-beam irradiance absorbed by sunlit leaves
    Icsun[i]<-ifelse(kb[i]==0,0,X[i]+Y[i]+Z[i]) # Irradiance absorbed by the sunlit fraction of the canopy
    Ic[i]<-ifelse(kb[i]==0,0,((1-pcb[i])*Ib0*(1-exp(-kb1[i]*LAIc))+(1-pcd)*Id0*(1-exp(-kd*LAIc)))*sinB[i]) # Total irradiance absorbed, per unit ground area
    Icsh[i]<-Ic[i]-Icsun[i] # Irradiance absorbed by the shaded fractions of canopy

    # Irradiance dependence of electron transport of sunlit and shaded fractions of canopies
    # bsh[i]<--(0.425*Icsh[i]+Jmax[i]) # Sum of PAR effectively absorbed by PSII and Jmax
    # csh[i]<-Jmax[i]*0.425*Icsh[i]
    #
    # a<-0.7 # Curvature of leaf response of electron transport to irradiance
    # Jsh1[i]<-(-bsh[i]+sqrt(bsh[i]^2-4*a*csh[i]))/(2*a) # Positive solution
    # Jsh2[i]<-(-bsh[i]-sqrt(bsh[i]^2-4*a*csh[i]))/(2*a) # Negative solution
    #
    # bsunsh[i]<--(0.425*Ic[i]+Jmax[i])
    # csunsh[i]<-0.425*Jmax[i]*Ic[i]
    # Jsunsh1[i]<-(-bsunsh[i]+sqrt(bsunsh[i]^2-4*a*csunsh[i]))/(2*a) # Positive solution
    # Jsunsh2[i]<-(-bsunsh[i]-sqrt(bsunsh[i]^2-4*a*csunsh[i]))/(2*a) # Negative solution

    Ratio1[i]<-ifelse(Lsun[i]==0,0,Icsh[i]/Ic[i])

    Ratio[i]<-ifelse(Ratio1[i]>0.4,0.4,Ratio1[i])

    # Gross photosynthesis rates for sunlit and shaded leaves
    Vjsh1[i]<-0.2*Vj[i]
    Vjsh2[i]<-Vj[i]*Ratio[i]
    Vjsh[i]<-ifelse(Vjsh1[i]>Vjsh2[i],Vjsh1[i],Vjsh2[i]) # Light-limited gross photosynthesis rates of shaded leaves, μmol m-2 s-1
    Vjsun[i]<-(Vj[i]-Vjsh[i]) # Light-limited gross photosynthesis rates of sunlit leaves, μmol m-2 s-1
    VcVjsun[i]<-ifelse(Vc[i]<Vjsun[i],Vc[i],Vjsun[i]) # Gross photosynthesis rates of sunlit leaves, μmol m-2 s-1
    VcVjsh[i]<-ifelse(Vc[i]<Vjsh[i],Vc[i],Vjsh[i]) # Gross photosynthesis rates of shaded leaves, μmol m-2 s-1


    # Net CO2 assimilation rate for sunlit and shaded leaves
    Rd[i]<-Vm[i]*0.015 # Leaf dark respiration, μmol m-2 s-1
    Asun[i]<-ifelse(VcVjsun[i]-Rd[i]<0,0,VcVjsun[i]-Rd[i]) # Net CO2 assimilation rate of sunlit leaves, μmol m-2 s-1
    Ash[i]<-ifelse(VcVjsh[i]-Rd[i]<0,0,VcVjsh[i]-Rd[i]) # Net CO2 assimilation rate of shaded leaves, μmol m-2 s-1
    #Asunsh[i]<-ifelse(Vc[i]<Vj[i],Vc[i],Vj[i]) # Just for test

    # Net CO2 assimilation rate for sunlit and shaded canopy
    Acsun[i]<-Asun[i]*Lsun[i] # Net CO2 assimilation rate for sunlit canopy, μmol m-2 s-1
    Acsh[i]<-Ash[i]*Lsh[i] # Net CO2 assimilation rate for shaded canopy, μmol m-2 s-1
    Acanopy[i]<-Acsun[i]+Acsh[i] # Net CO2 assimilation rate for canopy, μmol m-2 s-1

    # Stomatal model
    m<-Input_parameter$m # Coefficient
    g0<-Input_parameter$g0 # Initial stomatal conductance, m mol m-2 s-1
    S_swc<-Input_parameter$SWCs # Saturated soil volumetric moisture content at depth of 30 cm, %
    W_swc<-Input_parameter$SWCw # Wilting soil volumetric moisture content at depth of 30 cm, %
    VPD_close<-Input_parameter$VPD_close # The VPD at stomatal closure
    VPD_open<-Input_parameter$VPD_open # The VPD at stomatal opening

    gs[i] <- g0 + 100 * m * RH[i] * Acanopy[i] / Ca[i] # Stomatal conductance at leaf level, m mol m-2 s-1
    # gs[i]<-(g0+m*RH[i]*Acanopy[i]/Ca[i]*(max(0,min((SWC[i]-W_swc)/(S_swc-W_swc),1)))*(max(0,min((VPD_close-VPD1[i])/(VPD_close-VPD_open),1)))) # Stomatal conductance at leaf level, m mol m-2 s-1

    # Net CO2 assimilation rate for canopy by Leuning (1990)
    AV[i] <- (gs[i]*22.4/(8.314*(Tem[i]+273))*(Ca[i]-Ci[i])/1.6)*LAI/2

    #CJ4[i]<-ifelse(kb[i]==0,0,1)
    #fsun[i]<-exp(-kb[i]*LAIc)*CJ4[i]
    #fsh[i]<-1-fsun[i]


    # Gross primary production
    GPPsec[i]<-12/1000000*(Acanopy[i]+Rd[i])*1 # g C m-2 s-1


    # Ecosystem respiration
    # Respiration parameters
    Mf<-Input_parameter$Mf # Biomass density of for leaf, kg C m-2 day-1
    Ms<-Input_parameter$Ms # Biomass density of for sapwood, kg C m-2 day-1
    Mr<-Input_parameter$Mr # Biomass density of for root, kg C m-2 day-1

    rmf<-Input_parameter$rmf # Maintenance respiration coefficient for leaf
    rms<-Input_parameter$rms # Maintenance respiration coefficient for stem
    rmr<-Input_parameter$rmr # Maintenance respiration coefficient for root

    rgf<-Input_parameter$rgf # Growth respiration coefficient for leaf
    rgs<-Input_parameter$rgs # Growth respiration coefficient for sapwood
    rgr<-Input_parameter$rgr # Growth respiration coefficient for root

    raf<-Input_parameter$raf # Carbon allocation fraction for leaf
    ras<-Input_parameter$ras # Carbon allocation fraction for sapwood
    rar<-Input_parameter$rar # Carbon allocation fraction for root

    Q10<-Input_parameter$Q10 # Temperature sensitivity factor
    Tref<-Input_parameter$Tref # Base temperature for Q10, ℃

    # Maintenance respiration
    Rmf[i]<-Mf*rmf*Q10^((Tem[i]-Tref)/10) # Maintenance respiration for leaf, gC m-2 s-1
    Rms[i]<-Ms*rms*Q10^((Tem[i]-Tref)/10) # Maintenance respiration for stem, gC m-2 s-1
    Rmr[i]<-Mr*rmr*Q10^((Tem[i]-Tref)/10) # Maintenance respiration for root, gC m-2 s-1
    Rm[i]<-Rmf[i]+Rms[i]+Rmr[i] # Total maintenance respiration, gC m-2 s-1
    Rm30min[i]<-1000*Rm[i]/24/2 # Unit coversion, from gC m-2 30min-1

    # Growth respiration
    Rgf[i]<-rgf*raf*GPPsec[i] # Growth respiration for leaf, gC m-2 s-1
    Rgs[i]<-rgs*ras*GPPsec[i] # Growth respiration for sapwood, gC m-2 s-1
    Rgr[i]<-rgr*rar*GPPsec[i] # Growth respiration for root, gC m-2 s-1
    Rg[i]<-Rgf[i]+Rgs[i]+Rgr[i] # Total growth respiration, gC m-2 s-1
    Rg30min[i]<-1800*Rg[i] # Unit coversion, from gC m-2 30 min-1

    # Heterotrophic respiration
    R10<-1.5 # Soil respiration rate in 10 degrees Celsius
    Rh30min[i]<-(R10*Q10^((Tem[i]-10)/10))/48 # Heterotrophic respiration, gC m-2 30 min-1

    Re30min[i]<-Rh30min[i]+Rm30min[i]+Rg30min[i] # Ecosystem respiration, gC m-2 30 min-1

    # Carbon flux in ecosystem，碳通量计算
    GPP30min[i]<-ifelse(PPFD[i]==0,0,1800*GPPsec[i]) # Gross primary production, gC m-2 30min-1
    NPP30min[i]<--(GPP30min[i]-Rm30min[i]-Rg30min[i]) # Net primary production, gC m-2 30min-1
    NEP30min[i]<-GPP30min[i]-Re30min[i] # Net ecosystem production, gC m-2 30min-1

    # Evapotranspiration by Penman–Monteith model
    Cp<-1.013*1000 # Specific heat of the air, J kg-1 degrees Celsius-1
    Ve<-Input_variable$Vms # The wind speed at height Hw, m s-1
    Pdensity<-rep(1.29,x) # Air density, kg m-3
    r<-0.66/10 # Psychrometric constant, kPa degrees Celsius-1
    Hw<-rep(Input_parameter$Hw,x) # The height at wind measurement, m
    Hcanopy<-rep(Input_parameter$hc,x) # The average canopy height, m
    k<-rep(0.41,x) # Von karman’s constant

    gsms[i]<-gs[i]*22.4*10^-6*LAI/2 # Stomatal conductance at canopy level, m s-1
    Z0M[i]<-0.1*Hcanopy[i] # Roughness height for momentum transfer, m
    Z0V[i]<-0.5*Z0M[i] # Roughness height for vapor and heat transfer, m
    d[i]<-0.7*Hcanopy[i] # Zero plane displacement height, m
    svt[i]<-4098*(0.6108*exp(17.27*Tem[i]/(Tem[i]+237.3)))/(Tem[i]+237.3)^2 # The slope of the saturation vapor pressure against temperature curve, kPa degrees Celsius-1
    ra[i]<-log((Hw[i]-d[i])/Z0M[i])*log((Hw[i]-d[i])/Z0V[i])/(k[i]^2*Ve[i]) # Aerodynamic resistance, s m-1
    LES[i]<-(svt[i]*(Rn[i]-G[i])+Pdensity[i]*Cp*VPD1[i]/ra[i])/(svt[i]+r*(1+1/gsms[i]/ra[i])) # Latent heat, w m-2
    ETS[i]<-0.43*LES[i]/(597-0.564*Tem[i]) # Evapotranspiration, mm 30 min-1

    # print(paste(Acanopy[i],AV[i],sep = "--"))


    #Iteration condition of TRIPLEX_CW_Flux: (abs(Acanopy[i]-AV[i])>1)|(Acanopy[i]>3)
    while((abs(Acanopy[i]-AV[i])>1)|(Acanopy[i]>3)){

      # Iteration
      j <- j + 1
      Cippm[i] <- Cappm[i] - j # Setting initial value of intercellular CO2 concentration, assume:Cippm=Cappm
      Ca[i]<-Cappm[i]/1000000*1.013*100000 # Unit conversion: from ppm to Pa
      Ci[i]<-Cippm[i]/1000000*1.013*100000 # Unit conversion: from ppm to Pa

      Vm25<-Input_parameter$Vm25 # Maximum carboxylation rate at 25degrees Celsius, μmol m-2 s-1
      Rgas<-Input_parameter$Rgas # Molar gas constant, m3 Pa mol-1 K-1
      O2<-Input_parameter$O2 # Oxygen concentration in the atmosphere, Pa
      N<-Input_parameter$N # Leaf nitrogen content, %
      Nm<-Input_parameter$Nm # Maximum nitrogen content,%
      fN<-N/Nm # Nitrogen limitation term
      Error1[i]<-Ci[i]/Ca[i] # Just for check
      COp[i]<-1.92*10^-4*O2*1.75^((Tem[i]-25)/10) # CO2 concentration point without dark respiration, Pa
      Kc[i]<-30*2.1^((Tem[i]-25)/10) # Michaelis–Menten constants for CO2, Pa
      Ko[i]<-30000*1.2^((Tem[i]-25)/10) # Michaelis–Menten constants for O2, Pa
      K[i]<-Kc[i]*(1+O2/Ko[i]) # Function of enzyme kinetics, Pa
      fT[i]<-(1+exp((-220000+710*(Tem[i]+273))/(Rgas*(Tem[i]+273))))^-1 # Temperature limitation term
      Vm[i]<-Vm25*2.4^((Tem[i]-25)/10)*fT[i]*fN # Maximum carboxylation rate, μmol m-2 s-1
      Vc[i]<-Vm[i]*(Ci[i]-COp[i])/(Ci[i]+K[i]) # Rubisco-limited gross photosynthesis rates, μmol m-2 s-1

      J[i] <- SIFfull[i] * ((1 + Kdf) * qL[i] * PSmax[i]) / (1 - PSmax[i])

      #Jmax[i]<-29.1+1.64*Vm[i] # Light-saturated rate of electron transport in the photosynthetic carbon reduction cycle in leaf cells, μmol m-2 s-1
      #J[i]<-Jmax[i]*PPFD[i]/(PPFD[i]+2.1*Jmax[i]) # Electron transport rate, μmol m-2 s-1
      #J[i]<-Jmax[i]*PPFD[i]/(PPFD[i]+2.1*Jmax[i]) # Electron transport rate, μmol m-2 s-1

      Vj1[i]<-J[i]*(Ci[i]-COp[i])/(4.5*Ci[i]+10.5*COp[i])
      Vj[i]<-ifelse(Vj1[i]>0,Vj1[i],0) # Light-limited gross photosynthesis rates, μmol m-2 s-1
      # Vj[i]<-max(J[i]*(Ci[i]-COp[i])/(4.5*Ci[i]+10.5*COp[i]),0) # Light-limited gross photosynthesis rates, μmol m-2 s-1

      # Solar geometry
      Ls<-Input_parameter$Ls # Standard longitude of time zone
      Le<-Input_parameter$Le # Local longitude, degree
      latitude<-Input_parameter$latitude # Local latitude, degree
      lat<-latitude*3.14/180 # Convert angle to radian, radian
      LAI<-Input_parameter$LAI # Leaf area index, m2 m-2
      Ld[i]<-2*3.14*(Day[i]-1)/365 # Day angle
      Et[i]<-0.017+0.4281*cos(Ld[i])-7.351*sin(Ld[i])-3.349*cos(2*Ld[i])-9.731*sin(Ld[i]) # Equation of time, min
      to[i]<-12+(4*(Ls-Le)-Et[i])/60 # Solar noon, h
      h[i]<-3.14*(time[i]-to[i])/12 # Hour angle of sun, radians
      d[i]<--23.4*3.14/180*cos(2*3.14*(Day[i]+10)/365) # Solar declination angle, radians
      sinB[i]<-sin(lat)*sin(d[i])+cos(lat)*cos(d[i])*cos(h[i]) # Solar elevation angle, radians

      # Sun/Shade model
      # Leaf area index for sunlit and shaded of the canopy
      kb[i]<-ifelse(sinB[i]>0,0.5/sinB[i],0) # Beam radiation extinction coefficient of canopy
      LAIc<-LAI # Leaf area index of canopy, m2 m-2
      Lsun[i]<-ifelse(kb[i]==0,0,(1-exp(-kb[i]*LAIc))/kb[i]) # Sunlit leaf area index of canopy, m2 m-2
      Lsh[i]<-LAIc-Lsun[i] # Shaded leaf area index of canopy, m2 m-2

      # Canopy reflection coefficients
      lsc<-0.15 # Leaf scattering coefficient of PAR
      ph<-(1-sqrt(1-lsc))/(1+sqrt(1-lsc)) # Reflection coefficient of a canopy with horizontal leaves
      pcd<-0.036 # Canopy reflection coefficient for diffuse PAR
      I<-2083 # Total incident PAR μmol m-2 ground s-1
      fd<-0.159 # Fraction of diffuse irradiance
      Ib0<-I*(1-fd) # Beam irradiance
      Id0<-I*fd # Diffuse irradiance
      kd<-0.719 # Diffuse and scattered diffuse PAR extinction coefficient

      kb1[i]<-kb[i]*sqrt(1-lsc)
      kb2[i]<-0.46/sinB[i] # Beam and scattered beam PAR extinction coefficient
      pcb[i]<-1-exp(-2*ph*kb[i]/(1+kb[i])) # Canopy reflection coefficient for beam PAR


      # Absorbed irradiance of sunlit and shaded fractions of canopies
      X[i]<-Ib0*(1-lsc)*(1-exp(-kb[i]*LAIc))*sinB[i] # Direct-beam irradiance absorbed by sunlit leaves
      Y[i]<-Id0*(1-pcd)*(1-exp(-(kb[i]+kb1[i])*Lsun[i]))*(kd/(kd+kb[i]))*sinB[i] # Diffuse irradiance absorbed by sunlit leaves
      ZZ[i]<-Ib0*((1-pcb[i])*(1-exp(-(kb1[i]+kb[i])*Lsun[i]))*kb1[i]/(kb1[i]+kb[i])-(1-lsc)*(1-exp(-2*kb[i]*Lsun[i]))/2)*sinB[i]
      Z[i]<-ifelse(is.na(ZZ[i]),0,ZZ[i]) # Scattered-beam irradiance absorbed by sunlit leaves
      Icsun[i]<-ifelse(kb[i]==0,0,X[i]+Y[i]+Z[i]) # Irradiance absorbed by the sunlit fraction of the canopy
      Ic[i]<-ifelse(kb[i]==0,0,((1-pcb[i])*Ib0*(1-exp(-kb1[i]*LAIc))+(1-pcd)*Id0*(1-exp(-kd*LAIc)))*sinB[i]) # Total irratiance absorbed
      Icsh[i]<-Ic[i]-Icsun[i] # Irradiance absorbed by the shaded fractions of canopy

      # Irradiance dependence of electron transport of sunlit and shaded fractions of canopies
      # bsh[i]<--(0.425*Icsh[i]+Jmax[i]) # Sum of PAR effectively absorbed by PSII and Jmax
      # csh[i]<-Jmax[i]*0.425*Icsh[i]
      #
      # a<-0.7 # Curvature of leaf response of electron transport to irradiance
      # Jsh1[i]<-(-bsh[i]+sqrt(bsh[i]^2-4*a*csh[i]))/(2*a) # Positive solution
      # Jsh2[i]<-(-bsh[i]-sqrt(bsh[i]^2-4*a*csh[i]))/(2*a) # Negative solution
      #
      # bsunsh[i]<--(0.425*Ic[i]+Jmax[i])
      # csunsh[i]<-0.425*Jmax[i]*Ic[i]
      # Jsunsh1[i]<-(-bsunsh[i]+sqrt(bsunsh[i]^2-4*a*csunsh[i]))/(2*a) # Positive solution
      # Jsunsh2[i]<-(-bsunsh[i]-sqrt(bsunsh[i]^2-4*a*csunsh[i]))/(2*a) # Negative solution

      Ratio1[i]<-ifelse(Lsun[i]==0,0,Icsh[i]/Ic[i])

      Ratio[i]<-ifelse(Ratio1[i]>0.4,0.4,Ratio1[i])

      # Gross photosynthesis rates for sunlit and shaded leaves
      Vjsh1[i]<-0.2*Vj[i]#
      Vjsh2[i]<-Vj[i]*Ratio[i]#
      Vjsh[i]<-ifelse(Vjsh1[i]>Vjsh2[i],Vjsh1[i],Vjsh2[i]) # Light-limited gross photosynthesis rates of shaded leaves, umol m-2 s-1
      Vjsun[i]<-(Vj[i]-Vjsh[i]) # Light-limited gross photosynthesis rates of sunlit leaves, μmol m-2 s-1
      VcVjsun[i]<-ifelse(Vc[i]<Vjsun[i],Vc[i],Vjsun[i]) # Gross photosynthesis rates of sunlit leaves, μmol m-2 s-1
      VcVjsh[i]<-ifelse(Vc[i]<Vjsh[i],Vc[i],Vjsh[i]) # Gross photosynthesis rates of shaded leaves, μmol m-2 s-1


      # Net CO2 assimilation rate for sunlit and shaed leaves
      Rd[i]<-Vm[i]*0.015 # Leaf dark respiration, μmol m-2 s-1
      Asun[i]<-ifelse(VcVjsun[i]-Rd[i]<0,0,VcVjsun[i]-Rd[i]) # Net CO2 assimilation rate of sunlit leaves, μmol m-2 s-1
      Ash[i]<-ifelse(VcVjsh[i]-Rd[i]<0,0,VcVjsh[i]-Rd[i]) # Net CO2 assimilation rate of shaded leaves, μmol m-2 s-1
      #Asunsh[i]<-ifelse(Vc[i]<Vj[i],Vc[i],Vj[i]) # Just for test

      # Net CO2 assimilation rate for sunlit and shaded canopy
      Acsun[i]<-Asun[i]*Lsun[i] # Net CO2 assimilation rate for sunlit canopy, μmol m-2 s-1
      Acsh[i]<-Ash[i]*Lsh[i] # Net CO2 assimilation rate for shaded canopy, μmol m-2 s-1
      Acanopy[i]<-Acsun[i]+Acsh[i] # Net CO2 assimilation rate for canopy, μmol m-2 s-1

      # Stomatal model
      m<-Input_parameter$m # Coefficient
      g0<-Input_parameter$g0 # Initial stomatal conductance, m mol m-2 s-1
      S_swc<-Input_parameter$SWCs # Saturated soil volumetric moisture content at depth of 30 cm, %
      W_swc<-Input_parameter$SWCw # Wilting soil volumetric moisture content at depth of 30 cm, %
      VPD_close<-Input_parameter$VPD_close # The VPD at stomatal closure
      VPD_open<-Input_parameter$VPD_open # The VPD at stomatal opening
      gs[i]<-(g0+m*RH[i]*Acanopy[i]/Ca[i]*(max(0,min((SWC[i]-W_swc)/(S_swc-W_swc),1)))*(max(0,min((VPD_close-VPD1[i])/(VPD_close-VPD_open),1)))) # Stomatal conductance at leaf level, m mol m-2 s-1

      # Net CO2 assimilation rate for canopy by Leuning (1990)
      AV[i]<-(gs[i]*22.4/(8.314*(Tem[i]+273))*(Ca[i]-Ci[i])/1.6)*LAI/2

      #CJ4[i]<-ifelse(kb[i]==0,0,1)#
      #fsun[i]<-exp(-kb[i]*LAIc)*CJ4[i]#
      #fsh[i]<-1-fsun[i]#

      # Gross primary production
      GPPsec[i]<-12/1000000*(Acanopy[i]+Rd[i])*1# g C m-2 s-1

      # Ecosystem respiration
      # Respiration parameters
      Mf<-Input_parameter$Mf # Biomass density of for leaf, kg C m-2 day-1
      Ms<-Input_parameter$Ms # Biomass density of for sapwood, kg C m-2 day-1
      Mr<-Input_parameter$Mr # Biomass density of for root, kg C m-2 day-1

      rmf<-Input_parameter$rmf # Maintenance respiration coefficient for leaf
      rms<-Input_parameter$rms # Maintenance respiration coefficient for stem
      rmr<-Input_parameter$rmr # Maintenance respiration coefficient for root

      rgf<-Input_parameter$rgf # Growth respiration coefficient for leaf
      rgs<-Input_parameter$rgs # Growth respiration coefficient for sapwood
      rgr<-Input_parameter$rgr # Growth respiration coefficient for root

      raf<-Input_parameter$raf # Carbon allocation fraction for leaf
      ras<-Input_parameter$ras # Carbon allocation fraction for sapwood
      rar<-Input_parameter$rar # Carbon allocation fraction for root

      Q10<-Input_parameter$Q10 # Temperature sensitivity factor
      Tref<-Input_parameter$Tref # Base temperature for Q10

      # Maintenance respiration
      Rmf[i]<-Mf*rmf*Q10^((Tem[i]-Tref)/10) # Maintenance respiration for leaf, gC m-2 s-1
      Rms[i]<-Ms*rms*Q10^((Tem[i]-Tref)/10) # Maintenance respiration for stem, gC m-2 s-1
      Rmr[i]<-Mr*rmr*Q10^((Tem[i]-Tref)/10) # Maintenance respiration for root, gC m-2 s-1
      Rm[i]<-Rmf[i]+Rms[i]+Rmr[i] # Total maintenance respiration, gC m-2 s-1
      Rm30min[i]<-1000*Rm[i]/24/2 # Unit coversion, from gC m-2 30 min-1

      # Growth respiration
      Rgf[i]<-rgf*raf*GPPsec[i] # Growth respiration for leaf, gC m-2 s-1
      Rgs[i]<-rgs*ras*GPPsec[i] # Growth respiration for sapwood, gC m-2 s-1
      Rgr[i]<-rgr*rar*GPPsec[i] # Growth respiration for root, gC m-2 s-1
      Rg[i]<-Rgf[i]+Rgs[i]+Rgr[i] # Total growth respiration, gC m-2 s-1
      Rg30min[i]<-1800*Rg[i] # Unit coversion, from gC m-2 30 min-1

      # Heterotrophic respiration
      R10<-1.5# Soil respiration rate in 10 degrees Celsius
      Rh30min[i]<-(R10*Q10^((Tem[i]-10)/10))/48 # Heterotrophic respiration, gC m-2 30 min-1

      Re30min[i]<-Rh30min[i]+Rm30min[i]+Rg30min[i] # Ecosystem respiration, gC m-2 30 min-1

      # Carbon flux in ecosystem
      GPP30min[i]<-ifelse(PPFD[i]==0,0,1800*GPPsec[i]) # Gross primary production, gC m-2 30 min-1
      NPP30min[i]<--(GPP30min[i]-Rm30min[i]-Rg30min[i]) # Net primary production, gC m-2 30 min-1
      NEP30min[i]<-GPP30min[i]-Re30min[i] # Net ecosystem production, gC m-2 30 min-1

      # Evapotranspiration by Penman–Monteith model
      Cp<-1.013*1000 # Specific heat of the air, J kg-1 degrees Celsius-1
      Ve<-Input_variable$Vms # The wind speed at height Hw, m s-1
      Pdensity<-rep(1.29,x) # Air density, kg m-3
      r<-0.66/10 # Psychrometric constant, kPa ℃-1
      Hw<-rep(Input_parameter$Hw,x) # The height at wind measurement, m
      Hcanopy<-rep(Input_parameter$hc,x) # The average canopy height, m
      k<-rep(0.41,x) # Von karman’s constant

      gsms[i]<-gs[i]*22.4*10^-6*LAI/2 # Stomatal conductance at canopy level, m s-1
      Z0M[i]<-0.1*Hcanopy[i] # Roughness height for momentum transfer, m
      Z0V[i]<-0.5*Z0M[i] # Roughness height for vapor and heat transfer, m
      d[i]<-0.7*Hcanopy[i] # Zero plane displacement height, m
      svt[i]<-4098*(0.6108*exp(17.27*Tem[i]/(Tem[i]+237.3)))/(Tem[i]+237.3)^2 # The slope of the saturation vapor pressure against temperature curve, kPa degrees Celsius-1
      ra[i]<-log((Hw[i]-d[i])/Z0M[i])*log((Hw[i]-d[i])/Z0V[i])/(k[i]^2*Ve[i]) # Aerodynamic resistance, s m-1
      LES[i]<-(svt[i]*(Rn[i]-G[i])+Pdensity[i]*Cp*VPD1[i]/ra[i])/(svt[i]+r*(1+1/gsms[i]/ra[i])) # Latent heat, w m-2
      ETS[i]<-0.43*LES[i]/(597-0.564*Tem[i]) # Evapotranspiration, mm 30 min-1

      # Output condition of TRIPLEX_CW_Flux model: ((abs(Acanopy[i]-AV[i])<1)|(Acanopy[i]<3))
      if ((abs(Acanopy[i]-AV[i])<1)|(Acanopy[i]<3)){
        break
      }
      # print(j)

    }
    # print(J[i])
    print(i)
  }
  result<-data.frame(Input_variable, NEP30min, ETS, GPP30min, Re30min, SIFfull, PSmax, qL, J, Ci, COp)
  #View(result)----

  # Overall simulated plot----
  library(ggplot2)
  library(ggpubr)
  library(dplyr)

  regre <- lm(GPP30min ~ GPP, result)
  r <- formatC(summary(regre)$r.squared, format = "f", digits = 2)
  rmse <- formatC(sqrt(sum(residuals(regre)^2) / (nrow(result) - 2)), format = "f", digits = 2)
  ia <- formatC(1 - sum((result$GPP30min - result$GPP)^2) /
                  sum((abs(result$GPP30min - mean(result$GPP)) +
                         abs(result$GPP - mean(result$GPP)))^2), format = "f", digits = 2)
  # 获取线性回归方程的系数，并保留两位小数
  intercept <- formatC(coef(regre)[1], format = "f", digits = 2)
  slope <- formatC(coef(regre)[2], format = "f", digits = 2)
  ## 创建ggplot基础图
  p_overall <- ggplot(result, aes(x= GPP, y = GPP30min)) +
    geom_point(shape = 21, size = 3, color = 'blue', stroke = 1) + # 散点图，边缘蓝色，填充白色
    geom_smooth(method = 'lm', color = 'red', se = FALSE, size = 1) + # 线性拟合线
    geom_abline(intercept = 0, slope = 1, color = 'grey60', size = 1, linetype = 'dashed') + # 1:1对线线
    labs(x = expression(Observed~GPP~(g~C~m^-2~30~min^-1)),
         y = expression(Simulated~GPP~(g~C~m^-2~30~min^-1))) +  # 坐标轴标签

    theme_minimal(base_size = 15) +  # 设置主题
    theme(axis.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 15),
          panel.grid = element_line(),
          panel.border = element_rect(fill = NA, color = 'black',size = 0.5))
  ## 设置坐标轴范围和刻度
  min_val <- min(result$GPP, result$GPP30min) - 0.1
  max_val <- max(result$GPP, result$GPP30min) * 1.1
  p_overall <- p_overall +
    scale_x_continuous(limits = c(0, max_val), breaks = seq(0, max_val, by = 0.2)) +  # 限制坐标轴范围并设置刻度
    scale_y_continuous(limits = c(0, max_val), breaks = seq(0, max_val, by = 0.2))
  ## 添加r方，RMSE，和IA的注释
  p_overall <- p_overall +
    annotate("text", x = 0.05,
             y = max_val * 0.95,
             label = paste("y =", slope, "x +", intercept),size = 5, hjust = 0) +
    annotate("text", x = 0.05,
             y = max_val * 0.85,
             label = bquote(R^2 == .(r)), size = 5, hjust = 0) +
    annotate("text", x = 0.05,
             y = max_val * 0.75,
             label = paste("RMSE =", rmse), size = 5, hjust = 0) +
    annotate("text", x = 0.05,
             y = max_val * 0.65,
             label = paste("IA =", ia), size = 5, hjust = 0) +
    annotate("text", x = 0,
             y = max_val * 0.95,
             label = "(b)", size = 5, hjust = 0)
  print(p_overall)
  ggsave("./result/plot_overall.jpg", plot = p_overall, width = 8, height = 6)


  # Seasonal simulated plot----

  # 添加季节列
  output_SIF_season <- result %>%
    mutate(Season = case_when(
      Month %in% 3:5 ~ "Spring",
      Month %in% 6:8 ~ "Summer",
      Month %in% 9:11 ~ "Autumn",
      Month == 12 ~ "Winter",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Season)) %>%
    mutate(Season = factor(Season, levels = c("Spring", "Summer", "Autumn", "Winter")))

  # 创建ggplot基础图
  p_season <- ggplot(output_SIF_season, aes(x = GPP, y = GPP30min)) +
    geom_point(shape = 21, size = 3, color = 'blue', stroke = 1) + # 散点图，边缘蓝色，填充白色
    geom_smooth(method = 'lm', color = 'red', se = FALSE, size = 1) + # 线性拟合线
    geom_abline(intercept = 0, slope = 1, color = 'grey60', size = 1, linetype = 'dashed') + # 1:1对线线
    labs(x = expression(Observed~GPP~(g~C~m^-2~30~min^-1)),
         y = expression(Simulated~GPP~(g~C~m^-2~30~min^-1))) +  # 坐标轴标签
    theme_minimal(base_size = 15) +  # 设置主题
    theme(axis.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 15),
          panel.grid = element_line(),
          panel.border = element_rect(fill = NA, color = 'black', size = 0.5)) +
    facet_wrap(~ Season, nrow = 2, ncol = 2)  # 按季节分面

  # 设置坐标轴范围和刻度
  min_val <- min(output_SIF_season$GPP, output_SIF_season$GPP30min) - 0.1
  max_val <- max(output_SIF_season$GPP, output_SIF_season$GPP30min) * 1.1
  p_season <- p_season +
    scale_x_continuous(limits = c(0, max_val), breaks = seq(0, max_val, by = 0.2)) +  # 限制坐标轴范围并设置刻度
    scale_y_continuous(limits = c(0, max_val), breaks = seq(0, max_val, by = 0.2))

  # 计算斜率和截距，并保留两位小数
  season_stats <- output_SIF_season %>%
    group_by(Season) %>%
    summarise(
      slope = round(coef(lm(GPP30min ~ GPP))[2], 2),
      intercept = round(coef(lm(GPP30min ~ GPP))[1], 2),
      r_squared = round(summary(lm(GPP30min ~ GPP))$r.squared, 2),
      rmse = round(sqrt(mean(residuals(lm(GPP30min ~ GPP))^2)), 2),
      ia = round(1 - sum((GPP30min - GPP)^2) / sum((abs(GPP30min - mean(GPP)) + abs(GPP - mean(GPP)))^2), 2)
    )

  # 添加线性公式、r方，RMSE，和IA的注释
  season_stats <- season_stats %>%
    mutate(label = paste0("y = ", slope, "x + ", intercept,
                          "\nR² = ", r_squared,
                          "\nRMSE = ", rmse,
                          "\nIA = ", ia))
  # 在ggplot中添加文本标签
  p_season <- p_season +
    geom_text(data = season_stats, aes(x = 0.05, y = max_val * 0.95, label = label),
              hjust = 0, vjust = 1, size = 5, color = "black")

  print(p_season)
  ggsave("./result/plot_Seasonal.jpg", plot = p_season, width = 12, height = 10)

  # Diurnal dynamics plot----

  # 去除 time 为 5.5, 6, 6.5, 19, 19.5 的数据
  filtered_data <- result %>%
    filter(!time %in% c(5.5, 19.5))

  # 按 time 分类汇总，计算 GPP30min、GPP 和 SIFfull 的平均值及误差条范围
  daily_stats <- filtered_data %>%
    group_by(time) %>%
    summarise(
      GPP30min_avg = mean(GPP30min, na.rm = TRUE),
      GPP30min_sd = sd(GPP30min, na.rm = TRUE),  # 计算 GPP30min 的标准差
      GPP_avg = mean(GPP, na.rm = TRUE),
      GPP_sd = sd(GPP, na.rm = TRUE),  # 计算 GPP 的标准差
      SIFfull_avg = mean(SIFfull, na.rm = TRUE),
      SIFfull_sd = sd(SIFfull, na.rm = TRUE)  # 计算 SIFfull 的标准差
    ) %>%
    mutate(
      GPP30min_min = GPP30min_avg - GPP30min_sd,  # 平均值减去标准差
      GPP30min_max = GPP30min_avg + GPP30min_sd,  # 平均值加上标准差
      GPP_min = GPP_avg - GPP_sd,  # 平均值减去标准差
      GPP_max = GPP_avg + GPP_sd,  # 平均值加上标准差
      SIFfull_min = SIFfull_avg - SIFfull_sd,  # 平均值减去标准差
      SIFfull_max = SIFfull_avg + SIFfull_sd,   # 平均值加上标准差
      GPP_min = pmax(GPP_min, 0),  # 将 GPP的最小值限制为 0
      GPP30min_min = pmax(GPP30min_min, 0),  # 将 GPP30min 的最小值限制为 0
      SIFfull_min = pmax(SIFfull_min, 0)    # 将 SIFfull 的最小值限制为 0
    )

  # 绘制 GPP30min 的日内变化图（带误差条）
  plot_gpp30min <- ggplot(daily_stats, aes(x = time, y = GPP30min_avg)) +
    geom_line(aes(color = "Simulated GPP"), size = 1) +  # 平均值折线图
    geom_point(aes(color = "Simulated GPP", shape = "Simulated GPP"), size = 2) +  # 平均值点
    geom_errorbar(aes(ymin = GPP30min_min, ymax = GPP30min_max, color = "Simulated GPP"), width = 0.2) +  # 误差条
    scale_x_continuous(
      breaks = seq(6, max(daily_stats$time, na.rm = TRUE), by = 2),  # x刻度从6开始，间隔2
      name = "Hour of day, 2023"
    ) +
    scale_color_manual(values = c("Simulated GPP" = "red")) +
    scale_shape_manual(values = c("Simulated GPP" = 16)) +
    labs(
      x = "Hour of day, 2023",
      y = expression(Siumlated~GPP~(g~C~m^-2~day^-1)),
      color = "Legend",
      shape = "Legend"
    ) +
    theme_minimal() +
    theme(
      legend.position = c(0, 1),  # 图例放置在左上角
      legend.justification = c(0, 1),
      legend.title = element_blank(),
      legend.text = element_text(size = 15),
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15),
      panel.border = element_rect(color = "black", fill = NA)  # 添加边框
    ) +
    annotate("text", x = max(daily_stats$time, na.rm = TRUE), y = max(daily_stats$GPP30min_max, na.rm = TRUE) * 1.1,
             label = "(a)", color = "red", size = 5, hjust = 1)

  # 绘制 GPP 的日内变化图（带误差条）
  plot_gpp <- ggplot(daily_stats, aes(x = time, y = GPP_avg)) +
    geom_line(aes(color = "Observed GPP"), size = 1) +  # 平均值折线图
    geom_point(aes(color = "Observed GPP", shape = "Observed GPP"), size = 2) +  # 平均值点
    geom_errorbar(aes(ymin = GPP_min, ymax = GPP_max, color = "Observed GPP"), width = 0.2) +  # 误差条
    scale_x_continuous(
      breaks = seq(6, max(daily_stats$time, na.rm = TRUE), by = 2),  # x刻度从6开始，间隔2
      name = "Hour of day, 2023"
    ) +
    scale_color_manual(values = c("Observed GPP" = "blue")) +
    scale_shape_manual(values = c("Observed GPP" = 16)) +
    labs(
      x = "Hour of day, 2023",
      y = expression(Observed~GPP~(g~C~m^-2~30~min^-1)),
      color = "Legend",
      shape = "Legend"
    ) +
    theme_minimal() +
    theme(
      legend.position = c(0, 1),  # 图例放置在左上角
      legend.justification = c(0, 1),
      legend.title = element_blank(),
      legend.text = element_text(size = 15),
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15),
      panel.border = element_rect(color = "black", fill = NA)  # 添加边框
    ) +
    annotate("text", x = max(daily_stats$time, na.rm = TRUE), y = max(daily_stats$GPP_max, na.rm = TRUE) * 1.1,
             label = "(b)", color = "blue", size = 5, hjust = 1)

  # 绘制 GPP30min 和 GPP 的叠加图
  plot_gpp_combined <- ggplot(daily_stats, aes(x = time)) +
    geom_line(aes(y = GPP_avg, color = "Observed GPP"), size = 1) +
    geom_point(aes(y = GPP_avg, color = "Observed GPP", shape = "Observed GPP"), size = 2) +
    geom_errorbar(aes(ymin = GPP_min, ymax = GPP_max, color = "Observed GPP"), width = 0.2) +
    geom_line(aes(y = GPP30min_avg, color = "Simulated GPP"), size = 1) +
    geom_point(aes(y = GPP30min_avg, color = "Simulated GPP", shape = "Simulated GPP"), size = 2) +
    geom_errorbar(aes(ymin = GPP30min_min, ymax = GPP30min_max, color = "Simulated GPP"), width = 0.2) +
    scale_x_continuous(
      breaks = seq(6, max(daily_stats$time, na.rm = TRUE), by = 2),
      name = "Hour of day, 2023"
    ) +
    scale_color_manual(values = c( "Observed GPP" = "blue","Simulated GPP" = "red")) +
    scale_shape_manual(values = c( "Observed GPP" = 16,"Simulated GPP" = 16)) +
    labs(
      x = "Hour of day, 2023",
      y = expression(GPP~(g~C~m^-2~30~min^-1)),
      color = "Legend",
      shape = "Legend"
    ) +
    theme_minimal() +
    theme(
      legend.position = c(0, 1),
      legend.justification = c(0, 1),
      legend.title = element_blank(),
      legend.text = element_text(size = 15),
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15),
      panel.border = element_rect(color = "black", fill = NA)
    ) +
    annotate("text", x = max(daily_stats$time, na.rm = TRUE), y = max(daily_stats$GPP30min_max, na.rm = TRUE) * 1.1,
             label = "(a)", color = "black", size = 8, hjust = 1)

  # 计算 GPP30min 和 SIFfull 的最大值，并计算转换因子 coef
  max_gpp30min <- max(daily_stats$GPP30min_max, na.rm = TRUE)
  max_siffull <- max(daily_stats$SIFfull_max, na.rm = TRUE)
  coef <- max_gpp30min / max_siffull  # 转换因子，确保左右轴长度一致

  # 绘制 GPP30min 和 SIFfull 的叠加图（双轴）
  plot_gpp_siffull_combined <- ggplot(daily_stats, aes(x = time)) +
    geom_line(aes(y = SIFfull_avg * coef, color = "SIFfull"), size = 1) +
    geom_point(aes(y = SIFfull_avg * coef, color = "SIFfull", shape = "SIFfull"), size = 2) +
    geom_errorbar(aes(ymin = SIFfull_min * coef, ymax = SIFfull_max * coef, color = "SIFfull"), width = 0.2) +
    geom_line(aes(y = GPP30min_avg, color = "Simulated GPP"), size = 1) +
    geom_point(aes(y = GPP30min_avg, color = "Simulated GPP", shape = "Simulated GPP"), size = 2) +
    geom_errorbar(aes(ymin = GPP30min_min, ymax = GPP30min_max, color = "Simulated GPP"), width = 0.2) +
    scale_x_continuous(
      breaks = seq(6, max(daily_stats$time, na.rm = TRUE), by = 2),
      name = "Hour of day, 2023"
    ) +
    scale_y_continuous(
      name = NULL,
      limits = c(0, max_gpp30min * 1.1),
      sec.axis = sec_axis(~ . / coef, name = expression(SIFfull~(u~mol~m^-2~s^-1)))  # 右轴反向缩放
    ) +
    scale_color_manual(values = c( "SIFfull" = "skyblue", "Simulated GPP" = "red")) +
    scale_shape_manual(values = c("SIFfull" = 16, "Simulated GPP" = 16)) +
    labs(
      x = "Hour of day, 2023",
      color = "Legend",
      shape = "Legend"
    ) +
    theme_minimal() +
    theme(
      legend.position = c(0, 1),
      legend.justification = c(0, 1),
      legend.title = element_blank(),
      legend.text = element_text(size = 15),
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15),
      axis.title.y.right = element_text(angle = 90),  # 右轴标题从下到上排列
      panel.border = element_rect(color = "black", fill = NA)
    ) +
    annotate("text", x = max(daily_stats$time, na.rm = TRUE), y = max_gpp30min * 1.1,
             label = "(b)", color = "black", size = 8, hjust = 1)

  # 保存叠加图表
  ggsave("./result/plot_diurnal_GPP_Combined.jpg", plot = plot_gpp_combined, width = 8, height = 6)
  ggsave("./result/plot_diurnal_GPP_SIFfull_Combined.jpg", plot = plot_gpp_siffull_combined, width = 8, height = 6)
  print(plot_gpp_combined)
  print(plot_gpp_siffull_combined)

  # Environmental factors plot----
  # 加载必要的包
  library(egg)
  library(grid)

  # 按Month分类，对Ta、VPD、RH、Rn、SIF进行求和汇总
  monthly_summary <- result %>%
    group_by(Month) %>%
    summarise(
      Ta = mean(Ta, na.rm = TRUE),
      VPD = mean(VPDhpa, na.rm = TRUE),
      RH = mean(RH, na.rm = TRUE),
      Rn = mean(Rn, na.rm = TRUE),
      SIF = mean(SIF, na.rm = TRUE)
    )

  # 确保 Month 按月份大小升序排序
  monthly_summary$Month <- factor(monthly_summary$Month,levels = unique(monthly_summary$Month))

  # 绘制 Ta 的折线图（实心圆点）
  plot_ta <- ggplot(monthly_summary, aes(x = Month, y = Ta, group = 1)) +
    geom_point(aes(shape = "Ta"), size = 3, color = "black", na.rm = TRUE) +  # 改为黑色
    geom_line(aes(linetype = "Ta"), color = "black", na.rm = TRUE) +  # 改为黑色
    scale_shape_manual(values = c("Ta" = 16)) +
    scale_linetype_manual(values = c("Ta" = "solid")) +
    scale_y_continuous(name = "Ta (℃)", limits = c(0, max(monthly_summary$Ta) * 1.1)) +
    # scale_x_discrete(labels = function(x) gsub("月", "", x)) +  # 去掉月份中的“月”
    labs(x = "Month", shape = "Legend", linetype = "Legend") +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA),  # 添加边框
          legend.position = c(0.9, 0.95),  # 图例放在左上角
          legend.text = element_text(size = 15),  # 图例文本大小
          legend.title = element_blank(),  # 隐藏图例标题
          axis.text = element_text(size = 15),  # 调整坐标轴文本字号
          axis.title = element_text(size = 15),  # 调整轴标题字号
          axis.title.y.right = element_text(angle = 90))  # 右 y 轴标题从下到上排列

  # 计算 RH 和 VPD 的最大值，并计算转换因子 coef_rh_vpd
  max_rh <- max(monthly_summary$RH, na.rm = TRUE)
  max_vpd <- max(monthly_summary$VPD, na.rm = TRUE)
  coef_rh_vpd <- max_rh / max_vpd  # 计算 RH 和 VPD 的转换因子

  # 绘制 RH 和 VPD 的折线图（空心圆点和实心三角），实现双 y 轴
  plot_rh_vpd <- ggplot(monthly_summary, aes(x = Month)) +
    geom_point(aes(y = RH, shape = "RH"), size = 3, color = "black", na.rm = TRUE) +  # 改为黑色
    geom_line(aes(y = RH, linetype = "RH", group = 1), color = "black", na.rm = TRUE) +  # 改为黑色
    geom_point(aes(y = VPD * coef_rh_vpd, shape = "VPD"), size = 3, color = "black", na.rm = TRUE) +  # 改为黑色
    geom_line(aes(y = VPD * coef_rh_vpd, linetype = "VPD", group = 1), color = "black", na.rm = TRUE) +  # 改为黑色
    scale_shape_manual(values = c("RH" = 1, "VPD" = 17)) +
    scale_linetype_manual(values = c("RH" = "solid", "VPD" = "dashed")) +
    scale_y_continuous(
      name = "RH (%)",
      limits = c(0, max_rh * 1.1),
      sec.axis = sec_axis(~ . / coef_rh_vpd, name = "VPD (hPa)")  # 反向转换 VPD 的值
    ) +
    # scale_x_discrete(labels = function(x) gsub("月", "", x)) +  # 去掉月份中的“月”
    labs(x = "Month", shape = "Legend", linetype = "Legend") +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA),  # 添加边框
          legend.position = c(0.9, 0.95),  # 图例放在左上角
          legend.title = element_blank(),  # 隐藏图例标题
          legend.text = element_text(size = 15),  # 图例文本大小
          axis.text = element_text(size = 15),  # 调整坐标轴文本字号
          axis.title = element_text(size = 15),  # 调整轴标题字号
          axis.title.y.right = element_text(angle = 90))  # 右 y 轴标题从下到上排列

  # 计算 Rn 和 SIF 的最大值，并计算转换因子 coef
  max_rn <- max(monthly_summary$Rn, na.rm = TRUE)
  max_sif <- max(monthly_summary$SIF, na.rm = TRUE)
  coef <- max_rn / max_sif

  # 绘制 Rn 和 SIF 的折线图（空心矩形和*号），SIF 乘以转换因子 coef
  plot_rn_sif <- ggplot(monthly_summary, aes(x = Month)) +
    geom_point(aes(y = Rn, shape = "Rn"), size = 3, color = "black", na.rm = TRUE) +  # 改为黑色
    geom_line(aes(y = Rn, linetype = "Rn", group = 1), color = "black", na.rm = TRUE) +  # 改为黑色
    geom_point(aes(y = SIF * coef, shape = "SIF"), size = 3, color = "black", na.rm = TRUE) +  # 改为黑色
    geom_line(aes(y = SIF * coef, linetype = "SIF", group = 1), color = "black", na.rm = TRUE) +  # 改为黑色
    scale_shape_manual(values = c("Rn" = 0, "SIF" = 8)) +
    scale_linetype_manual(values = c("Rn" = "solid", "SIF" = "dotted")) +
    scale_y_continuous(
      name = expression(Rn ~ (W~m^-2)),
      limits = c(0, max_rn * 1.1),
      sec.axis = sec_axis(~ . / coef, name = expression(SIF ~ (mW~m^-2~nm^-1~sr^-1)))
    ) +
    # scale_x_discrete(labels = function(x) gsub("月", "", x)) +  # 去掉月份中的“月”
    labs(x = "Month", shape = "Legend", linetype = "Legend") +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA),  # 添加边框
          legend.position = c(0.9, 0.95),  # 图例放在左上角
          legend.title = element_blank(),  # 隐藏图例标题
          legend.text = element_text(size = 15),  # 图例文本大小
          axis.text = element_text(size = 15),  # 调整坐标轴文本字号
          axis.title = element_text(size = 15),  # 调整轴标题字号
          axis.title.y.right = element_text(angle = 90))  # 右 y 轴标题从下到上排列


  # 设置绘图区大小的函数
  save_plot_with_fixed_panel_size <- function(plot, filename, panel_width, panel_height, width, height, dpi = 300) {
    # 使用 egg::set_panel_size 调整绘图区大小
    adjusted_plot <- egg::set_panel_size(
      plot,
      width = unit(panel_width, "in"),
      height = unit(panel_height, "in")
    )

    # 创建一个新的绘图设备
    png(filename = filename, width = width, height = height, units = "in", res = dpi)

    # 绘制调整后的图
    grid.draw(adjusted_plot)

    # 关闭绘图设备
    dev.off()
  }

  # 设置绘图区大小和图片大小
  panel_width <- 6.5  # 绘图区宽度（单位：英寸）
  panel_height <- 5  # 绘图区高度（单位：英寸）
  image_width <- 8  # 图片总宽度（单位：英寸）
  image_height <- 6  # 图片总高度（单位：英寸）

  print(plot_ta)
  # 保存 Ta 图表
  save_plot_with_fixed_panel_size(
    plot = plot_ta,
    filename = "./result/plot_environment_Ta.png",
    panel_width = panel_width,
    panel_height = panel_height,
    width = image_width,
    height = image_height
  )

  print(plot_rh_vpd)
  # 保存 RH 和 VPD 图表
  save_plot_with_fixed_panel_size(
    plot = plot_rh_vpd,
    filename = "./result/plot_environment_RH_VPD.png",
    panel_width = panel_width,
    panel_height = panel_height,
    width = image_width,
    height = image_height
  )

  print(plot_rn_sif)
  # 保存 Rn 和 SIF 图表
  save_plot_with_fixed_panel_size(
    plot = plot_rn_sif,
    filename = "./result/plot_environment_Rn_SIF.png",
    panel_width = panel_width,
    panel_height = panel_height,
    width = image_width,
    height = image_height
  )


}

