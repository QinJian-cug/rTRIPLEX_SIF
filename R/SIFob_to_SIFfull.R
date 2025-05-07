library(pracma)
RC_PSII <- read.csv('../data/mat2csv/RC_PSII.csv') # Read the standard database of chlorophyll fluorescence by relative path 

# SIFob convert to SIFfull according to Liu(2022)
fPS2 <- Input_parameter$fPS2 # The contribution of PSII fluorescence to SIF
fPAR[i] <- 0.79 * (0.516 * WDRVI[i] + 0.726) # The radio of absorbed photosynthesis active radiation of canopy
fesc[i] <- NDVI[i] * NIR[i] / fPAR[i] # The escaped probability of SIF photon from leaf level to canopy level
svd_result <- svd(RC_PSII) # Using SVD method to decompose the standard database of chlorophyll fluorescence
fc <- svd_result$v # Function can convert SIF from PSII wavelength to total wavelength by using SVD 

SIFps2[i] <- SIFob[i] * fPS2 # Convert observed SIF to SIF in PSII
SIFps[i] <- SIFps2[i] / (0.9 * fesc[i]) # Downscale SIF from canopy scale to photosystem scale
SIFfc_vector[i] <- SIFps[i] / fc[760-640+1,1] * fc # Upscale the single point SIF to full wavelength SIF by using vector calculate
SIFfc[i] <- sum(SIFfc_vector[i]) # Accumulate all the full wavelength SIF in the vector to a whole SIF  
SIFfull[i] <- SIFfc[i] / (6.63 * 3 * 6002 * 10^3) # Convert the units by using simplified constants including Plank constant, light speed, Avogadro constant