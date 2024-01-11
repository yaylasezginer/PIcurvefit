Matlab script to analyze data from Photosynthesis-Irradiance experiments. Fit a curve to data to describe photosynthesis rates as a function of light level. 
Curves are based on photosynthesis irradiance model proposed by Webb et al. (1974), modified by Platt et al. (Date) to include a photoinhibition term. 
Units are not specific to allow script to be used for different photosynthesis metrics (i.e. carbon fixation, oxygen evolution, or electron transport rates)

Inputs:
light - irradiance data vector. 
pvar - photosynthesis rate measurements. data vector with length equal to light. 

Outputs: 
l - unique light levels used during the experiment. data vector
ps - average photosynthetic rate at each light level. data vector with equal length to l
ps_err - standard error of the measured photosynthetic rate at each light level. data vector with equal length to l
pmax - maximum photosynthetic rate. If photoinhibition is present, this is the theoretical maximum in the absence of photoinhibition. scalar. Derived from curve fit. 
ci_pmax - 95% confidence interval associated with pmax. 
ek - saturating light level. scalar. derived from curve fit. 
ci_ek - 95% confidence interval associated with ek
alpha - the initial light-limited slope/increase in photosynthesis with respect to light. scalar. derived from curve fit. 
ci_a - 95% confidence interval associated with alpha. 
R2 - R^2 of curve fit
mused - photosynthesis-irradiance model used (no photoinhibition model vs. photoinhibition model) the program automatically picks the model that provides the lowest R2
modelfit - matlab object containing fit parameters. 

