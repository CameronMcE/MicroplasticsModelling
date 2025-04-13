import numpy as np

source_data = np.load("emissions_UKESM.npy") #path to source emissions, shape is [month, size bin, lat, lon]

#Using 10-25 micron bin to extrapolate down with the chosen alpha value

#Constants
x_ub_ref = 25 #range of source data size bins
x_lb_ref = 10
alpha = 1.81 #from Leusch et al 2023, for airborne microplastics. Within uncertainty limits to match higher bins. 

#geometric means for bins in nanometres
gm_aitken = np.sqrt(5*50) #from UM size bins Aitken mode
gm_accum = np.sqrt(50*500) #from UM size bins Accum mode 
gm_coarse = np.sqrt(500*2500) #from UM size bins Coarse mode
gm_super = np.sqrt(2500*250000) #lower lim is UM size bin for super-coarse mode, upper lim is emissions data upper limit
gm_base = np.sqrt(10000*25000) #from Nikos' data of size bin to extrapolate

#radius in m
r_aitken = (gm_aitken/2)/1e9
r_accum = (gm_accum/2)/1e9
r_coarse = (gm_coarse/2)/1e9
r_super = (gm_super/2)/1e9
r_base = (gm_base/2)/1e9

#mass in kg
rho = 1000 #density of MPs, kg/m3
vol_aitken = (4/3)*np.pi*np.power(r_aitken,3.0); mass_aitken = (rho * vol_aitken) 
vol_accum = (4/3)*np.pi*np.power(r_accum,3.0); mass_accum = (rho * vol_accum) 
vol_coarse = (4/3)*np.pi*np.power(r_coarse,3.0); mass_coarse = (rho * vol_coarse) 
vol_super = (4/3)*np.pi*np.power(r_super,3.0); mass_super = (rho * vol_super)
vol_base = (4/3)*np.pi*np.power(r_base,3.0); mass_base = (rho * vol_base)

#size range limits of prediced MPS (aitken, accumulation, coarse)
lb_pred = np.array((0.005,0.05,0.5))
ub_pred = np.array((0.05,0.5,2.5))

#converting to number flux number = mass flux / mass of particle
C = source_data[:,1,:,:]/ mass_base #selecting the 10 - 25 micron bin

#applying the power law from Leusch et al. (2022)
i=0; x_lb_pred = lb_pred[i]; x_ub_pred = ub_pred[i]
n_pred_aitken = C * ((x_ub_pred - x_lb_pred)/(x_ub_ref - x_lb_ref)) * np.power(((x_ub_pred * x_lb_pred)/(x_ub_ref * x_lb_ref)),(-alpha/2.0))

i=1; x_lb_pred = lb_pred[i]; x_ub_pred = ub_pred[i]
n_pred_accum = C * ((x_ub_pred - x_lb_pred)/(x_ub_ref - x_lb_ref)) * np.power(((x_ub_pred * x_lb_pred)/(x_ub_ref * x_lb_ref)),(-alpha/2.0))

i=2; x_lb_pred = lb_pred[i]; x_ub_pred = ub_pred[i]
n_pred_coarse = C * ((x_ub_pred - x_lb_pred)/(x_ub_ref - x_lb_ref)) * np.power(((x_ub_pred * x_lb_pred)/(x_ub_ref * x_lb_ref)),(-alpha/2.0))

#Convert number fluxes to mass fluxes in kg/m2/s using mass of microplastics
aitken_mass = n_pred_aitken * mass_aitken
accum_mass = n_pred_accum*mass_accum 
coarse_mass = n_pred_coarse * mass_coarse
