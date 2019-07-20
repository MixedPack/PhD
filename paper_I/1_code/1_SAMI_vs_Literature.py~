from scipy.optimize import minimize
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import polyfit
import math
import pandas as pd

#*****************************Input values**********************************************
Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table_DWARFS.csv')	
FCC_index = Galaxies['FCC']M_r = Galaxies['M_r(mag)']M_g = Galaxies['M_g(mag)']R_e = Galaxies['R_e(arcsec)'] * 0.09550 * 1000		#in pclog_stel_mass_r = Galaxies['log10(M_*/M_sun)']log_stel_mass_V = log_stel_mass_r - 0.4*0.16g_r = Galaxies['g'] - Galaxies['r']ERR_M_r = Galaxies['ERR_m_r(mag)']ERR_M_g = Galaxies['ERR_m_g(mag)']ERR_R_e = Galaxies['ERR_R_e(arcsec)']ERR_log_stMass_r = Galaxies['ERR_log10(M_*)']ERR_log_stMass_V = ERR_log_stMass_r
Kinematics = pd.read_csv('../2_pipeline/1_Kinematics_Table/Kinematics_Table_DWARFS.csv')
Dispersion = Kinematics['Sigma']
ERR_Disp = Kinematics['Error.Sigma']

#*****************************Calculation***********************************************L_sun = 3.848*math.pow(10,26) #wattM_r_sun = 4.64M_V_sun = 4.83Mass_sun = 2*math.pow(10,30) #kgM_V = -0.565*g_r - 0.016 + M_gERR_M_V = 0.001*g_r - 0.001 + ERR_M_gDyn_Mass_half = 930 * np.power(Dispersion,2) * R_e	#M_1/2 half-light ratio dynamical mass [M_Sun]Lumin_half_V = 0.5 * np.power(10, 0.4*(M_V_sun - M_V)) #L_1/2 half-light ration luminosity [L_Sun]Lumin_half_r = 0.5 * np.power(10, 0.4*(M_r_sun - M_r))ML_V = Dyn_Mass_half/Lumin_half_VML_r = Dyn_Mass_half/Lumin_half_rStell_Mass_r =  np.power(10,log_stel_mass_r)					#stellar mass in solar mass unitERR_M_dyn = 930 * np.sqrt(np.power(2*Dispersion*ERR_Disp*R_e,2)+np.power(np.power(Dispersion,2)*ERR_R_e,2))ERR_L_half_V = -0.4*Lumin_half_V*np.log(10)*ERR_M_VERR_L_half_r = -0.4*Lumin_half_r*np.log(10)*ERR_M_rERR_ML_V = (Dyn_Mass_half/Lumin_half_V) * np.sqrt(np.power(ERR_M_dyn/Dyn_Mass_half,2) + np.power(ERR_L_half_V/Lumin_half_V,2))ERR_ML_r = (Dyn_Mass_half/Lumin_half_r) * np.sqrt(np.power(ERR_M_dyn/Dyn_Mass_half,2) + np.power(ERR_L_half_r/Lumin_half_r,2))#*****************************Falcon-Barroso et al. 2011***********************************************Falcon11 = pd.read_csv('../0_data/Literature/FalconBarroso11_Spec.csv')F11_Disper = Falcon11['\sigma_e(km/s)']F11_ML_stel_V = Falcon11['M\L*_V(M/Lsun)']ERR_F11_ML_stel_V = abs(Falcon11['M\L*_V+ERR'] - Falcon11['M\L*_V-ERR'])/2Falcon11_2 = pd.read_csv('../0_data/Literature/FalconBarroso11_Phot.csv')F11_Re = 1000 * math.pi * Falcon11_2['R_{eV}(arcsec)'] * Falcon11_2['D(Mpc)'] / 648 #in pcF11_M_V = Falcon11_2['M_V(mag)']ERR_F11_Re = 1000 * math.pi * Falcon11_2['R_{eV}_ERR'] * Falcon11_2['D(Mpc)'] / 648ERR_F11_M_V = Falcon11_2['M_V_ERR']ERR_F11_Disper = 0 #for the time beingF11_DMass_half = 930 * np.power(F11_Disper,2) * F11_ReF11_Lumin_half_V = 0.5 * np.power(10, 0.4*(M_V_sun - F11_M_V))F11_ML_V = F11_DMass_half/F11_Lumin_half_VERR_F11_M_dyn = 930 * np.sqrt(np.power(2*F11_Disper*ERR_F11_Disper*F11_Re,2)+np.power(np.power(F11_Disper,2)*ERR_F11_Re,2))ERR_F11_Lum_V = -0.4*F11_Lumin_half_V*np.log(10)*ERR_F11_M_VERR_F11_ML_V = (F11_DMass_half/F11_Lumin_half_V) * np.sqrt(np.power(ERR_F11_M_dyn/F11_DMass_half,2) + np.power(ERR_F11_Lum_V/F11_Lumin_half_V,2))F11_log_stMass_V = np.log10(F11_ML_stel_V*(2*F11_Lumin_half_V))ERR_F11_log_stMass_V = np.sqrt(np.power(ERR_F11_ML_stel_V/(np.log(10)*F11_ML_stel_V),2)+np.power(ERR_F11_Lum_V/(np.log(10)*F11_Lumin_half_V),2))
#*****************************Toloba et al. 2014 (2 papers)**********************************Toloba14 = pd.read_csv('../0_data/Literature/Toloba14_Masses.csv')T14_log_DMass = Toloba14['log(M_e,dyn)[M_sun]']T14_log_StMass_r = Toloba14['log(M_e,*)[M_sun]']T14_log_StMass_V = T14_log_StMass_r - 0.4*0.16T14_ML_r = Toloba14['(M/L)_dyn,r[M_sun/L_sun,r]']T14_ML_V = T14_ML_r * (10**0.16)Toloba14_2 = pd.read_csv('../0_data/Literature/Toloba14_PhotKinem.csv')T14_M_r = Toloba14_2['M_r[mag]']T14_M_V = 0.16 + T14_M_r 	# Girardi et al. 2004 and Toloba et al. 2014T14_L_half_V = 0.5 * np.power(10, 0.4*(M_V_sun - T14_M_V))ERR_T14_log_DMass = Toloba14['log(M_e,dyn)_ERR']ERR_T14_log_StMass_r = Toloba14['log(M_e,*)_ERR']ERR_T14_log_StMass_V = ERR_T14_log_StMass_rERR_T14_log_ML_V = (Toloba14['(M/L)_dyn,r_ERR']* (10**0.16))/(np.log(10)*T14_ML_V)ERR_T14_L_half_V = 0 #for the time being
#*****************************Plotting**************************************************
Figure = plt.figure(figsize=(20,7))

ax = Figure.add_subplot(131)

kwargs = dict(histtype='stepfilled', alpha=0.3, density=True, bins=[-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13])
plt.hist(M_V, **kwargs)
plt.hist(F11_M_V, **kwargs)
plt.hist(T14_M_V, **kwargs)

ax = Figure.add_subplot(132)
kwargs = dict(histtype='stepfilled', alpha=0.3, density=True, bins=[7,8,9,10,11,12])
plt.hist(log_stel_mass_V, **kwargs)
plt.hist(F11_log_stMass_V , **kwargs)
plt.hist(T14_log_StMass_V, **kwargs)

Figure.savefig("../2_pipeline/2_SAMI_vs_Literature/SAMIvsLiter.pdf")

