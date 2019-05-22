from scipy.optimize import minimize
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import polyfit
import math
import pandas as pd

#x_axis = "log_Sigma"
x_axis = "log_Re"
name = 1
residuals = 0
reverse = 1

#*****************************Input values**********************************************
Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table.csv')	
FCC_index = Galaxies['FCC']
M_r = Galaxies['M_r(mag)']
mu_r = Galaxies['mu_r(mag/arcsec2)'] 	#surface brightness calculated with Re in arcsec
R_eff = Galaxies['R_e(arcsec)'] * 0.092 			#in kpc
I_e = 10**(-0.4*(mu_r - 27))

Kinematics = pd.read_csv('../2_pipeline/1_Kinematics_Table/Kinematics_Table.csv')
Velocity = Kinematics['Velocity']
Dispersion = Kinematics['Sigma']
Vel_error = Kinematics['Error.Velocity']
Disp_error = Kinematics['Error.Sigma']

k1 = (np.log10(Dispersion**2) + np.log10(R_eff))/ np.sqrt(2)
k2 = (np.log10(Dispersion**2) + 2*np.log10(I_e) - np.log10(R_eff))/ np.sqrt(6)
k3 = (np.log10(Dispersion**2) - np.log10(I_e) - np.log10(R_eff))/ np.sqrt(3)

#*****************************Toloba et al. (2011)**********************************
#Toloba et all (2011)
Toloba = pd.read_csv('../0_data/Literature/Toloba11.csv')
TDisper = Toloba['Dispersion']
TR_eff = Toloba['R_e'] * 0.0727 #in kpc
Tmu_V_e = Toloba['mu_V_e']
ColB_V = Toloba['E(B-V)']

Colg_V = 0.630*ColB_V - 0.124
Colg_r = (-1 * Colg_V + 0.016)/(-0.565)
Colr_V = Colg_V - Colg_r
Tmu_r_e = Tmu_V_e + Colr_V
TI_e = 10**(-0.4*(Tmu_r_e - 27))
l
Tk1 = (np.log10(TDisper**2) + np.log10(TR_eff))/ np.sqrt(2)
Tk2 = (np.log10(TDisper**2) + 2*np.log10(TI_e) - np.log10(TR_eff))/ np.sqrt(6)
Tk3 = (np.log10(TDisper**2) - np.log10(TI_e) - np.log10(TR_eff))/ np.sqrt(3)

#*****************************Plotting kappa space*********************************
fig = plt.figure()
s1 = fig.add_axes([0.1, 0.5, 0.8, 0.4])
s1.plot(k1,k2,'b.')
s1.plot(Tk1,Tk2,'gx')
plt.ylabel("$\kappa_2$")
plt.text(1.5, 4.2, "r_band")
s2 = fig.add_axes([0.1, 0.1, 0.8, 0.4])
s2.plot(k1,k3,'b.')
s2.plot(Tk1,Tk3,'gx')
plt.xlabel("$\kappa_1$")
plt.ylabel("$\kappa_3$")
plt.savefig("../2_pipeline/3_FP_kappa+Toloba/FP_kappa+Toloba.pdf")
