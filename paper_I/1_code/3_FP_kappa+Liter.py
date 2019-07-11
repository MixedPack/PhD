from scipy.optimize import minimizeimport numpy as npfrom astropy.io import fitsimport matplotlib.pyplot as pltfrom numpy.polynomial.polynomial import polyfitimport mathimport pandas as pd#x_axis = "log_Sigma"x_axis = "log_Re"name = 1residuals = 0reverse = 1exclude = 1
#*****************************Input values**********************************************if exclude ==1:
	Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table_EXC.csv')	
else:
	Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table.csv')FCC_index = Galaxies['FCC']M_r = Galaxies['M_r(mag)']mu_r = Galaxies['mu_r(mag/arcsec2)'] 	#surface brightness calculated with Re in arcsecR_eff = Galaxies['R_e(arcsec)'] * 0.092 			#in kpcI_e = 10**(-0.4*(mu_r - 27))if exclude ==1:
	Kinematics = pd.read_csv('../2_pipeline/1_Kinematics_Table/Kinematics_Table_EXC.csv')else:
	Kinematics = pd.read_csv('../2_pipeline/1_Kinematics_Table/Kinematics_Table.csv')
Velocity = Kinematics['Velocity']Dispersion = Kinematics['Sigma']Vel_error = Kinematics['Error.Velocity']Disp_error = Kinematics['Error.Sigma']k1 = (np.log10(Dispersion**2) + np.log10(R_eff))/ np.sqrt(2)k2 = (np.log10(Dispersion**2) + 2*np.log10(I_e) - np.log10(R_eff))/ np.sqrt(6)k3 = (np.log10(Dispersion**2) - np.log10(I_e) - np.log10(R_eff))/ np.sqrt(3)#*****************************Toloba et al. (2011)**********************************#Toloba et all (2011)Toloba = pd.read_csv('../0_data/Literature/Toloba11.csv')TDisper = Toloba['Dispersion']TR_eff = Toloba['R_e'] * 0.0727 #in kpcTmu_V_e = Toloba['mu_V_e']ColB_V = Toloba['E(B-V)']Colg_V = 0.630*ColB_V - 0.124Colg_r = (-1 * Colg_V + 0.016)/(-0.565)Colr_V = Colg_V - Colg_rTmu_r_e = Tmu_V_e + Colr_VTI_e = 10**(-0.4*(Tmu_r_e - 27))Tk1 = (np.log10(TDisper**2) + np.log10(TR_eff))/ np.sqrt(2)Tk2 = (np.log10(TDisper**2) + 2*np.log10(TI_e) - np.log10(TR_eff))/ np.sqrt(6)Tk3 = (np.log10(TDisper**2) - np.log10(TI_e) - np.log10(TR_eff))/ np.sqrt(3)
#*****************************Falcon-Barroso et al. 2011***********************************************Falcon11 = pd.read_csv('../0_data/Literature/FalconBarroso11_Spec.csv')F11_Disper = Falcon11['\sigma_e(km/s)']
Falcon11_2 = pd.read_csv('../0_data/Literature/FalconBarroso11_Phot.csv')F11_R_eff = Falcon11_2['R_{eV}(arcsec)']  #in arcsecF11_m_V = Falcon11_2['M_V(mag)'] + 30.88F11_mu_V = Falcon11_2['\mu_{eV}(mag arcsec^{-2})']
F11_mu_r = F11_mu_V - 0.16
FI_e = 10**(-0.4*(F11_mu_r - 27))

Fk1 = (np.log10(F11_Disper**2) + np.log10(F11_R_eff))/ np.sqrt(2)Fk2 = (np.log10(F11_Disper**2) + 2*np.log10(FI_e) - np.log10(F11_R_eff))/ np.sqrt(6)Fk3 = (np.log10(F11_Disper**2) - np.log10(FI_e) - np.log10(F11_R_eff))/ np.sqrt(3)
#*****************************Plotting kappa space*********************************fig = plt.figure()s1 = fig.add_axes([0.1, 0.5, 0.8, 0.4])s1.plot(k1,k2,'b.')s1.plot(Tk1,Tk2,'y*')
s1.plot(Fk1,Fk2,'g+')plt.ylabel("$\kappa_2\propto log((M/L){I_e}^3)$",fontsize=12)plt.text(1.5, 4.2, "r_band",fontsize=12)s2 = fig.add_axes([0.1, 0.1, 0.8, 0.4])s2.plot(k1,k3,'b.')s2.plot(Tk1,Tk3,'y*')
s2.plot(Fk1,Fk3,'g+')plt.xlabel("$\kappa_1\propto log(M)$",fontsize=12)plt.ylabel("$\kappa_3\propto log(M/L)$",fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=10)
if exclude == 1:
	plt.savefig("../2_pipeline/3_FP_kappa+Liter/FP_kappa+Liter_EXC.pdf")else:
	plt.savefig("../2_pipeline/3_FP_kappa+Liter/FP_kappa+Liter.pdf")