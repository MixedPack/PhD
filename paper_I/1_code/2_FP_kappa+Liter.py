from scipy.optimize import minimizeimport numpy as npfrom astropy.io import fitsimport matplotlib.pyplot as pltfrom numpy.polynomial.polynomial import polyfitimport mathimport pandas as pdexclude = 0
dwarf = 1
F3D = 0#*****************************Input values**********************************************if exclude ==1:
	Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table_EXC.csv')	
elif dwarf == 1:
	Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table_DWARF.csv')		
else:
	Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table.csv')FCC_index = Galaxies['FCC']M_r = Galaxies['M_r(mag)']mu_r = Galaxies['mu_r(mag/arcsec2)'] 	#surface brightness calculated with Re in arcsecR_eff = Galaxies['R_e(arcsec)'] * 0.092 			#in kpcI_e = 10**(-0.4*(mu_r - 27))if exclude ==1:
	Kinematics = pd.read_csv('../2_pipeline/0_Kinematics_Table/Kinematics_Table_EXC.csv')elif dwarf == 1:
	Kinematics = pd.read_csv('../2_pipeline/0_Kinematics_Table/Kinematics_Table_DWARF.csv')
else:
	Kinematics = pd.read_csv('../2_pipeline/0_Kinematics_Table/Kinematics_Table.csv')
Velocity = Kinematics['Velocity']Dispersion = Kinematics['Sigma']Vel_error = Kinematics['Error.Velocity']Disp_error = Kinematics['Error.Sigma']

#*******Replacing Giant SAMIs with Fornax3D results*******
if F3D == 1:
	F3D = pd.read_csv('../0_data/Literature/F3D.csv')	Fobj = F3D['OBJ']	Fsigma = F3D['sigma_Re']	for j in range(len(FCC_index)):		for k in range(len(Fobj)):			if FCC_index[j] == Fobj[k]:				Dispersion[j] = Fsigma[k]
				k1 = (np.log10(Dispersion**2) + np.log10(R_eff))/ np.sqrt(2)k2 = (np.log10(Dispersion**2) + 2*np.log10(I_e) - np.log10(R_eff))/ np.sqrt(6)k3 = (np.log10(Dispersion**2) - np.log10(I_e) - np.log10(R_eff))/ np.sqrt(3)
#*****************************Toloba et al. (2011)**********************************#Toloba et all (2011)Toloba = pd.read_csv('../0_data/Literature/Toloba11.csv')TDisper = Toloba['Dispersion']TR_eff = Toloba['R_e'] * 0.0727 #in kpcTmu_V_e = Toloba['mu_V_e']ColB_V = Toloba['E(B-V)']Colg_V = 0.630*ColB_V - 0.124Colg_r = (-1 * Colg_V + 0.016)/(-0.565)Colr_V = Colg_V - Colg_rTmu_r_e = Tmu_V_e + Colr_VTI_e = 10**(-0.4*(Tmu_r_e - 27))Tk1 = (np.log10(TDisper**2) + np.log10(TR_eff))/ np.sqrt(2)Tk2 = (np.log10(TDisper**2) + 2*np.log10(TI_e) - np.log10(TR_eff))/ np.sqrt(6)Tk3 = (np.log10(TDisper**2) - np.log10(TI_e) - np.log10(TR_eff))/ np.sqrt(3)
#*****************************Falcon-Barroso et al. 2011***********************************************Falcon11 = pd.read_csv('../0_data/Literature/FalconBarroso11_Spec.csv')F11_Disper = Falcon11['\sigma_e(km/s)']
Falcon11_2 = pd.read_csv('../0_data/Literature/FalconBarroso11_Phot.csv')F11_R_eff = math.pi * Falcon11_2['R_{eV}(arcsec)'] * Falcon11_2['D(Mpc)'] / 648 #in kpcF11_m_V = Falcon11_2['M_V(mag)'] + 30.88F11_mu_V = Falcon11_2['\mu_{eV}(mag arcsec^{-2})']
F11_mu_r = F11_mu_V - 0.16
FI_e = 10**(-0.4*(F11_mu_r - 27))

Fk1 = (np.log10(F11_Disper**2) + np.log10(F11_R_eff))/ np.sqrt(2)Fk2 = (np.log10(F11_Disper**2) + 2*np.log10(FI_e) - np.log10(F11_R_eff))/ np.sqrt(6)Fk3 = (np.log10(F11_Disper**2) - np.log10(FI_e) - np.log10(F11_R_eff))/ np.sqrt(3)
#*****************************Plotting kappa space*********************************fig = plt.figure()s1 = fig.add_axes([0.2, 0.5, 0.7, 0.4])s1.plot(k1,k2,'b.')s1.plot(Tk1,Tk2,'y*')
s1.plot(Fk1,Fk2,'g+')plt.xticks((2,3))
plt.ylabel("$\kappa_2\propto log((M/L){I_e}^3)$",fontsize=12)plt.text(1.5, 4.2, "r_band",fontsize=12)s2 = fig.add_axes([0.2, 0.1, 0.7, 0.4])s2.plot(k1,k3,'b.')s2.plot(Tk1,Tk3,'y*')
s2.plot(Fk1,Fk3,'g+')
plt.yticks((-0.25,0.00,0.25,0.50))plt.xlabel("$\kappa_1\propto log(M)$",fontsize=12)plt.ylabel("$\kappa_3\propto log(M/L)$",fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=10)
if exclude == 1:
	plt.savefig("../2_pipeline/2_FP_kappa+Liter/FP_kappa+Liter_EXC.pdf")elif dwarf == 1:
	plt.savefig("../2_pipeline/2_FP_kappa+Liter/FP_kappa+Liter_DWARF.pdf")
else:
	plt.savefig("../2_pipeline/2_FP_kappa+Liter/FP_kappa+Liter.pdf")