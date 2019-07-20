import numpy as npimport matplotlib.pyplot as pltfrom astropy.io import fitsimport csvimport pyfitsimport mathfrom numpy.polynomial.polynomial import polyfitfrom mpl_toolkits.axes_grid1 import make_axes_locatableimport pandas as pd

exclude = 1#*************************SAMI Input values*********************************************if exclude == 1 :
	Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table_EXC.csv')	else:
	Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table.csv')
FCC_index = Galaxies['FCC']M_r = Galaxies['M_r(mag)']mu_r = Galaxies['mu_r(mag/arcsec2)'] 	#surface brightness calculated with Re in arcsecR_eff = Galaxies['R_e(arcsec)']if exclude == 1 :
	Kinematics = pd.read_csv('../2_pipeline/0_Kinematics_Table/Kinematics_Table_EXC.csv')else:
	Kinematics = pd.read_csv('../2_pipeline/0_Kinematics_Table/Kinematics_Table.csv')
Dispersion = Kinematics['Sigma']Disp_error = Kinematics['Error.Sigma']#*****************************Toloba et al. (2011)**********************************#Toloba et all (2011)Toloba = pd.read_csv('../0_data/Literature/Toloba11.csv')TDisper = Toloba['Dispersion']TR_eff = Toloba['R_e'] * 0.0727 #in kpcTmu_V_e = Toloba['mu_V_e']ColB_V = Toloba['E(B-V)']TM_V = Toloba['M_V']Colg_V = 0.630*ColB_V - 0.124Colg_r = (-1 * Colg_V + 0.016)/(-0.565)Colr_V = Colg_V - Colg_rTM_r = TM_V + Colr_VTmu_r_e = Tmu_V_e + Colr_V
#*****************************Falcon-Barroso et al. 2011***********************************************Falcon11 = pd.read_csv('../0_data/Literature/FalconBarroso11_Spec.csv')F11_Disper = Falcon11['\sigma_e(km/s)']
Falcon11_2 = pd.read_csv('../0_data/Literature/FalconBarroso11_Phot.csv')F11_Re = math.pi * Falcon11_2['R_{eV}(arcsec)'] * Falcon11_2['D(Mpc)'] / 648 *1000#in kpcF11_M_V = Falcon11_2['M_V(mag)']F11_mu_V = Falcon11_2['\mu_{eV}(mag arcsec^{-2})']
ERR_F11_Re = math.pi * Falcon11_2['R_{eV}_ERR'] * Falcon11_2['D(Mpc)'] / 648 *1000ERR_F11_M_V = Falcon11_2['M_V_ERR']ERR_F11_Disper = 0 #for the time being
ERR_F11_mu_V = Falcon11_2['\mu_{eV}_ERR']
F11_M_r = F11_M_V - 0.16F11_mu_r = F11_mu_V - 0.16#****************Plotting Faber_Jackson relation****************************************cm = plt.cm.get_cmap('spring')x = M_ry = np.log10(Dispersion)y_err = (Disp_error)/(Dispersion*np.log(10))b, m = polyfit(x, y, 1)fitline = b + m * xTx = TM_rTy = np.log10(TDisper)Fx = F11_M_r
Fy = np.log10(F11_Disper)
fig = plt.figure(figsize = (20,7))ax = fig.add_subplot(1, 1, 1)plt.errorbar(x, y, yerr=y_err, fmt='.', ecolor='green')plt.scatter(x, y, s=R_eff*5, c=mu_r, cmap=cm)cbar = plt.colorbar()cbar.ax.invert_yaxis()cbar.set_label('$\mu_e$', rotation=-270, fontsize=16)cbar.ax.tick_params(labelsize=16)
plt.plot(Tx, Ty,'b*')
plt.plot(Fx, Fy, 'g+')plt.plot(x, fitline, c='grey')textstr = 'log($\sigma_{15}$) = ' + str("%.2f" % m) + ' $M_r$ ' + str("%.2f" % b)plt.title(textstr, y=1.05, fontsize=18)major_ticks = np.arange(-14,-23,-1)ax.set_xticks(major_ticks)#for i in range(len(FCC_index)):#	plt.annotate(str(int(FCC_index[i])), (x[i], y[i]), fontsize=8)plt.ylabel('log($\sigma_{15}$) (km/s)', fontsize=18)plt.ylim(-0.5, 3.5)plt.xlabel('$M_r$ (mag)', fontsize=18)plt.tick_params(axis='both', which='major', labelsize=16)
if exclude == 1:
	fig.savefig('../2_pipeline/2_FJ+Liter/Faber_Jackson+Liter_EXC.pdf')else:
	fig.savefig('../2_pipeline/2_FJ+Liter/Faber_Jackson+Liter.pdf')