from scipy.optimize import minimize
reverse_axes = 1

mu_r = Galaxies['mu_r(mag/arcsec2)']
T14_M_r = Toloba14['M_r[mag]']
ERR_T14_Disp = Toloba14['sigma_e_ERR']
ERR_T14_R_e = 0 # for the time being
#Toloba14_2 = pd.read_csv('../0_data/Literature/Toloba14_kinematics.csv')
#T14_epsilon = Toloba14_2['epsilon_e']
T14_mu_r = Toloba14['M_r[mag]'] + 31.08 + 2.5 * np.log10(math.pi*np.power(Toloba14['R_e-r[arcsec]'],2)) + 2.5 * np.log10(2)
ERR_T14_mu_r = ERR_T14_M_r + (5 / np.log(10)) * (ERR_T14_R_e / Toloba14['R_e-r[arcsec]'])

Falcon11_2 = pd.read_csv('../0_data/Literature/FalconBarroso11_Phot.csv')
ERR_F11_Re = math.pi * Falcon11_2['R_{eV}_ERR'] * Falcon11_2['D(Mpc)'] / 648
ERR_F11_mu_V = Falcon11_2['\mu_{eV}_ERR']
F11_M_r = F11_M_V - 0.16

	plt.plot(Y, fitline, c='grey')
		plt.plot(np.log10(T14_Disp), output.x[0]*T14_mu_r+output.x[1]*np.log10(T14_R_e)+output.x[2], 'r*', label= 'Toloba et al. 2014')
		plt.plot(np.log10(F11_Disper), output.x[0]*F11_mu_r+output.x[1]*np.log10(F11_Re)+output.x[2],'gX', label= 'Falcon-Barroso et al. 2011')
		plt.plot(np.log10(T14_R_e), output.x[0]*T14_mu_r+output.x[1]*np.log10(T14_Disp)+output.x[2], 'r*', label= 'Toloba et al. 2014')
		plt.plot(np.log10(F11_Re), output.x[0]*F11_mu_r+output.x[1]*np.log10(F11_Disper)+output.x[2],'gX', label= 'Falcon-Barroso et al. 2011')
	elif x_axis == "log_Re":
	if x_axis == "log_Sigma":

else:
	plt.plot(fitline, Y, c='grey')
		plt.plot(output.x[0]*T14_mu_r+output.x[1]*np.log10(T14_R_e)+output.x[2], np.log10(T14_Disp), 'r*', label= 'Toloba et al. 2014')
		plt.plot(output.x[0]*F11_mu_r+output.x[1]*np.log10(F11_Re)+output.x[2], np.log10(F11_Disper), 'gX', label= 'Falcon-Barroso et al. 2011')
		plt.plot(output.x[0]*T14_mu_r+output.x[1]*np.log10(T14_Disp)+output.x[2], np.log10(T14_R_e), 'r*', label= 'Toloba et al. 2014')
		plt.plot(output.x[0]*F11_mu_r+output.x[1]*np.log10(F11_Disper)+output.x[2],np.log10(F11_Re), 'gX', label= 'Falcon-Barroso et al. 2011')
		plt.ylim(0.5, 3.0)	
	elif x_axis == "log_Re":
	if x_axis == "log_Sigma":


#********************************Residuals********************************************
		return np.abs(Y[i] - fit(X[i], output.x))/np.sqrt(output.x[0]**2+output.x[1]**2+1)
	
	
	
	