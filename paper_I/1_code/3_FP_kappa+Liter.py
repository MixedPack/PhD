from scipy.optimize import minimize

	Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table_EXC.csv')	
else:
	Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table.csv')
	Kinematics = pd.read_csv('../2_pipeline/1_Kinematics_Table/Kinematics_Table_EXC.csv')
	Kinematics = pd.read_csv('../2_pipeline/1_Kinematics_Table/Kinematics_Table.csv')
Velocity = Kinematics['Velocity']
#*****************************Falcon-Barroso et al. 2011***********************************************
Falcon11_2 = pd.read_csv('../0_data/Literature/FalconBarroso11_Phot.csv')
F11_mu_r = F11_mu_V - 0.16
FI_e = 10**(-0.4*(F11_mu_r - 27))

Fk1 = (np.log10(F11_Disper**2) + np.log10(F11_R_eff))/ np.sqrt(2)

s1.plot(Fk1,Fk2,'g+')
s2.plot(Fk1,Fk3,'g+')
plt.tick_params(axis='both', which='major', labelsize=10)

	plt.savefig("../2_pipeline/3_FP_kappa+Liter/FP_kappa+Liter_EXC.pdf")
	plt.savefig("../2_pipeline/3_FP_kappa+Liter/FP_kappa+Liter.pdf")