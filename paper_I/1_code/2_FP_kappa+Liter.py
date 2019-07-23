from scipy.optimize import minimize
dwarf = 1
F3D = 0
	Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table_EXC.csv')	
elif dwarf == 1:
	Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table_DWARF.csv')		
else:
	Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table.csv')
	Kinematics = pd.read_csv('../2_pipeline/0_Kinematics_Table/Kinematics_Table_EXC.csv')
	Kinematics = pd.read_csv('../2_pipeline/0_Kinematics_Table/Kinematics_Table_DWARF.csv')
else:
	Kinematics = pd.read_csv('../2_pipeline/0_Kinematics_Table/Kinematics_Table.csv')
Velocity = Kinematics['Velocity']

#*******Replacing Giant SAMIs with Fornax3D results*******
if F3D == 1:
	F3D = pd.read_csv('../0_data/Literature/F3D.csv')
				

#*****************************Falcon-Barroso et al. 2011***********************************************
Falcon11_2 = pd.read_csv('../0_data/Literature/FalconBarroso11_Phot.csv')
F11_mu_r = F11_mu_V - 0.16
FI_e = 10**(-0.4*(F11_mu_r - 27))

Fk1 = (np.log10(F11_Disper**2) + np.log10(F11_R_eff))/ np.sqrt(2)

s1.plot(Fk1,Fk2,'g+')
plt.ylabel("$\kappa_2\propto log((M/L){I_e}^3)$",fontsize=12)
s2.plot(Fk1,Fk3,'g+')
plt.yticks((-0.25,0.00,0.25,0.50))
plt.tick_params(axis='both', which='major', labelsize=10)

	plt.savefig("../2_pipeline/2_FP_kappa+Liter/FP_kappa+Liter_EXC.pdf")
	plt.savefig("../2_pipeline/2_FP_kappa+Liter/FP_kappa+Liter_DWARF.pdf")
else:
	plt.savefig("../2_pipeline/2_FP_kappa+Liter/FP_kappa+Liter.pdf")