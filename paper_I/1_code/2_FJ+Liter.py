

exclude = 1
	Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table_EXC.csv')	
	Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table.csv')
FCC_index = Galaxies['FCC']
	Kinematics = pd.read_csv('../2_pipeline/0_Kinematics_Table/Kinematics_Table_EXC.csv')
	Kinematics = pd.read_csv('../2_pipeline/0_Kinematics_Table/Kinematics_Table.csv')
Dispersion = Kinematics['Sigma']
#*****************************Falcon-Barroso et al. 2011***********************************************
Falcon11_2 = pd.read_csv('../0_data/Literature/FalconBarroso11_Phot.csv')
ERR_F11_Re = math.pi * Falcon11_2['R_{eV}_ERR'] * Falcon11_2['D(Mpc)'] / 648 *1000
ERR_F11_mu_V = Falcon11_2['\mu_{eV}_ERR']
F11_M_r = F11_M_V - 0.16
Fy = np.log10(F11_Disper)

plt.plot(Tx, Ty,'b*')
plt.plot(Fx, Fy, 'g+')
if exclude == 1:
	fig.savefig('../2_pipeline/2_FJ+Liter/Faber_Jackson+Liter_EXC.pdf')
	fig.savefig('../2_pipeline/2_FJ+Liter/Faber_Jackson+Liter.pdf')