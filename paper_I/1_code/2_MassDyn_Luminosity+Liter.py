from scipy.optimize import minimize
#*****************************Input values**********************************************
if dwarf ==1 :
	Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table_DWARF.csv')	
	Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table_EXC.csv')
FCC_index = Galaxies['FCC']
	Kinematics = pd.read_csv('../2_pipeline/0_Kinematics_Table/Kinematics_Table_DWARF.csv')
	Kinematics = pd.read_csv('../2_pipeline/0_Kinematics_Table/Kinematics_Table_EXC.csv')
Dispersion = Kinematics['Sigma']
if dwarf == 0:
		log_L_r = np.log10(Lumin_half_r[l])
		log_Mdyn_r = np.log10(Dyn_Mass_half[l])
		ERR_log_Mdyn_r = ERR_M_dyn[l]/(Dyn_Mass_half[l]*np.log(10))	
		M_stelL = Stell_Mass_r[l]/Lumin_half_r[l]
		dict={'FCC':FCC_index[l],
		ListDict.append(dict.copy())

	Figure.savefig("../2_pipeline/2_MassDyn_Luminosity+Liter/DyM-L+Liter_DWARF.pdf", bbox_inches='tight')
	Figure.savefig("../2_pipeline/2_MassDyn_Luminosity+Liter/DyM-L+Liter.pdf", bbox_inches='tight')
	
bx.tick_params(axis='both', which='major', labelsize=18)


	Figure2.savefig("../2_pipeline/2_MassDyn_Luminosity+Liter/DyM_StM+Liter_DWARF.pdf")
	Figure2.savefig("../2_pipeline/2_MassDyn_Luminosity+Liter/DyM_StM+Liter.pdf")