
x_axis = "sigma"
name = 0
exclude = 1

#*************************SAMI Input values*********************************************
	Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table.csv')	
FCC_index = Galaxies['FCC']
	Kinematics = pd.read_csv('../2_pipeline/0_Kinematics_Table/Kinematics_Table.csv')
Dispersion = Kinematics['Sigma']
	y = M_r
if x_axis == "M_r":
	x = M_r
	plt.errorbar(x, y, xerr=x_err, fmt='.', ecolor='green')
	plt.errorbar(x, y, yerr=y_err, fmt='.', ecolor='green')
plt.scatter(x, y, s=R_eff*5, c=mu_r, cmap=cm)
if x_axis == "sigma":
	plt.gca().invert_yaxis()
if x_axis == "M_r":
	plt.gca().invert_xaxis() 
plt.plot(x, fitline, c='grey')
	for i in range(len(FCC_index)):	
if x_axis == "sigma":
	textstr = ' $M_r$ =' + str("%.2f" % m) + 'log($\sigma_{15}$)' + str("%.2f" % b)
if x_axis == "M_r":
	textstr = 'log($\sigma_{15}$) = ' + str("%.2f" % m) + ' $M_r$ ' + str("%.2f" % b)
	
plt.tick_params(axis='both', which='major', labelsize=16)
		fig.savefig('../2_pipeline/1_FJ/Faber_Jackson_x-axis=sigma_EXC.pdf')
		fig.savefig('../2_pipeline/1_FJ/Faber_Jackson_x-axis=M_r_EXC.pdf')
	if x_axis == "sigma":
		fig.savefig('../2_pipeline/1_FJ/Faber_Jackson_x-axis=sigma.pdf')
		fig.savefig('../2_pipeline/1_FJ/Faber_Jackson_x-axis=M_r.pdf')
#************************For test - FJ inputs: sigma************************************