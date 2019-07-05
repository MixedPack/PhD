from scipy.optimize import minimizeimport numpy as npfrom astropy.io import fitsimport matplotlib.pyplot as pltfrom numpy.polynomial.polynomial import polyfitimport mathimport pandas as pd#*****************************Input values**********************************************Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table.csv')	FCC_index = Galaxies['FCC']g = Galaxies['g']r = Galaxies['r']i = Galaxies['i']u = Galaxies['u']Kinematics = pd.read_csv('../2_pipeline/1_Kinematics_Table/Kinematics_Table.csv')Dispersion = Kinematics['Sigma']#*****************************3d Data***********************************************fig, axs = plt.subplots(4, 1, sharex=True, figsize=(10,15))
fig.subplots_adjust(hspace=0)  # Remove horizontal space between axes

axs[0].scatter(np.log10(Dispersion), g-r)
axs[0].set_ylabel('$g-r$ (mag)',fontsize=20)
axs[0].tick_params(axis='both', which='major', labelsize=18)
axs[1].scatter(np.log10(Dispersion), u-g)
axs[1].set_ylabel('$u-g (mag)$',fontsize=20)
axs[1].tick_params(axis='both', which='major', labelsize=18)
axs[2].scatter(np.log10(Dispersion), r-i)
axs[2].set_ylabel('$r-i$ (mag)',fontsize=20)
axs[2].tick_params(axis='both', which='major', labelsize=18)
axs[3].scatter(np.log10(Dispersion), i-g)
axs[3].set_ylabel('$i-g$ (mag)',fontsize=20)axs[3].set_xlabel('log($\sigma_{15})$ (km/s)',fontsize=20)
axs[3].tick_params(axis='both', which='major', labelsize=18)

#fig =  plt.figure()#ax = fig.add_axes([0.1, 0.7, 0.8, 0.85])#plt.scatter(Dispersion, g-r)#plt.ylabel('$g-r$')#ax = fig.add_axes([0.1, 0.5, 0.8, 0.65])#plt.scatter(Dispersion, u-g)#plt.ylabel('$u-g$')#ax = fig.add_axes([0.1, 0.3, 0.8, 0.45])#plt.scatter(Dispersion, r-i)#plt.ylabel('$r-i$')#ax = fig.add_axes([0.1, 0.1, 0.8, 0.25])#plt.scatter(Dispersion, i-g)plt.savefig("../2_pipeline/3_color_Disp/Color_Disp.pdf")