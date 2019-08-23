import csvfrom astropy.io import fitsimport numpy as npfrom pathlib import Pathfrom astropy.io import asciiimport matplotlib.pyplot as plt
#%precision 2	#print(GalaxyList[0].keys()) #the columns' namesdef import_sigma(index,ratio):		dictS = {}	dictE = {}
	fileS = "../0_data/SAMI/Diff_Apertures/Re_" + str(int(ratio)) + "/fcc" + str(int(index)) + "_blue_results_simple_onebin_" + str(int(ratio)) + ".fits"	fileE = "../0_data/SAMI/Diff_Apertures/Re_" + str(int(ratio)) + "/fcc" + str(int(index)) + "_blue_results_mc_onebin_" + str(int(ratio)) + ".fits"	my_file = Path(fileS)	if my_file.exists():
		hdulist = fits.open(fileS)		table = hdulist[1].data		sigma = np.mean(table["SPXF"])		dictS.update({"Sigma/"+str(int(ratio)): "%.4f" % sigma})	else:		dictS.update({"Sigma/"+str(int(ratio)): np.NaN})			my_file = Path(fileE)	if my_file.exists():		hdulist = fits.open(fileE)		table = hdulist[1].data		Serror= np.mean(table['DSp'])		dictE.update({"Error.Sigma/"+str(int(ratio)): "%.4f" % Serror})	else:		dictE.update({"Error.Sigma/"+str(int(ratio)): np.NaN})			return (dictS, dictE)with open('../2_pipeline/0_Galaxies_Table/Galaxies_Table_DWARF.csv') as csvfile:	listdict_info = list(csv.DictReader(csvfile))ListDictGal = []for d in listdict_info:	dictgal = {'FCC':d['FCC'], 'Re(arcsec)':d['R_e(arcsec)'], 'mu_r':d['mu_r(mag/arcsec2)']}	dictgal.update(import_sigma(d['FCC'],1)[0])	dictgal.update(import_sigma(d['FCC'],1)[1])
	dictgal.update(import_sigma(d['FCC'],2)[0])	dictgal.update(import_sigma(d['FCC'],2)[1])
	dictgal.update(import_sigma(d['FCC'],4)[0])	dictgal.update(import_sigma(d['FCC'],4)[1])
	dictgal.update(import_sigma(d['FCC'],6)[0])	dictgal.update(import_sigma(d['FCC'],6)[1])	ListDictGal.append(dictgal)keys = ListDictGal[0].keys()
with open('../2_pipeline/1_Kinematics_DiffApert/Kinematics_DiffApert_DWARF.csv', 'w') as output_file:    dict_writer = csv.DictWriter(output_file, keys)    dict_writer.writeheader()    dict_writer.writerows(ListDictGal)

#*******************************Plotting***************************************#

Figure = plt.figure(figsize=(15, 30))
ax = Figure.add_subplot(111)

for d in listdict_info:
	x=[]; y=[]; Ey=[]	
	for i in [6,4,2,1]:		
		y_dummy = import_sigma(d['FCC'],i)[0]
		Ey_dummy = import_sigma(d['FCC'],i)[1]
		x_dummy = float(d['R_e(arcsec)'])/i	
		key = y_dummy.keys()		
		yy = y_dummy.values(); y.append(float(yy[0]))
		Eyy = Ey_dummy.values(); Ey.append(float(Eyy[0]))
		x.append(x_dummy)	
	ax.plot(x, y)
	plt.annotate(str(d['FCC']), (x[0], y[0]), fontsize=8)
	ax.errorbar(x, y, xerr=0, yerr=Ey, fmt='o', ecolor='grey')
			
ax.set_xlabel("radius (arcsec)")
ax.set_ylabel("Aperture Sigma")

Figure.savefig("../2_pipeline/1_Kinematics_DiffApert/Kinematics_DiffApert_DWARF.pdf")
