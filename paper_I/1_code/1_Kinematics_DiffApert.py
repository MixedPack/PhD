import csv

	fileS = "../0_data/SAMI/Diff_Apertures/Re_" + str(int(ratio)) + "/fcc" + str(int(index)) + "_blue_results_simple_onebin_" + str(int(ratio)) + ".fits"
		hdulist = fits.open(fileS)
	dictgal.update(import_sigma(d['FCC'],2)[0])
	dictgal.update(import_sigma(d['FCC'],4)[0])
	dictgal.update(import_sigma(d['FCC'],6)[0])
with open('../2_pipeline/1_Kinematics_DiffApert/Kinematics_DiffApert_DWARF.csv', 'w') as output_file:

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