from astropy.io import fitsfrom astropy.io import asciifrom plotbin import display_binsfrom plotbin import display_pixelsfrom plotbin import plot_velfieldfrom plotbin import display_bins_generatorsimport numpy as npimport matplotlib.pyplot as pltimport pandas as pdimport mathimport csvimport matplotlib.gridspec as gridspec#************************FUNCTIONS******************************def import_file(index,type=None):		folder = "../../SAMI/Blue_Cubes/results_values/Binned/Simple/Velocity/"	file = folder + "fcc" + str(int(index)) + "_blue.fits"	hdulist = fits.open(file)	table = hdulist[1].data	binned_output_V = {"XBIN": table["XBIN"], "YBIN": table["YBIN"], "FLUXBIN": table["FLUXBIN"], 	           "VPXF": table["VPXF"], "SPXF": table["SPXF"]}		folder = "../../SAMI/Blue_Cubes/results_values/Binned/Simple/Velocity/"	file = folder + "fcc" + str(int(index)) + "_blue_unbinned.fits"	hdulist = fits.open(file)	table = hdulist[1].data	unbinned_output_V = {"X": table["X"], "Y": table["Y"], "FLUX": table["FLUX"]}		folder = "../../SAMI/Blue_Cubes/results_values/Binned/Simple/Dispersion/"	file = folder + "fcc" + str(int(index)) + "_blue.fits"	hdulist = fits.open(file)	table = hdulist[1].data	binned_output_S = {"XBIN": table["XBIN"], "YBIN": table["YBIN"], "FLUXBIN": table["FLUXBIN"], 	           "VPXF": table["VPXF"], "SPXF": table["SPXF"]}		folder = "../../SAMI/Blue_Cubes/results_values/Binned/Simple/Dispersion/"	file = folder + "fcc" + str(int(index)) + "_blue_unbinned.fits"	hdulist = fits.open(file)	table = hdulist[1].data	unbinned_output_S = {"X": table["X"], "Y": table["Y"], "FLUX": table["FLUX"]}		if type == None or type == "V":		return (unbinned_output_V, binned_output_V)	elif type == "S":		return (unbinned_output_S, binned_output_S)def get_MinMax_v1(galaxy_index,type=None):		range = ascii.read('./Min_Max.txt')		ID = np.where(range["FCC"]==galaxy_index)[0]	vmin = range["VMIN"][ID]	vmax = range["VMAX"][ID]	smin = range["SMIN"][ID]	smax = range["SMAX"][ID]	if type == None or type == "V":		return (vmin, vmax)	elif type == "S":		return (smin, smax)def get_MinMax_v2(galaxy_index,type=None):		if type == None or type == "V":		data_dict = import_file(galaxy_index, type="V")		binned_data = data_dict[1] 		Bin_Velocity = binned_data["VPXF"] - np.median(binned_data["VPXF"])		RangeMin = np.median(Bin_Velocity)-2*np.std(Bin_Velocity)		RangeMax = np.median(Bin_Velocity)+2*np.std(Bin_Velocity)		Good_Velocity = [Bin_Velocity[i] for i in range(len(Bin_Velocity)) if Bin_Velocity[i] >= RangeMin and Bin_Velocity[i] <= RangeMax]		return (np.median(Good_Velocity)-2*np.std(Good_Velocity), np.median(Good_Velocity)+2*np.std(Good_Velocity))	if type == "S":		data_dict = import_file(galaxy_index, type="S")		binned_data = data_dict[1] 		Bin_Dispersion = binned_data["SPXF"]		RangeMin = np.median(Bin_Dispersion)-2.5*np.std(Bin_Dispersion)		RangeMax = np.median(Bin_Dispersion)+2.5*np.std(Bin_Dispersion)		Good_Dispersion = [Bin_Dispersion[i] for i in range(len(Bin_Dispersion)) if Bin_Dispersion[i] >= RangeMin and Bin_Dispersion[i] <= RangeMax]		return (np.median(Good_Dispersion)-2.5*np.std(Good_Dispersion), np.median(Good_Dispersion)+2.5*np.std(Good_Dispersion))	def velocity_map(galaxy_index):		data_dict = import_file(galaxy_index, type="V")	unbinned_data = data_dict[0]	binned_data = data_dict[1]	VMIN = get_MinMax_v2(galaxy_index,type="V")[0]	VMAX = get_MinMax_v2(galaxy_index,type="V")[1]	ax = display_bins_generators.display_bins_generators(binned_data["XBIN"], binned_data["YBIN"], 	                                     binned_data["VPXF"]-np.median(binned_data["VPXF"]), unbinned_data["X"], unbinned_data["Y"],	                                     vmin=VMIN, vmax=VMAX)	return axdef dispersion_map(galaxy_index):		data_dict = import_file(galaxy_index, type="S")	unbinned_data = data_dict[0]	binned_data = data_dict[1]	SMIN = get_MinMax_v2(galaxy_index,type="S")[0]		#np.mean(binned_data ["SPXF"])-2*np.std(binned_data ["SPXF"])	SMAX = get_MinMax_v2(galaxy_index,type="S")[1]	ax = display_bins_generators.display_bins_generators(binned_data["XBIN"], binned_data["YBIN"], 	                                   binned_data["SPXF"], unbinned_data["X"], unbinned_data["Y"], 	                                   vmin=SMIN, vmax=SMAX)	return ax#************************MAIN_BODY******************************Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table.csv')FCC_index = Galaxies['FCC']total_rows = math.ceil(len(FCC_index)/5)figureV = plt.figure(2, figsize=(20,40))gs = gridspec.GridSpec(10*2-1,5*2-1, height_ratios=[1, 0.4, 1, 0.4, 1, 0.4, 1, 0.4, 1, 0.4, 1, 0.4, 1, 0.4, 1, 0.4, 1, 0.4, 1], 						width_ratios=[1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1])#,                       #left=0.05, right=0.95, bottom=0.08, top=0.93, wspace=0.02, hspace=0.03)for n in range(len(FCC_index)):	R = 2 * (n//5)	C = 2 * round(5*(n/5-n//5))	ax = plt.subplot(gs[R,C])	ax = velocity_map(FCC_index[n])	plt.title("FCC"+str(int(FCC_index[n])), fontsize=8)#,y=1.15)figureV.savefig("./Velocity_map.pdf")figureS = plt.figure(2, figsize=(20,40))gs = gridspec.GridSpec(10*2-1,5*2-1, height_ratios=[1, 0.4, 1, 0.4, 1, 0.4, 1, 0.4, 1, 0.4, 1, 0.4, 1, 0.4, 1, 0.4, 1, 0.4, 1], 						width_ratios=[1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1])for n in range(len(FCC_index)):	R = 2 * (n//5)	C = 2 * round(5*(n/5-n//5))	ax = plt.subplot(gs[R,C])	ax = dispersion_map(FCC_index[n])	plt.title("FCC"+str(int(FCC_index[n])), fontsize=8)figureS.savefig("./Dispersion_map.pdf")