from astropy.io import fitsfrom astropy.io import asciifrom plotbin import display_binsfrom plotbin import display_pixelsfrom plotbin import plot_velfieldfrom plotbin import display_bins_generatorsimport numpy as npimport matplotlib.pyplot as pltimport pandas as pdimport mathimport csvimport matplotlib.gridspec as gridspecfrom pathlib import Path
from astropy.wcs import WCS
from astropy.visualization import (ZScaleInterval, ImageNormalize, LinearStretch)
#************************FUNCTIONS******************************def import_file(index,type=None):		file = "../0_data/SAMI/fcc" + str(int(index)) + "_blue_results_simple_binned(V).fits"	hdulist = fits.open(file)	table = hdulist[1].data	binned_output_V = {"XBIN": table["XBIN"], "YBIN": table["YBIN"], "FLUXBIN": table["FLUXBIN"], 	           "VPXF": table["VPXF"], "SPXF": table["SPXF"]}		file = "../0_data/SAMI/fcc" + str(int(index)) + "_blue_unbinned.fits"	hdulist = fits.open(file)	table = hdulist[1].data	unbinned_output = {"X": table["X"], "Y": table["Y"], "FLUX": table["FLUX"]}		file = "../0_data/SAMI/fcc" + str(int(index)) + "_blue_results_simple_binned(S).fits"	hdulist = fits.open(file)	table = hdulist[1].data	binned_output_S = {"XBIN": table["XBIN"], "YBIN": table["YBIN"], "FLUXBIN": table["FLUXBIN"], 	           "VPXF": table["VPXF"], "SPXF": table["SPXF"]}		if type == None or type == "V":		return (unbinned_output, binned_output_V)	elif type == "S":		return (unbinned_output, binned_output_S)
def get_MinMax_v1(galaxy_index,type=None):		range = ascii.read('../0_data/tmp/V&S_MinMax.txt')		ID = np.where(range["FCC"]==galaxy_index)[0]	vmin = range["VMIN"][ID]	vmax = range["VMAX"][ID]	smin = range["SMIN"][ID]	smax = range["SMAX"][ID]	if type == None or type == "V":		return (vmin, vmax)	elif type == "S":		return (smin, smax)
def get_MinMax_v2(galaxy_index,type=None):		if type == None or type == "V":		data_dict = import_file(galaxy_index, type="V")		binned_data = data_dict[1] 		Bin_Velocity = binned_data["VPXF"] - np.median(binned_data["VPXF"])		RangeMin = np.median(Bin_Velocity)-2*np.std(Bin_Velocity)		RangeMax = np.median(Bin_Velocity)+2*np.std(Bin_Velocity)		Good_Velocity = [Bin_Velocity[i] for i in range(len(Bin_Velocity)) if Bin_Velocity[i] >= RangeMin and Bin_Velocity[i] <= RangeMax]		return (np.median(Good_Velocity)-2*np.std(Good_Velocity), np.median(Good_Velocity)+2*np.std(Good_Velocity))	if type == "S":		data_dict = import_file(galaxy_index, type="S")		binned_data = data_dict[1] 		Bin_Dispersion = binned_data["SPXF"]		RangeMin = np.median(Bin_Dispersion)-2.5*np.std(Bin_Dispersion)		RangeMax = np.median(Bin_Dispersion)+2.5*np.std(Bin_Dispersion)		Good_Dispersion = [Bin_Dispersion[i] for i in range(len(Bin_Dispersion)) if Bin_Dispersion[i] >= RangeMin and Bin_Dispersion[i] <= RangeMax]		return (np.median(Good_Dispersion)-2.5*np.std(Good_Dispersion), np.median(Good_Dispersion)+2.5*np.std(Good_Dispersion))

	def velocity_map(galaxy_index):		data_dict = import_file(galaxy_index, type="V")	unbinned_data = data_dict[0]	binned_data = data_dict[1]	VMIN = get_MinMax_v2(galaxy_index,type="V")[0]	VMAX = get_MinMax_v2(galaxy_index,type="V")[1]	ax = display_bins_generators.display_bins_generators(binned_data["XBIN"], binned_data["YBIN"], 	                                     binned_data["VPXF"]-np.median(binned_data["VPXF"]), unbinned_data["X"], unbinned_data["Y"],	                                     vmin=VMIN, vmax=VMAX)
	#plt.contour(ax, colors='black')	
	return axdef dispersion_map(galaxy_index):		data_dict = import_file(galaxy_index, type="S")	unbinned_data = data_dict[0]	binned_data = data_dict[1]	SMIN = get_MinMax_v2(galaxy_index,type="S")[0]		#np.mean(binned_data ["SPXF"])-2*np.std(binned_data ["SPXF"])	SMAX = get_MinMax_v2(galaxy_index,type="S")[1]	ax = display_bins_generators.display_bins_generators(binned_data["XBIN"], binned_data["YBIN"], 	                                   binned_data["SPXF"], unbinned_data["X"], unbinned_data["Y"], 	                                   vmin=SMIN, vmax=SMAX)	#plt.plot(unbinned_data["X"], unbinned_data["Y"], '.')
	return ax
def photometry_image(figure, galaxy_index, X0_SAMI, Y0_SAMI, X0_FDS, Y0_FDS, Re):
	for i in np.arange(1,33):		
		field_file = '../0_data/Literature/FDS_Photometry_images/r'+str(int(i))+'_cropped.fits.fz'
		my_file = Path(field_file)		if my_file.exists():
			hdu_list = fits.open(field_file)
			image_data = hdu_list[1].data
			HD = hdu_list[1].header
			wcsHeader = WCS(hdu_list[1].header)
			field_RAmax, field_DECmin = wcsHeader.wcs_pix2world(0,0,0)		#RA is plotted inversely 55->54
			field_RAmin, field_DECmax = wcsHeader.wcs_pix2world(HD['NAXIS1']-1,HD['NAXIS2']-1,0)					
			#print(wcsHeader.wcs_world2pix(X0_FDS, Y0_FDS,0))			
			if X0_FDS >= field_RAmin and X0_FDS <= field_RAmax and Y0_FDS >= field_DECmin and Y0_FDS <= field_DECmax:		
				print("Found FCC"+str(int(galaxy_index))+" in Field NO."+str(int(i)))				
				#print("Field"+str(i), "RAmin:"+str(field_RAmin), "RAmax:"+str(field_RAmax), "DECmin:"+str(field_DECmin), "DECmax:"+str(field_DECmax))							
				px = figure.add_subplot(141, projection=wcsHeader)				
				FDS_pix = wcsHeader.wcs_world2pix(X0_FDS, Y0_FDS, 0)
				v = [FDS_pix[0]-Re*4/0.2, FDS_pix[0]+Re*4/0.2, FDS_pix[1]-Re*4/0.2, FDS_pix[1]+Re*4/0.2]
				px.axis(v,transform=px.transAxes)
				norm = ImageNormalize(image_data, ZScaleInterval(), stretch=LinearStretch())
				px.imshow(image_data, origin='lower', norm=norm, cmap='gray')			
				px.coords[0].set_ticklabel(size=18)
				px.coords[1].set_ticklabel(size=18) 
				px.set_title("FCC"+str(int(galaxy_index)),fontweight="bold", size=24)
				px.set_xlabel("RA (ddmmss)", fontsize = 20)
				px.set_ylabel("DEC (ddmmss)", fontsize = 20)				
				SAMI_pix = wcsHeader.wcs_world2pix(X0_SAMI, Y0_SAMI, 0)
				circle=plt.Circle((SAMI_pix[0], SAMI_pix[1]), 37.5, color='r', fill=False)
				px.add_artist(circle)


			else:
				continue
			hdu_list.close()

	return
	
def SetAxis(xx):
	#xx.set_title(str, fontsize=22) # Title
	xx.set_ylabel('arcsec', fontsize = 20) # Y label
	xx.set_xlabel('arcsec', fontsize = 20) # X label
	xx.tick_params(axis='both', which='major', labelsize=18)	return

def SetColBar(xx,str):
	cbar = plt.colorbar(xx, fraction = 0.0465)
	cbar.set_label(str, rotation=-270, fontsize=20)
	cbar.ax.tick_params(labelsize=18)
	return#************************MAIN_BODY******************************CENT = pd.read_csv('../2_pipeline/1_V&S_Maps/tmp/centers.csv')
FCC_index = CENT['FCC']

for n in range(len(FCC_index)):	
	Figure = plt.figure(figsize=(40,7))	
	
	X0_SAMI = CENT.SAMI_RA[CENT['FCC'] == FCC_index[n]]
	Y0_SAMI = CENT.SAMI_DEC[CENT['FCC'] == FCC_index[n]]	X0_FDS = CENT.FDS_RA[CENT['FCC'] == FCC_index[n]]
	Y0_FDS = CENT.FDS_DEC[CENT['FCC'] == FCC_index[n]]
	Re = CENT.Re[CENT['FCC'] == FCC_index[n]]
	print(n, FCC_index[n], X0_SAMI.values, Y0_SAMI.values, X0_FDS.values, Y0_FDS.values, Re.values)
	photometry_image(Figure, FCC_index[n], X0_SAMI.values, Y0_SAMI.values, X0_FDS.values, Y0_FDS.values, Re.values)
	
	vx = Figure.add_subplot(142)
	SetAxis(vx)
	vx = velocity_map(FCC_index[n])
	SetColBar(vx, '$V_{los}$ (km/s)')
	
	sx = Figure.add_subplot(143)
	SetAxis(sx)	
	sx = dispersion_map(FCC_index[n])	SetColBar(sx, '$\sigma$ (km/s)')
		
	M_r = CENT.M_r[CENT['FCC'] == FCC_index[n]]
	mu_r = CENT.mu_r[CENT['FCC'] == FCC_index[n]]
	Sersic = CENT.n[CENT['FCC'] == FCC_index[n]]
	string = "$M_r$ = " + str(M_r.values[0]) + "\n\n$R_e$ = " + str(Re.values[0]) + "\n\n$\mu_r$ = " + str(mu_r.values[0]) + "\n\n n = " + str(Sersic.values[0])
	plt.text(20, 0, string, fontsize=24 )
		Figure.savefig("../2_pipeline/1_V&S_Maps/" + str(int(FCC_index[n])) + "Velocity_map.pdf")