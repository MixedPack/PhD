from astropy.io import fitsfrom astropy.io import asciifrom plotbin import display_binsfrom plotbin import display_pixelsfrom plotbin import plot_velfieldfrom plotbin import display_bins_generatorsimport numpy as npimport matplotlib.pyplot as pltimport pandas as pdimport mathimport csvimport matplotlib.gridspec as gridspecfrom pathlib import Pathfrom astropy.wcs import WCSfrom astropy.visualization import (ZScaleInterval, ImageNormalize, LinearStretch)from astropy.nddata import Cutout2D
#************************FUNCTIONS******************************def photometry_image(galaxy_index, X0_FDS, Y0_FDS, Re):	if int(galaxy_index) == 83:
		field_file = '../0_data/Literature/FDS_Photometry_images/r19_cropped.fits.fz'	
		hdu_list = fits.open(field_file)		image_data = hdu_list[1].data		HD = hdu_list[1].header		hdu_list.close()
		wcsHeader = WCS(hdu_list[1].header)		field_RAmax, field_DECmin = wcsHeader.wcs_pix2world(0,0,0)		#RA is plotted inversely 55->54		field_RAmin, field_DECmax = wcsHeader.wcs_pix2world(HD['NAXIS1']-1,HD['NAXIS2']-1,0)							FDS_pix = wcsHeader.wcs_world2pix(X0_FDS, Y0_FDS, 0)		
		position = (FDS_pix[0],FDS_pix[1])    #(x,y)hdu
		size = (Re*8/0.2,Re*8/0.2) 			  # pixels (y,x)
		print(size)		
		cut_image=Cutout2D(image_data, position, size, wcs=wcsHeader, copy=True)					
		hdu=fits.PrimaryHDU(cut_image.data)
		hdu.header.update(cut_image.wcs.to_header())
		outfile = "../2_pipeline/1_V&S_Maps/tmp/" + str(int(galaxy_index)) + "cut_image.fits"
		hdu.writeto(outfile, overwrite=True)	
	
	else:	
		for i in np.arange(1,33):						
			field_file = '../0_data/Literature/FDS_Photometry_images/r'+str(int(i))+'_cropped.fits.fz'			my_file = Path(field_file)			if my_file.exists():				hdu_list = fits.open(field_file)				image_data = hdu_list[1].data				HD = hdu_list[1].header				hdu_list.close()
				wcsHeader = WCS(hdu_list[1].header)				field_RAmax, field_DECmin = wcsHeader.wcs_pix2world(0,0,0)		#RA is plotted inversely 55->54				field_RAmin, field_DECmax = wcsHeader.wcs_pix2world(HD['NAXIS1']-1,HD['NAXIS2']-1,0)													if X0_FDS >= field_RAmin and X0_FDS <= field_RAmax and Y0_FDS >= field_DECmin and Y0_FDS <= field_DECmax:							print("Found FCC"+str(int(galaxy_index))+" in Field NO."+str(int(i)))																	FDS_pix = wcsHeader.wcs_world2pix(X0_FDS, Y0_FDS, 0)		
					position = (FDS_pix[0],FDS_pix[1])    #(x,y)hdu
					size = (Re*8/0.2,Re*8/0.2) 			  # pixels (y,x)
					cut_image=Cutout2D(image_data, position, size, wcs=wcsHeader, copy=True)				
					#TEST1				
					#norm = ImageNormalize(image_data, ZScaleInterval(), stretch=LinearStretch())								
					#plt.imshow(image_data, origin='lower', norm=norm, cmap='gray')
					#cut_image.plot_on_original(color='white')				
					#TEST2				
					#plt.imshow(cut_image.data, origin='lower', norm=norm, cmap='gray')
					#plt.show()
					hdu=fits.PrimaryHDU(cut_image.data)
					hdu.header.update(cut_image.wcs.to_header())
					outfile = "../2_pipeline/1_V&S_Maps/tmp/" + str(int(galaxy_index)) + "cut_image.fits"
					hdu.writeto(outfile, overwrite=True)
				else:					continue	return
#************************MAIN_BODY******************************CENT = pd.read_csv('../2_pipeline/1_V&S_Maps/tmp/centers.csv')FCC_index = CENT['FCC']for n in range(len(FCC_index)):		X0_FDS = CENT.FDS_RA[CENT['FCC'] == FCC_index[n]]	Y0_FDS = CENT.FDS_DEC[CENT['FCC'] == FCC_index[n]]	Re = CENT.Re[CENT['FCC'] == FCC_index[n]]	print(n, FCC_index[n], X0_FDS.values, Y0_FDS.values, Re.values)	photometry_image(FCC_index[n], X0_FDS.values, Y0_FDS.values, Re.values)	