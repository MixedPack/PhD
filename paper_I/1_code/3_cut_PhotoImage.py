from astropy.io import fits

		field_file = '../0_data/Literature/FDS_Photometry_images/r19_cropped.fits.fz'	
		hdu_list = fits.open(field_file)
		wcsHeader = WCS(hdu_list[1].header)
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
			field_file = '../0_data/Literature/FDS_Photometry_images/r'+str(int(i))+'_cropped.fits.fz'
				wcsHeader = WCS(hdu_list[1].header)
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

