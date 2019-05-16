import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import csv
import pyfits
import math

#The list of Good SAMI_2015&2016 galaxies - with their FCC and FDS names
with open( '../0_data/Literature/FCC_FDS.csv', 'r' ) as f:
    name = csv.reader(f)
    name_list = list(map(tuple, name))

size = len(name_list)

#Online Vizier Catalog - Ferguson
hdulist_Ferguson = pyfits.open('../0_data/Literature/FergusonCatalog.fit')

#FDS's Catalog (2018) (faint and giant Fornax galaxies)
hdulist_FDS = pyfits.open('../0_data/Literature/FDSCatalog_ALL.fits')
tbdata_FDS = hdulist_FDS[1].data
#FDS catalog with errors (only dwarfs)
hdulist_FDSerr = pyfits.open('../0_data/Literature/FDSCatalog_DwarfErr.fit')
tbdata_FDSerr = hdulist_FDSerr[1].data

COV_mr_logRe2 = -0.00528375

#importing information from FDS's catalog
M_r = np.zeros(size)
R_eff = np.zeros(size)
mu_r = np.zeros(size)
axis_ratio = np.zeros(size)
RA = np.zeros(size)
DEC = np.zeros(size)
log_mass = np.zeros(size)
ERR_m_r = np.zeros(size)
ERR_mu_r = np.zeros(size)
ERR_R_eff = np.zeros(size)

#stellar_mass = np.zeros(size)
u = np.zeros(size)
g = np.zeros(size)
r = np.zeros(size)
i = np.zeros(size)

for l in range(size):
	draft = name_list[l]
	name_FCC = 'FCC'+str(draft[0])
	name_FDS = 'FDS'+str(draft[1])+'_DWARF'+str("{0:0=3d}".format(int(draft[2])))
	
	for j in range(len(tbdata_FDS.field('target'))):
		if tbdata_FDS.field('target')[j]==name_FDS:
			RA[l] = tbdata_FDS.field('RA')[j]
			DEC[l] = tbdata_FDS.field('DEC')[j]
			M_r[l] = tbdata_FDS.field('r_mag')[j] - 31.0
			R_eff[l] = tbdata_FDS.field('Reff')[j]
			axis_ratio[l] = tbdata_FDS.field('axis_ratio')[j]
			mu_r[l] = tbdata_FDS.field('r_mag')[j] + 2.5 * np.log10(math.pi*axis_ratio[l]*math.pow(R_eff[l],2)) + 2.5 * np.log10(2)
			u[l] = tbdata_FDS.field('u')[j]
			g[l] = tbdata_FDS.field('g')[j]
			r[l] = tbdata_FDS.field('r')[j]
			i[l] = tbdata_FDS.field('i')[j]
			log_mass[l] = 1.15 + 0.70*(g[l]-i[l]) - 0.4*M_r[l] + 0.4*(r[l]-i[l])
			#stellar_mass[l] = math.pow(10,log_mass[l]) # in solar mass unit
	for k in range(len(tbdata_FDSerr.field('Target'))):	
		if tbdata_FDSerr.field('Target')[k]==name_FDS:
			ERR_m_r[l] = tbdata_FDSerr.field('e_rmag')[k]
			ERR_R_eff[l] = tbdata_FDSerr.field('e_reff')[k]
			R_eff[l] = tbdata_FDSerr.field('reff')[k]
			ERR_mu_r[l] = ERR_m_r[l] + (5/np.log(10))*(ERR_R_eff[l]/R_eff[l])+2.5*COV_mr_logRe2
			
#writing a csv output
ListDict = []
for l in range(size):
	draft = name_list[l]
	dict={'FCC':draft[0], 'FDS':draft[1], '_DWARF':"{0:0=3d}".format(int(draft[2])), 'RA(deg)':RA[l], 'DEC(deg)':DEC[l],
			 'M_r(mag)':"%.2f" % M_r[l], 'ERR_m_r(mag)':"%.2f" % ERR_m_r[l] ,
			'mu_r(mag/arcsec2)':"%.2f" % mu_r[l], 'ERR_mu_r(mag/arcsec2)':"%.2f" % ERR_mu_r[l] , 
			'R_e(arcsec)':"%.2f" % R_eff[l], 'ERR_R_e(arcsec)':"%.2f" % ERR_R_eff[l] , 
			'axis_ratio':"%.2f" % axis_ratio[l], 'log10(M_*/M_sun)':"%.4f" % log_mass[l],  
			'u':"%.4f" % u[l], 'g':"%.4f" % g[l], 'r':"%.4f" % r[l], 'i':"%.4f" % i[l]}
	ListDict.append(dict.copy())

with open('../2_pipeline/0_Galaxies_Table/Galaxies_Table.csv', 'w') as output_file:
    dict_writer = csv.DictWriter(output_file, dict.keys())
    dict_writer.writeheader()
    dict_writer.writerows(ListDict)

#******************************ERRORS - Seperate File*************************************