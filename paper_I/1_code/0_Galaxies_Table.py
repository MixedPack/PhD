import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import csv
import pyfits
import math
from numpy.polynomial.polynomial import polyfit

#The list of Good SAMI_2015&2016 galaxies - with their FCC and FDS names
#with open( '../0_data/Literature/FCC_FDS.csv', 'r' ) as f:
with open( '../0_data/tmp/NotChosen_FCC_FDS.csv', 'r' ) as f:
    name = csv.reader(f)
    name_list = list(map(tuple, name))

size = len(name_list)

#Online Vizier Catalog - Ferguson
hdulist_Ferguson = pyfits.open('../0_data/Literature/FergusonCatalog.fit')

#FDS's Catalog (2018) (faint and giant Fornax galaxies with their errors)
hdulist_FDS = pyfits.open('../0_data/Literature/FDSCatalog_AllError.fits')
tbdata_FDS = hdulist_FDS[1].data

COV_mr_logRe2 = -0.00528375

#importing information from FDS's catalog
M_r = np.zeros(size)
M_g = np.zeros(size)
R_eff = np.zeros(size)
mu_r = np.zeros(size)
axis_ratio = np.zeros(size)
RA = np.zeros(size)
DEC = np.zeros(size)
log_mass = np.zeros(size)
ERR_m_r = np.zeros(size)
ERR_m_g = np.zeros(size)
ERR_mu_r = np.zeros(size)
ERR_R_eff = np.zeros(size)
ERR_log_mass = np.zeros(size)
sersic_index = np.zeros(size)

#stellar_mass = np.zeros(size)
u = np.zeros(size)
g = np.zeros(size)
r = np.zeros(size)
i = np.zeros(size)

for l in range(size):

	draft = name_list[l]
	name_FCC = 'FCC'+str(draft[0])
	name_FDS = 'FDS'+str(draft[1])+'_DWARF'+str("{0:0=3d}".format(int(draft[2])))
	
	for j in range(len(tbdata_FDS.field('Target'))):
		if tbdata_FDS.field('Target')[j]==name_FDS:
			RA[l] = tbdata_FDS.field('RA')[j]
			DEC[l] = tbdata_FDS.field('DEC')[j]
			M_r[l] = tbdata_FDS.field('r_mag')[j] - 31.0
			ERR_m_r[l] = tbdata_FDS.field('r_mag_e')[j]
			M_g[l] = tbdata_FDS.field('g_mag')[j] - 31.0
			ERR_m_g[l] = tbdata_FDS.field('g_mag_e')[j]
			R_eff[l] = tbdata_FDS.field('reff')[j]
			ERR_R_eff[l] = tbdata_FDS.field('reff_e')[j]
			axis_ratio[l] = tbdata_FDS.field('arat')[j]
			mu_r[l] = tbdata_FDS.field('r_mag')[j] + 2.5 * np.log10(math.pi*axis_ratio[l]*math.pow(R_eff[l],2)) + 2.5 * np.log10(2)
			ERR_mu_r[l] = ERR_m_r[l] + (5 / np.log(10)) * (ERR_R_eff[l] / R_eff[l]) + 2.5 * COV_mr_logRe2
			u[l] = tbdata_FDS.field('u')[j]
			g[l] = tbdata_FDS.field('g')[j]
			r[l] = tbdata_FDS.field('r')[j]
			i[l] = tbdata_FDS.field('i')[j]
			ERR_g = tbdata_FDS.field('g_e')[j]
			ERR_r = tbdata_FDS.field('r_e')[j]
			ERR_i = tbdata_FDS.field('i_e')[j]
			log_mass[l] = 1.15 + 0.70*(g[l]-i[l]) - 0.4*M_r[l] + 0.4*(r[l]-i[l])
			ERR_log_mass[l] = np.sqrt(0.49*(ERR_g**2+ERR_i**2) + 0.16*ERR_m_r[l]**2 + 0.16*(ERR_r**2+ERR_i**2))/(log_mass[l]*np.log(10))
			sersic_index[l] = tbdata_FDS.field('n')[j]

xFit = [g[k]-r[k] for k in range(size) if u[k]-g[k] > 0.0]
yFit = [u[k]-g[k] for k in range(size) if u[k]-g[k] > 0.0]
b, m = polyfit(xFit, yFit, 1)
Fitline = [b + m * xFit[k] for k in range(len(xFit))]
plt.plot(xFit, yFit, '.')
plt.plot(xFit, Fitline, c='grey')

for l in range(size):
	draft = name_list[l]
	name_FCC = 'FCC' + str(draft[0])
	if name_FCC=='FCC37' or name_FCC=='FCC46' or name_FCC=='FCC33' or name_FCC=='FCC29':
			u[l] = g[l] + b + (g[l]-r[l])*m
			plt.plot(g[l]-r[l], u[l]-g[l], 'x')
#plt.show()

#writing a csv output
ListDict = []
for l in range(size):
	draft = name_list[l]
	dict={'FCC':draft[0], 'FDS':draft[1], '_DWARF':"{0:0=3d}".format(int(draft[2])), 'RA(deg)':RA[l], 'DEC(deg)':DEC[l],
			 'M_r(mag)':"%.2f" % M_r[l], 'ERR_m_r(mag)':"%.2f" % ERR_m_r[l], 'M_g(mag)':"%.2f" % M_g[l], 'ERR_m_g(mag)':"%.2f" % ERR_m_g[l] ,
			'mu_r(mag/arcsec2)':"%.2f" % mu_r[l], 'ERR_mu_r(mag/arcsec2)':"%.2f" % ERR_mu_r[l] , 
			'R_e(arcsec)':"%.2f" % R_eff[l], 'ERR_R_e(arcsec)':"%.2f" % ERR_R_eff[l] , 
			'axis_ratio':"%.2f" % axis_ratio[l], 'Sersic_index':"%.2f" % sersic_index[l],'log10(M_*/M_sun)':"%.4f" % log_mass[l], 'ERR_log10(M_*)':"%.4f" % ERR_log_mass[l],
			'u':"%.4f" % u[l], 'g':"%.4f" % g[l], 'r':"%.4f" % r[l], 'i':"%.4f" % i[l]}
	ListDict.append(dict.copy())

#with open('../2_pipeline/0_Galaxies_Table/Galaxies_Table.csv', 'w') as output_file:
with open('../2_pipeline/0_Galaxies_Table/tmp/NotChosen_Galaxies_Table.csv', 'w') as output_file:
    dict_writer = csv.DictWriter(output_file, dict.keys())
    dict_writer.writeheader()
    dict_writer.writerows(ListDict)

#******************************ERRORS - Seperate File*************************************