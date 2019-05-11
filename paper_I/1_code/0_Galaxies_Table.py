import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import csv
import pyfits
import math

#The list of Good SAMI_2015&2016 galaxies - with their FCC and FDS names
with open( '../0_data/FCC_FDS_1516.csv', 'r' ) as f:
    name = csv.reader(f)
    name_list = list(map(tuple, name))
with open('../0_data/FCC_FDS_18.csv', 'r') as f:
    name2 = csv.reader(f)
    name_list2 = list(map(tuple, name2))

name_list.extend(name_list2)
size = len(name_list)

#Online Vizier Catalog - Ferguson
hdulist_Ferguson = pyfits.open('./FergusonCatalog.fit')

#Aku's Catalog (2018) - Fornax Cluster Galaxies
hdulist_Aku = pyfits.open('./fds_cluster.fits')
tbdata_Aku = hdulist_Aku[1].data


#importing information from Aku's catalog
abs_mag_r = np.zeros(size)
R_eff = np.zeros(size)
surface_bright = np.zeros(size)
axis_ratio = np.zeros(size)
RA = np.zeros(size)
DEC = np.zeros(size)
log_mass = np.zeros(size)
#stellar_mass = np.zeros(size)
u = np.zeros(size)
g = np.zeros(size)
r = np.zeros(size)
i = np.zeros(size)

for l in range(size):
	draft = name_list[l]
	name_FCC = 'FCC'+str(draft[0])
	name_FDS = 'FDS'+str(draft[1])+'_DWARF'+str("{0:0=3d}".format(int(draft[2])))
	
	for j in range(len(tbdata_Aku.field('target'))):
		if tbdata_Aku.field('target')[j]==name_FDS:
			RA[l] = tbdata_Aku.field('RA')[j]
			DEC[l] = tbdata_Aku.field('DEC')[j]
			abs_mag_r[l] = tbdata_Aku.field('r_mag')[j] - 31.0
			R_eff[l] = tbdata_Aku.field('Reff')[j]
			axis_ratio[l] = tbdata_Aku.field('axis_ratio')[j]
			surface_bright[l] = tbdata_Aku.field('r_mag')[j] + 2.5 * np.log10(math.pi*axis_ratio[l]*math.pow(R_eff[l],2)) + 2.5 * np.log10(2)
			u[l] = tbdata_Aku.field('u')[j]
			g[l] = tbdata_Aku.field('g')[j]
			r[l] = tbdata_Aku.field('r')[j]
			i[l] = tbdata_Aku.field('i')[j]
			log_mass[l] = 1.15 + 0.70*(g[l]-i[l]) - 0.4*abs_mag_r[l] + 0.4*(r[l]-i[l])
			#stellar_mass[l] = math.pow(10,log_mass[l]) # in solar mass unit

#writing a csv output
ListDict = []
for l in range(size):
	draft = name_list[l]
	dict={'FCC':draft[0], 'FDS':draft[1], '_DWARF':"{0:0=3d}".format(int(draft[2])), 'RA(deg)':RA[l], 
			'DEC(deg)':DEC[l], 'M_r(mag)':"%.2f" % abs_mag_r[l], 'mu_r(mag/arcsec2)':"%.2f" % surface_bright[l], 
			'R_e(arcsec)':"%.2f" % R_eff[l], 'axis_ratio':"%.2f" % axis_ratio[l], 
			'log10(M_*/M_sun)':"%.4f" % log_mass[l], 'u':"%.4f" % u[l], 'g':"%.4f" % g[l], 
			'r':"%.4f" % r[l], 'i':"%.4f" % i[l]}
	ListDict.append(dict.copy())

with open('Galaxies_Table.csv', 'w') as output_file:
    dict_writer = csv.DictWriter(output_file, dict.keys())
    dict_writer.writeheader()
    dict_writer.writerows(ListDict)

