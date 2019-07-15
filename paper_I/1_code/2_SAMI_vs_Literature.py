from scipy.optimize import minimize
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import polyfit
import math
import pandas as pd

#*****************************Input values**********************************************
Galaxies = pd.read_csv('../2_pipeline/0_Galaxies_Table/Galaxies_Table_DWARFS.csv')	
FCC_index = Galaxies['FCC']
Kinematics = pd.read_csv('../2_pipeline/1_Kinematics_Table/Kinematics_Table_DWARFS.csv')
Dispersion = Kinematics['Sigma']
ERR_Disp = Kinematics['Error.Sigma']

#*****************************Calculation***********************************************
#*****************************Toloba et al. 2014 (2 papers)**********************************
#*****************************Plotting**************************************************
Figure = plt.figure(figsize=(20,7))

ax = Figure.add_subplot(131)

kwargs = dict(histtype='stepfilled', alpha=0.3, density=True, bins=[-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13])
plt.hist(M_V, **kwargs)
plt.hist(F11_M_V, **kwargs)
plt.hist(T14_M_V, **kwargs)

ax = Figure.add_subplot(132)
kwargs = dict(histtype='stepfilled', alpha=0.3, density=True, bins=[7,8,9,10,11,12])
plt.hist(log_stel_mass_V, **kwargs)
plt.hist(F11_log_stMass_V , **kwargs)
plt.hist(T14_log_StMass_V, **kwargs)

Figure.savefig("../2_pipeline/2_SAMI_vs_Literature/SAMIvsLiter.pdf")
