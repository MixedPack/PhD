from scipy.optimize import minimize
g_r = Galaxies['g'] - Galaxies['r']

F11_ML_stel_V = Falcon11['M\L*_V(M/Lsun)']
Falcon11_2 = pd.read_csv('../0_data/Literature/FalconBarroso11_Phot.csv')
F11_log_stMass_V = np.log10(F11_ML_stel_V*(2*F11_Lumin_half_V))
ERR_F11_log_stMass_V = np.sqrt(np.power(ERR_F11_ML_stel_V/(np.log(10)*F11_ML_stel_V),2)+np.power(ERR_F11_Lum_V/(np.log(10)*F11_Lumin_half_V),2))

#writing a csv output
index = Falcon11['Galaxy']
T14_log_StMass_V = T14_log_StMass_r - 0.4*0.16
T14_ML_V = T14_ML_r * (10**0.16)
Toloba14_2 = pd.read_csv('../0_data/Literature/Toloba14_PhotKinem.csv')
T14_L_half_V = 0.5 * np.power(10, 0.4*(M_V_sun - T14_M_V))
ERR_T14_log_DMass = Toloba14['log(M_e,dyn)_ERR']
ERR_T14_log_StMass_V = ERR_T14_log_StMass_r
ERR_T14_L_half_V = 0 #for the time being

#*****************************Plotting**************************************************
#DynMass vs. Lum_1/2
plt.plot(np.log10(Dyn_Mass_half),np.log10(Lumin_half_V),'r.')
plt.plot(np.log10(W_DMass),np.log10(W_Lumin),'bx')
plt.plot(np.log10(F11_DMass_half),np.log10(F11_Lumin_half_V),'g+')
plt.plot(T14_log_DMass, np.log10(T14_L_half_V),'y*')
ax.set_xlabel("log($(M_{dyn})_{1/2}$) [$M_\odot$]", fontsize = 20)
ax.tick_params(axis='both', which='major', labelsize=18)

#DynMass vs. M/L
plt.plot(np.log10(Dyn_Mass_half),np.log10(ML_V),'r.')
plt.plot(np.log10(W_DMass),np.log10(W_ML),'bx')
plt.plot(np.log10(F11_DMass_half),np.log10(F11_ML_V),'g+')
plt.plot(T14_log_DMass, np.log10(T14_ML_V),'y*')
bx.set_xlabel("log($(M_{dyn})_{1/2}$) [$M_\odot$]", fontsize = 20)
bx.tick_params(axis='both', which='major', labelsize=18)

#Lum_1/2 vs. M/L
plt.plot(np.log10(Lumin_half_V),np.log10(ML_V),'r.')
plt.plot(np.log10(W_Lumin),np.log10(W_ML),'bx')
plt.plot(np.log10(F11_Lumin_half_V),np.log10(F11_ML_V),'g+')
plt.plot(np.log10(T14_L_half_V), np.log10(T14_ML_V),'y*')
cx.set_xlabel("log($L_{1/2}$) [$L_\odot$]", fontsize = 20)

plt.plot(log_stel_mass_V, np.log10(Dyn_Mass_half),'r.')
plt.plot(T14_log_StMass_V, T14_log_DMass,'y*')
plt.plot(F11_log_stMass_V, np.log10(F11_DMass_half),'g+')
ax.set_xlabel("log($M_{stel}$) [$M_\odot$]")
#M_* vs. M_dyn/M_*
plt.plot(log_stel_mass_V, np.log10(Dyn_Mass_half)-log_stel_mass_V,'r.')
plt.plot(T14_log_StMass_V, T14_log_DMass-T14_log_StMass_V,'y*')
#Falcon-Barroso et al. 2011 (SAURON))
plt.plot(F11_log_stMass_V, np.log10(F11_DMass_half)-F11_log_stMass_V,'g+')
bx.set_xlabel("log($M_{stel}$) [$M_\odot$]")