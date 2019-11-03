from scipy.optimize import minimize
Sfit_params, Spcov = optm.curve_fit(parabola, Sx1, Sy)
plt.plot(X_fit, Y_fit, 'k-')
#plt.plot(X2_fit, Y2_fit, 'k-.')
ax.set_ylabel("log($(M_{dyn})_{1/2}$ / $L_{1/2}$)  [$M_\odot/L_\odot$]", fontsize = 20)
ymin = -Sfit_params[1]**2/(4*Sfit_params[0])+Sfit_params[2]
ax.annotate("[%.3f"%xmin+",%.3f"%ymin+"]", xy=(xmin, ymin), xytext=(xmin-0.25, ymin+0.9),
            arrowprops=dict(arrowstyle="->"),fontsize = 16)
#xmin2 = -Sfit_params2[1]/(2*Sfit_params2[0])
#ymin2 = -Sfit_params2[1]**2/(4*Sfit_params2[0])+Sfit_params2[2]
#ax.annotate("[%.3f"%xmin2+",%.3f"%ymin2+"] WFT", xy=(xmin2, ymin2), xytext=(xmin2, ymin2+1.2),
#            arrowprops=dict(arrowstyle="->"),fontsize = 16)

plt.plot(X_fit, Y_fit, 'k-')
#plt.plot(X2_fit, Y2_fit, 'k-.')
bx.set_xlabel("log($(M_{dyn})_{1/2}$) [$M_\odot$]", fontsize = 20)
ymin = -DMfit_params[1]**2/(4*DMfit_params[0])+DMfit_params[2]
bx.annotate("[%.3f"%xmin+",%.3f"%ymin+"]", xy=(xmin, ymin), xytext=(xmin-0.75, ymin+0.9),
            arrowprops=dict(arrowstyle="->"),fontsize = 16)
#xmin2 = -DMfit_params2[1]/(2*DMfit_params2[0])
#ymin2 = -DMfit_params2[1]**2/(4*DMfit_params2[0])+DMfit_params2[2]
#bx.annotate("[%.3f"%xmin2+",%.3f"%ymin2+"] WFT", xy=(xmin2, ymin2), xytext=(xmin2, ymin2+1.2),
#            arrowprops=dict(arrowstyle="->"),fontsize = 16)
PB_WSTF, = plt.plot(X_fit, Y_fit, 'k-')
#PB_WTF, = plt.plot(X2_fit, Y2_fit, 'k-.')
cx.set_xlabel("log($L_{1/2}$) [$L_\odot$]", fontsize = 20)
ymin = -Lfit_params[1]**2/(4*Lfit_params[0])+Lfit_params[2]
cx.annotate("[%.3f"%xmin+",%.3f"%ymin+"]", xy=(xmin, ymin), xytext=(xmin-1.75, ymin+0.9),
            arrowprops=dict(arrowstyle="->"),fontsize = 16)
#xmin2 = -Lfit_params2[1]/(2*Lfit_params2[0])
#ymin2 = -Lfit_params2[1]**2/(4*Lfit_params2[0])+Lfit_params2[2]
#cx.annotate("[%.3f"%xmin2+",%.3f"%ymin2+"] WFT", xy=(xmin2, ymin2), xytext=(xmin2-1.0, ymin2+1.2),
#            arrowprops=dict(arrowstyle="->"),fontsize = 16)

plt.legend(fontsize=16)