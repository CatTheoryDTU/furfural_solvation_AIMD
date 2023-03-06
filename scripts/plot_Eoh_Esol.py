###Esol-Eoh

import matplotlib.pyplot as  plt
import numpy as np
import pandas as pd

fig,ax = plt.subplots(figsize =(8, 6))
ax.grid(True)
###raw data

surfaces = ['Au','Cu','Pd','Pt','Rh']
Eaimd = [[-0.57,-0.46],[-0.11,-0.32],[-1.19,-0.77],[-1.15,-1.00],[-1.19,-1.08,-1.01]]
Evac = [-0.84,-0.93,-1.67,-1.91,-2.08]
Ew = [-0.30,-0.39,-0.42,-0.41,-0.45]
Eoh = [1.21,0.56,0.74,0.68,0.35]


Eaimd_avg = [sum(i)/len(i) for i in Eaimd]
Es_v = np.subtract(Eaimd_avg, Evac)
print(Eaimd_avg)
y_error = [np.std(i) for i in Eaimd]
print(y_error)

colors = ['gold','brown','royalblue','silver','lightblue']
for i, surf in enumerate(surfaces):
    ax.errorbar(Eoh[i], Es_v[i], yerr = y_error[i], fmt='o',capsize=10,color=colors[i],markersize=15)

a1,b1 = np.polyfit(Eoh,Es_v,1)
variance = np.var(Es_v)
residuals = np.var([(a1*xx + b1 - yy)  for xx,yy in zip(Eoh,Es_v)])
Rsqr = np.round(1-residuals/variance, decimals=2)

xx = np.linspace(-5,5,100)
ax.plot(xx, a1*xx+b1, color='black',ls='--',lw=3,label = str(round(a1,2))+'x + '+str(round(b1,2))+' = y, R$^2$ = '+str(round(Rsqr,2)))

ax.set_ylabel('$\Delta E_{sol}$ (eV)', fontsize=30)
ax.tick_params(labelsize=30)
ax.set_xlim(0,1.5)
ax.set_ylim(0,1.5)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.legend(fontsize=20,loc=1, frameon = False)
fig.tight_layout()
