###compare binding energies of furfural: vacuum, implicit solvation, aimd, bond-additivity, exp.

import matplotlib.pyplot as  plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FormatStrFormatter

fig = plt.figure(figsize=(10,12))
gs = fig.add_gridspec(2, hspace=0.1,height_ratios=[1.5,1])
ax, ax1 = gs.subplots(sharex=True)

###raw data
surfaces = ['Au','Cu','Pd','Pt','Rh']
Eaimd = [[-0.57,-0.46,-0.10],[-0.11,-0.32,-0.40],[-1.13,-1.19,-0.77],[-1.05,-1.00,-0.64],[-1.19,-1.08, -1.01]]
Evac = [-0.84,-0.93,-1.67,-1.91,-2.08]
Eimp = [-1.08,-1.15,-1.93,-2.23,-2.29] 
Ehex = [-0.06,0.21,-0.67,-1.08,-1.16] 
Eexp = [None, None,None, -0.89,-0.88]
Eba = [None, None,None,-0.96,-0.69]

Eaimd_avg = [sum(i)/len(i) for i in Eaimd]
print(Eaimd_avg)
y_error = [np.std(i) for i in Eaimd]
print(y_error)

colors = ['gold','brown','silver','lightblue']
xx = range(5)

Esol_aimd = [Eaimd_avg[i]-Evac[i] for i in xx]
print(Esol_aimd)
Esol_imp = [Eimp[i]-Evac[i] for i in xx]
Esol_hex = [Ehex[i]-Evac[i] for i in xx]


##adsorption energy
ax.scatter(xx, Evac, s=200, marker = 'x',color='blue',label='Vacuum, this work')
ax.scatter(xx, Eimp, s=200, marker = 'x',color='red',label='Implicit solvation, this work')

ax.errorbar(xx, Eaimd_avg, yerr = y_error, fmt='o',capsize=5,color='black',markersize=10,label='AIMD, this work')
ax.scatter(xx, Eba, s=200, marker = '^',color='brown',label='Bond-additivity model')
ax.scatter(xx, Eexp, s=200, marker = '*',color='green',label='Experiment')

##solvation energy
ax1.errorbar(xx, Esol_aimd, yerr = y_error, fmt='o',capsize=5,color='black',markersize=10)

xxx = np.linspace(-1,10,100)
yyy = xxx*0
ax1.plot(xxx,yyy,'--',color='black')

#plt
ax.set_ylabel('$\Delta E$ (eV)', fontsize=25)
ax1.set_ylabel('$\Delta E_{solv}$ (eV)', fontsize=25)
ax.set_xticks(xx)
ax.set_xticklabels(surfaces)

ax.tick_params(labelsize=20)
ax1.tick_params(labelsize=20)
ax.set_xlim(-0.5, 4.5)
ax.set_ylim(-2.5, 0.5)
ax1.set_ylim(-0.5, 1.5)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.legend(fontsize=14,frameon=False,loc=1,ncol=1)
plt.show()
plt.tight_layout()
