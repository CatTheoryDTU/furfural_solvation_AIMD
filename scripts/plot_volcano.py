###1D volcano 1bar H2, 0.0001bar FCHO (liquid), 300K
import numpy as np
from scipy.interpolate import make_interp_spline,interp1d

fig,ax=plt.subplots(figsize=(9,6))
##mkm results
surfaces = ['Ru','Rh','Pt','Ni','Pd','Cu','Au','Ag']
raw_data = [(-2.46,4.581483345280234e-15),(-2.18,3.782165749918155e-10),(-1.91,1.2685487599132238e-05),(-1.87,0.0002294043364320653),(-1.67,0.0015799689400214395),(-0.93,0.4524887316430913),(-0.84,0.0002216166681327576),(-0.79,6.929966740731143e-07)]
X = [x[0] for x in raw_data]
Y = [np.log10(y[1]) for y in raw_data]
fit = np.polyfit(X, Y, 2)
f = np.poly1d(fit)

# surfaces_ = ['Ag(111)','Co(0001)','Ni(111)','Pd(111)','Ru(0001)']
# Ew_ = [-0.28,-0.43,-0.46,-0.41,-0.53]

xx = np.linspace(-5,5,1000)
ax.plot(xx, f(xx), color='black')

##metal
ax.scatter(X,Y,color='black',s=100)

##with_solvation
surfaces_solv = ['Au','Cu','Pd','Pt','Rh']
Eaimd = [[-0.57,-0.46,-0.10],[-0.11,-0.32,-0.40],[-1.13,-1.19,-0.77],[-1.05,-1.00,-0.64],[-1.19,-1.08, -1.01]]
Evac = [-0.84,-0.93,-1.67,-1.91,-2.08]
Ew = [-0.30,-0.39,-0.42,-0.41,-0.45]

Eaimd_avg = [sum(i)/len(i) for i in Eaimd]
Es_v = np.subtract(Eaimd_avg, Evac)
print(Eaimd_avg)
x_error = [np.std(i) for i in Eaimd]
print(x_error)

colors = ['gold','brown','royalblue','silver','lightblue']
for i, surf in enumerate(surfaces_solv):
    ax.errorbar(Eaimd_avg[i], f(Eaimd_avg[i]), xerr = x_error[i], fmt='o',capsize=5,color='c',markersize=12)

for i,s in enumerate(surfaces):
    ax.annotate(s,(X[i]+0.02,Y[i]+0.03),fontsize=20)
    
for i,s in enumerate(surfaces_solv):
    ax.annotate(s,(Eaimd_avg[i]+0.02,f(Eaimd_avg[i])+0.08),fontsize=20,color='c')
    

ax.scatter(10,10,label='Vacuum',color='black', s = 80)
ax.scatter(10,10,label='Aqueous',color='c', s = 80)
ax.set_ylim(-16,2)
ax.set_xlim(-2.6,-0)
ax.tick_params(labelsize=25)
ax.set_xlabel("$\Delta E_{FCHO}$ (eV)", fontsize = 30)
ax.set_ylabel("log(TOF (s$^{-1}$))", fontsize = 30)
ax.legend(fontsize = 20, loc='lower center',frameon =False)

fig.tight_layout()
