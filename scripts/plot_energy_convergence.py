"""Report the dipole moment with the right errors."""
import json
import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint   
import collections

def check_convergence(data, stdev=0.02):
    """Check if the dipole data is converged.
    If the last 1000 fs have a stdev less than std, then 
    I will consider it converged. If it doesn't have at least 
    2000 ps of runtime the trajectory will by default 
    be considered not converged."""
    if len(data) < 10000:
        # Run at least 2000 ps for me to consider this trajectory
        return False
    last_ps = data[-1000:]
    if np.std(last_ps) < stdev:
        return True
    else:
        return False

metals = ['Au_111','Cu_111','Pd_111','Pt_111','Rh_111']
metal_labels = ['Au(111)','Cu(111)','Pd(111)','Pt(111)','Rh(111)']
ads = ['clean','fur']
runs = [1,2,3]
lss = ['-','--']

if __name__ == '__main__':
    colors = ['silver','grey','black']
    fig = plt.figure(figsize=(14,6))
    gs = fig.add_gridspec(ncols=5, nrows=1, wspace=0)
    axs = gs.subplots(sharex=True, sharey=True)
    #plt.ylabel(r'Energy (eV)',fontsize=20) 
    #plt.xlabel(r'Time (ps)',fontsize=20)
    CUTOFF = 1000 # Remove these many enteries from the data
    for im, m in enumerate(metals):
        axs[im].plot([], [], label=metal_labels[im], c = 'black') #add legends
        for ia, a in enumerate(ads):
            for i, run in enumerate(runs):
            # This json file has all the timestep and dipole data
                with open(m+'/'+a+'-40w-'+str(run)+'/1e-5/'+'output/test.json', 'r') as handle:
                    data = json.load(handle)
    # Store the converged dipole moments here
                data_conv = {}
                for structure, results in data.items():
                    print(f'Parsing data for {structure}')
                    print(structure)
                    energies = np.array(results)[1]

        # Check if the dipole moments are converged
                    if check_convergence(energies, stdev=0.005):
                        converged_energies = energies[-1]
                        print(f'Converged energies: {converged_energies}')
                        data_conv[structure] = converged_energies 
                        results = np.array(results) 
                    # remove the first CUTOFF entries

                        energies = results[1][CUTOFF:]
                        timestep = results[0][CUTOFF:]
                        axs[im].plot(timestep, energies-energies[0], lw = 2, color=colors[i], ls = lss[ia])
                #with open('output/converged_energies_'+str(run)+'.json', 'w') as handle:
                #json.dump(data_conv, handle, indent=4)

            # print the sorted energies
                #sorted_energies = collections.OrderedDict(sorted(data_conv.items()))
            # for a in ax:
            axs[im].set_xlim([min(timestep), 51])
        #ax.set_ylim([min(energies)-0.25, max(energies)+0.25])
            axs[im].set_ylim([-2,2])
            axs[im].legend(loc=1, frameon=False,fontsize=20) 
            axs[im].tick_params(labelsize=20)
    fig.text(0.03, 0.5, r'Energy (eV)', va='center', ha='center',rotation='vertical', fontsize=20)
    fig.text(0.5, 0.02, r'Time (ps)', va='center', ha='center',  fontsize=20)
    fig.savefig('output/convergence_energies_all.png', dpi=300, transparent =True)
