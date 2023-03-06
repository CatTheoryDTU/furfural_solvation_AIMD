#!/usr/bin/env python

import matplotlib.pyplot as  plt
import numpy as np
from ase.data import chemical_symbols as symbols
import sys
#sys.path.append("/Users/mianle/Documents/ClusterConnect/TOOLS &PACKAGE/SimSoliqTools")
from simsoliq.io import init_mdtraj
from simsoliq.plotting.standard_plots import plot_density
from simsoliq.analyze.density import isolate_solvent_density, get_peak_integral, \
    get_average_solvent_bulk_density
from simsoliq.mdtraj_average import average_densities

metals = ['Au','Cu','Pt','Rh']
#metals = ['Au']
runs = ['2','3']
ads = ['clean','furfural']


path1 = '/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/'
path2 = '/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/'
path3 = '/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/'
path4 = '/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/'
path5 = '/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/'
path6 = '/restart/restart/restart/restart/restart/restart/restart/restart/restart/'
path7 = '/restart/restart/restart/restart/restart/restart/restart/restart/'
path8 = '/restart/restart/restart/restart/restart/restart/restart/'

paths = [path1,path2,path3,path4,path5,path6,path7,path8]

if __name__ == "__main__":

    ##################################################################
    ### this example demonstrates the handling of density profiles ###
    ### with simsoliq:                                             ###
    ###     (1) automated creation of density plots of solvent     ###
    ###     (2) automated analysis of density of solvent           ###
    ###     (3) averaging of density from different trajectories   ###
    ###     (4) handling and plotting of raw densities from mdtraj ###
    ##################################################################


    

    ###     (3) averaging of density from different trajectories   ###
    ##################################################################
    
    # read trajectories (use restart as individual traj)
    for m in metals:
        fig = plt.figure(figsize=(8,6))
        plt.ylabel('density (mol/cm$^3$)',fontsize=20, labelpad=30)
        plt.xlabel('$\mathrm{\AA}$', fontsize=20, labelpad=30)
        plt.xticks([])
        plt.yticks([])
        gs = fig.add_gridspec(2, hspace=0)
        axs = gs.subplots(sharex=True, sharey=True)
        print(m)
        for ia, a in enumerate(ads):
            trj = []
            for r in runs:
                path0 = '/p/project/pra121/nitish/projects/1_solvation/'+m+'_111/'+a+'-40w-'+r+'/1e-5'
                for p in paths:
                    trj.append(init_mdtraj(path0+p+"vasprun.xml", fmat='vasp'))
            # (only one composition included here)
            av_dens = average_densities(trj, tstart=0)
            #print(len(av_dens))
            for comp in av_dens:
                c_dens = isolate_solvent_density(av_dens[comp])
                binc = c_dens['binc']
                hist = c_dens['hists']
                axs[ia].plot(binc,hist['Osolv'],color='r',lw=3, label='Ow')
                axs[ia].plot(binc,hist['Hsolv'],color='c',lw=3, label='Hw')
                #integrate the first water layer density
                x_range = [i for i, b in enumerate(binc) if b >= 2.5 and b <= 4.0]
                y_range = [hist['Osolv'][i] for i in x_range]
                y_min = min(y_range)
                y_min_index = y_range.index(y_min)
                x_min = x_range[y_min_index]
                fill_x = x_range[:y_min_index]
                y_range_o = hist['Osolv'][:y_min_index]
                y_range_h = hist['Hsolv'][:y_min_index]
                o_int = np.trapz(fill_x, y_range_o)
                h_int = np.trapz(fill_x, y_range_h)
                print('Integrated O and H density in 1st layer are '+str(round(o_int,6))+', '+str(round(h_int,6)))
                axs[ia].fill_between(fill_x,y_range_h,[x*0 for x in fill_x],color='c',alpha = 0.8)
                axs[ia].fill_between(fill_x,y_range_o,[x*0 for x in fill_x],color='r',alpha = 0.8)
                axs[ia].set_ylim(0,0.27)
                axs[ia].set_xlim(0,15.0)
                axs[ia].tick_params(labelsize=15)
                if ia == 0:
                    axs[ia].legend(frameon=False, fontsize=15)
        fig.tight_layout()
        fig.savefig('%s_density_integral_1st_wlayer.png'%m)
                
                

    

