import os
import numpy as np
from ase import db
from ase.io import read, write
from ase import Atoms

def init_traj(filepath, filename, fmat='vasp'):
    fpath = filepath
    fname = filename ##e.g., vasprun.xml
    raw_data = read(fpath+'/'+fname,':')
    return raw_data

def get_raw_traj(raw_data):
    raw_traj=[]
    for rd in raw_data:
        raw_traj.append(rd.get_positions())
    return raw_traj

def count_time(raw_traj, step):
    runtime = len(raw_traj)*step
    return runtime

def get_top_slab_mean_z_positions(raw_traj, N_s, n_s):
    mean_z_top = []
    for traj in raw_traj:
        surf_pos_z = traj[:N_s,2]
        surf_pos_z_top = np.sort(surf_pos_z)[-1*n_s:]
        mean_z_top.append(np.mean(surf_pos_z_top))
    return mean_z_top

def get_solvent_z_positions(raw_traj, N_w, N_s): ###use water here
    sol_pos_z = []
    for traj in raw_traj:
        sol_pos_z.append(traj[N_s:(N_s+N_w),2])
    return sol_pos_z

def get_ads_mean_z_positions(raw_traj, N_ads):
    mean_z_ads = []
    for traj in raw_traj:
        ads_pos_z = traj[-1*N_ads:,2]
        mean_z_ads.append(np.mean(ads_pos_z))
    return mean_z_ads


def get_adsorbed_solvent(raw_traj, n_s, dist): ##N_w_ads/N_s = N_w_ads site-1
    N_w_ads = []
    sol_pos_z = get_solvent_z_positions(raw_traj,N_w,N_s)
    mean_z_top = get_top_slab_mean_z_positions(raw_traj,N_s, n_s)
    for i, spz in enumerate(sol_pos_z):
        count = 0
        for s in spz:
            if np.abs(s - mean_z_top[i]) <= dist:
                count += 1
            else:
                count += 0
        N_w_ads.append(count)
    #print(N_w_ads)
    return np.mean(N_w_ads)/n_s

def get_hbonds_ads(raw_traj,raw_data, N_s, N_w, N_ads, h_dist): #e.g., h_dist = 3 for HO--H
    hbonds_ads = []
    for traj in raw_traj:
        ads_traj = traj[-1*N_ads:]
        ads_symb = raw_data[0].get_chemical_symbols()[-1*N_ads:] # list(ads_traj.symbols)
        #list_ads = [i for i, elem in enumerate(ads_elem) if elem in ('H', 'O')]
        O_traj = traj[N_s:N_s+N_w]
        H_traj = traj[N_s+N_w:N_s+2*N_w]
        count = 0
        for i, symb in enumerate(ads_symb):
            if symb == 'O':
                for h_pos in H_traj:
                    if h_dist - 0.05 <= np.linalg.norm(ads_traj[i] - h_pos) <= h_dist + 0.05:
                        count += 1
            if symb == 'H':
                for o_pos in O_traj:
                    if h_dist - 0.05 <= np.linalg.norm(ads_traj[i] - o_pos) <= h_dist + 0.05:
                        count += 1
            else:
                count += 0
        hbonds_ads.append(count)
    #print(hbonds_ads)
    return np.mean(hbonds_ads)

##main script
#fp = '/p/project/pra121/nitish/projects/1_solvation/Rh_111/clean-40w-1/1e-5/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart'
#fn = 'vasprun.xml'

surfaces = ['Au_111','Cu_111','Pt_111','Rh_111']
N_s, n_s, N_w, N_ads, dist, h_dist = 64, 16, 40, 11, 2.50, 3.0

for s in surfaces:
    for i in range(1,4):
        fp = '/p/project/pra121/nitish/projects/1_solvation/'+s+'/clean-40w-'+str(i)+'/1e-5/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart'
        fn = 'vasprun.xml'
        raw_data = init_traj(fp, fn)
        print('Runtime is '+ str(count_time(raw_data, 1))+' fs')
        raw_traj = get_raw_traj(raw_data)[:]
        output = get_adsorbed_solvent(raw_traj, n_s, dist)
        print('Cumulative Number of Adsorbed Water on '+s+'_'+str(i)+' is '+str(round(output,3))+' /site')

paths = ['/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart',
        '/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart',
        '/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart']

for s in surfaces:
    for i in range(0,3):
        fp = '/p/project/pra121/nitish/projects/1_solvation/'+s+'/fur-40w-3/1e-5'+paths[i]
        fn = 'vasprun.xml'
        raw_data = init_traj(fp, fn)
        print('Runtime is '+ str(count_time(raw_data, 1))+' fs')
        raw_traj = get_raw_traj(raw_data)[:]
        output = get_adsorbed_solvent(raw_traj, n_s, dist)
        print('Cumulative Number of Adsorbed Water with adsorbate on '+s+'_'+str(i)+' is '+str(round(output,4))+' /site')



#raw_data = init_traj(fp, fn)
#runtime = count_time(raw_data, 1)
#print('Runtime is '+ str(runtime)+' fs')

#raw_traj = get_raw_traj(raw_data)[:]
#N_s, n_s, N_w, N_ads, dist, h_dist = 64, 16, 40, 11, 2.55, 3.0

#mean_z_top = get_top_slab_mean_z_positions(raw_traj, N_s, n_s)
#print('mean z positions of top surface is '+str(mean_z_top))

#sol_pos_z = get_solvent_z_positions(raw_traj,N_w,N_s)
#print(sol_pos_z[0])

#output = get_adsorbed_solvent(raw_traj, n_s, dist)
#print('Cumulative Number of Adsorbed Water is '+str(round(output,3))+' /site')

#hbonds_ads = get_hbonds_ads(raw_traj,raw_data, N_s, N_w, N_ads, h_dist)
#print('Cumulative Number of hbonds for adsorbate is '+str(hbonds_ads)+' /ads')
