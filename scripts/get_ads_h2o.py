import os
import numpy as np
from ase import db
from ase.io import read, write
from ase import Atoms

#functions
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

def get_elements(raw_data):
    all_elements = raw_data[0].symbols
    return all_elements

def count_time(raw_traj, step):
    runtime = len(raw_traj)*step
    return runtime

def get_top_slab_mean_z_positions(traj, n_s, slab_list):###
    slab_traj = [traj[i] for i in slab_list]
    surf_pos_z = np.array(slab_traj)[:,2] #get z coordinates for all slab atoms
    surf_pos_top_z = np.sort(surf_pos_z)[-1*n_s:] # get the top-layer slab z coordinates
    mean_top_z = np.mean(surf_pos_top_z)
    return mean_top_z

def get_solvent_traj(raw_data, traj, ads_list, slab_list):
    all_elements = get_elements(raw_data)
    solvent_traj = [traj[i] for i in np.arange(len(traj)) if i not in slab_list if i not in ads_list]   
    solvent_symb = [all_elements[i] for i in np.arange(len(all_elements)) if i not in slab_list if i not in ads_list]   
    return solvent_traj, solvent_symb

def get_h2o_separate_traj(raw_data, traj, N_w, ads_list, slab_list):
    solvent_traj, solvent_symb = get_solvent_traj(raw_data, traj, ads_list, slab_list)
    
    o_h2o_traj = []
    h_h2o_traj = []
    for i, s in enumerate(solvent_symb):
        if s == 'O':
            o_h2o_traj.append(solvent_traj[i])
        else:
            h_h2o_traj.append(solvent_traj[i])

    return solvent_traj, o_h2o_traj, h_h2o_traj

def get_adsorbed_h2o(raw_data,traj, N_w, n_s, ads_list, slab_list, dist): ##N_w_ads/N_s = N_w_ads site-1
    """
    retrieve the cumulative number of adsorbed solvents by calculating the
    z-direction distance between solvent molecule and top-slab plane;
    here, we use O in H2O to estimate the position of water molecules;
    """
    #print('Start count the adsorbed waters!')
    solvent_traj, o_h2o_traj, h_h2o_traj = get_h2o_separate_traj(raw_data, traj, N_w, ads_list, slab_list)
    sol_pos_z = [t[2] for t in o_h2o_traj]
    mean_top_z = get_top_slab_mean_z_positions(traj, n_s, slab_list)
    count = 0
    for spz in sol_pos_z:
        if abs(abs(spz - mean_top_z) - dist) <= 0.05:
            count += 1
        else:
            count += 0
    return count/n_s ##normalized to per site

##########
surfaces = ['Pd_111']
#surfaces = ['Au_111','Cu_111','Pt_111','Rh_111']
N_s, n_s, N_w, N_ads, dist, h_dist = 64, 16, 40, 11, 2.5, 2.55

paths = [#'/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart',]
        '/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart',]
        #'/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart',
        #'/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart',
        #'/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart/restart']

##w/o furfural
for s in surfaces:
    #break
    slab_list = np.arange(N_s)
    ads_list = []
    for i in [1,3]:
        fp_0 = '/p/project/pra121/nitish/projects/1_solvation/'+s+'/clean-40w-'+str(i)+'/1e-5'
        fn = 'vasprun.xml'
        raw_data = []
        #make a large init_traj of all sampled paths
        for p in paths:
            fp = fp_0 + p
            raw_data += init_traj(fp, fn)
        print('Runtime is '+ str(count_time(raw_data, 1))+' fs')
        raw_traj = get_raw_traj(raw_data)[:]
        output_all = []
        for traj in raw_traj:
            N_w_ads = get_adsorbed_h2o(raw_data,traj, N_w, n_s, ads_list, slab_list, dist)
            output_all.append(N_w_ads)
        output_avg = np.mean(output_all)
        print('Cumulative Number of Adsorbed Water on '+s+'_'+str(i)+' is '+str(round(output_avg,4))+' /site')

##w/ furfural
for s in surfaces:
    break
    slab_list = np.arange(N_s)
    ads_list = list(np.arange(N_s, N_s+2))+list(np.arange(N_s+2+N_w,N_s+2+N_w+4))+list(np.arange(N_s+2+N_w+4+N_w*2,N_s+2+N_w+4+N_w*2+5)) ##C5H4O2
    for i in range(2,3):
        fp_0 = '/p/project/pra121/nitish/projects/1_solvation/'+s+'/fur-40w-'+str(i)+'/1e-5'
        fn = 'vasprun.xml'
        raw_data = []
        #make a large init_traj of all sampled paths
        for p in paths:
            fp = fp_0 + p
            raw_data += init_traj(fp, fn)
        print('Runtime is '+ str(count_time(raw_data, 1))+' fs')
        raw_traj = get_raw_traj(raw_data)[:]
        output_all = []
        for traj in raw_traj:
            N_w_ads = get_adsorbed_h2o(traj, N_w, n_s, ads_list, slab_list, dist)
            output_all.append(N_w_ads)
        output_avg = np.mean(output_all) 
        print('Cumulative Number of Adsorbed Water with adsorbate on '+s+'_'+str(i)+' is '+str(round(output_avg,4))+' /site')


