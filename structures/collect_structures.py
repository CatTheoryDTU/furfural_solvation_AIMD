#!/usr/bin/python

from ase import Atom, Atoms
from ase.io import read, write
from pathlib import Path
import os

metals = ['Au','Cu','Pd','Pt','Rh']
ads = ['clean','fur']
runs = range(1,4)

raw_dir = '/p/project/pra121/nitish/projects/1_solvation/'

for m in metals:
    for a in ads:
        for r in runs:
            str_dir = raw_dir + m + '_111/' + a + '-40w-' + str(r) + '/1e-5/'
            atoms = read(str_dir + 'POSCAR')
            atoms.write(m + '/' + m + '_' + a + '_' + str(r) + '.cif')




