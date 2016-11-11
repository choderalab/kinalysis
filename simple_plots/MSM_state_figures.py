import matplotlib
matplotlib.use('Agg')

import numpy as np
import mdtraj as md
from msmbuilder import dataset
from glob import glob
import re

import mdtraj as md
import seaborn as sns
import numpy as np

import sys
import os

####
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")
sns.set_context("poster")

import json
import argparse

# run inputs.py before running this script

# DEFINE YOUR INPUTS

files = "trajectories-ck2/CK2*.pdb"
trajectories = dataset.MDTrajDataset(files)
protein = 'CK2'
project = '11406'

#### END DEFINE INPUTS ####


# Make shukla plot

with open('KER_hbond.json', 'r') as fp:
    KER_hbond = json.load(fp)
with open('Aloop_def.json', 'r') as fp:
    Aloop_def = json.load(fp)

def shukla_coords(trajectories,KER,Aloop,SRC2):

    difference = []
    rmsd = []

    for traj in trajectories:

        # append difference
        k295e310 = md.compute_contacts(traj, [KER[0]])
        e310r409 = md.compute_contacts(traj, [KER[1]])
        difference.append(10*(e310r409[0] - k295e310[0])) # 10x because mdtraj is naturally in nm

        # append rmsd
        Activation_Loop_SRC2 = SRC2.top.select("backbone and (resid %s to %s)" %(Aloop[0],Aloop[1]))
        Activation_Loop_kinase = traj.top.select("backbone and (resid %s to %s)" %(Aloop[0],Aloop[1]))

        SRC2_cut = SRC2.atom_slice(Activation_Loop_SRC2)
        traj_cut = traj.atom_slice(Activation_Loop_kinase)

        rmsd.append(10*(md.rmsd(traj_cut,SRC2_cut,frame=0))) # 10x because mdtraj is naturaly in nm

    return [rmsd, difference]

SRC2 = md.load("%s_2SRC_A.pdb" %protein)

[rmsd_states,difference_states] = shukla_coords(trajectories,KER_hbond[protein],Aloop_def[protein],SRC2)

import matplotlib
colors = matplotlib.cm.hsv(np.linspace(0, 1, len(rmsd_states)))

for i in range(len(rmsd_states)):
    plt.plot(rmsd_states[i],difference_states[i][:,0],'o',color=colors[i],alpha=0.4,label='state %s'%i)

plt.xlabel('RMSD Activation Loop ($\AA$)')
plt.ylabel('d(E310-R409) - d(K295-E310) ($\AA$)')
plt.ylim(-20,20)
plt.xlim(0,10)
plt.legend()
plt.title('%s sims - %s' % (protein,project) )

plt.savefig('plot_states_%s.png' %project)

# make dfg plot

with open('DFG.json', 'r') as fp:
    DFG = json.load(fp)

def DFG_dihedral_states(trajectories,def_DFG):

    dihedral = []

    for traj in trajectories:
        dihedral.append(md.compute_dihedrals(traj,[def_DFG]))

    return [dihedral]

[dihedral] = DFG_dihedral_states(trajectories, DFG[protein])

print dihedral

np.save('dihedral_states_%s.npy'%protein,dihedral)

import math
dihedral_rotate = [ np.array( [A-(2*math.pi) if A >= 1.9 else A for A in run ] ) for run in dihedral ]

import matplotlib
colors = matplotlib.cm.hsv(np.linspace(0, 1, len(dihedral)))

for i in range(len(dihedral)):
    for k in range(len(dihedral[i])):
        if k == 0:
            plt.axvline(dihedral_rotate[i][k],color=colors[i],alpha=0.4,label='state %s'%i)
        else:
            plt.axvline(dihedral_rotate[i][k],color=colors[i],alpha=0.4)

plt.xlabel('Dihedral (radians)')
plt.ylabel('Occupancy')
plt.xlim(-5,2)
plt.legend(loc=0)

plt.title('%s sims DFG states - %s'%(protein,project))

plt.savefig('plot_states_DFG_%s.png'%project)



