# This script allows you to make plots according to the coordinates in Fig 2 of Shukla et al Nature Communications 2014.
#
# Right now the syntax is pretty silly but for now this is how it works:
#
#     python plotting_Shukla.py <project> <kinase> 
#  Example: python plotting_Shukla.py 11401 'SRC'
#
# The only kinases currently availabe are 'SRC', 'ABL', and 'DDR1', but it should be simple to add your own.
#
# Made by Sonya Hanson Jan 2016, for better or worse.

# import libraries
import matplotlib
matplotlib.use('Agg')

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import sys

from msmbuilder import dataset

import seaborn as sns
sns.set_style("whitegrid")
sns.set_context("poster")

# Define project.
project = sys.argv[1]

# Define kinase.
kinase = sys.argv[2]

# load trajectories
trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged/no-solvent/%s/*.h5" %project)

# import 2SRC structure to compare to
SRC2 = md.load("%s_2SRC_A.pdb" %kinase)

# Define hydrogen bond coordinates (0-indexed)
KER_hbond = { 'SRC' : [[28,43],[43,142]],
              'ABL' : [[29,44],[44,144]],
              'DDR1': [[51,68],[68,185]]}

# Define Activation loop (resid)
Aloop_def = { 'SRC': [138,158],
          'ABL': [140,160],
          'DDR1': [181,201]}

def shukla_coords(trajectories,KER,Aloop,SRC2):

    difference = []
    rmsd = []

    for traj in trajectories:

        # append difference
        k295e310 = md.compute_contacts(traj, [KER[0]])
        e310r409 = md.compute_contacts(traj, [KER[1]])
        difference.append(10*(e310r409[0] - k295e310[0])) # 10x because mdtraj is naturally in nm

        # append rmsd
        Activation_Loop_SRC2 = SRC2.top.select("backbone and (resid %s to %s)" %(138,158))
        Activation_Loop_kinase = traj.top.select("backbone and (resid %s to %s)" %(Aloop[0],Aloop[1]))

        SRC2_cut = SRC2.atom_slice(Activation_Loop_SRC2)
        traj_cut = traj.atom_slice(Activation_Loop_kinase)

        rmsd.append(10*(md.rmsd(traj_cut,SRC2_cut,frame=0))) # 10x because mdtraj is naturaly in nm

    # flatten list of arrays
    flattened_difference = np.asarray([val for sublist in difference for val in sublist])
    flattened_rmsd = np.asarray([val for sublist in rmsd for val in sublist])

    return [flattened_rmsd, flattened_difference]

# generate data
[rmsd,difference] = shukla_coords(trajectories,KER_hbond[kinase],Aloop_def[kinase],SRC2)

#save rmsd and difference data
np.save('rmsd_%s.npy'%kinase,rmsd)
np.save('difference_%s.npy'%kinase,difference)

#plot
sns.kdeplot(rmsd,difference[:,0],shade=True,log=True,cmap="Reds",shade_lowest=False)

plt.xlabel('RMSD Activation Loop ($\AA$)')
plt.ylabel('d(E310-R409) - d(K295-E310) ($\AA$)')
plt.ylim(-20,20)
plt.xlim(0,10)
plt.title('%s sims - %s' % (kinase,project) )

plt.savefig('plotting_Shukla_%s_%s.png' % (kinase,project),dpi=1000)
