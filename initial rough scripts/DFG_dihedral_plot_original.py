import matplotlib
matplotlib.use('Agg')

import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from msmbuilder import dataset

# load trajectories

Abl_trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged/no-solvent/10472/*.h5")
Src_trajectories = dataset.MDTrajDataset("/cbio/jclab/projects/fah/fah-data/munged/no-solvent/10471/*.h5")

# define DFG dihedral ( this is from Roux umbrella sampling paper and are AlaCbeta, AlaCalpha, AspCalpha, AspCgamma)

Abl_DFG = [2257,2255,2265,2270]
Src_DFG = [2190,2188,2198,2203]

def DFG_dihedral(trajectories,def_DFG):

    dihedral = []

    for traj in trajectories:

        dihedral.append(md.compute_dihedrals(traj,[def_DFG]))

    flattened_dihedral = np.asarray([val for sublist in dihedral for val in sublist])

    return [flattened_dihedral]

[Abl_dihedral] = DFG_dihedral(Abl_trajectories, Abl_DFG)
[Src_dihedral] = DFG_dihedral(Src_trajectories, Src_DFG)

np.save('Abl_dihedral.npy',Abl_dihedral)
np.save('Src_dihedral.npy',Src_dihedral)

import math

Abl_rotate =  [A-(2*math.pi) if A >= 1.9 else A for A in Abl_dihedral]
Src_rotate =  [S-(2*math.pi) if S >= 1.9 else S for S in Src_dihedral]

sns.distplot(Abl_rotate, color="r",label="Abl")
sns.distplot(Src_rotate, color="b",label="Src")
plt.xlabel('Dihedral (radians)')
plt.ylabel('Occupancy')
plt.legend()

plt.savefig('abl_src_DFG_dihedral_hist.eps',format='eps',dpi=1000)

