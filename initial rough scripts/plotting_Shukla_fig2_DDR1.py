# import libraries

import matplotlib
matplotlib.use('Agg')

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import glob

from msmbuilder import dataset

import seaborn as sns
sns.set_style("whitegrid")
sns.set_context("poster")

# load trajectories

#trajectories = dataset.MDTrajDataset("trajectories/*.h5")

# Load trajectory with ensembler models

t_models = md.load("models/full_ensembler/modelstraj.xtc", top = "models/full_ensembler/modelstraj-topol.pdb")

# Load new starting models

#start = md.load(glob.glob("models/*.pdb"))

active_files = ["models/1IR3.pdb","models/1K3A.pdb","models/1Y57.pdb","models/2F4J.pdb","models/3DQW.pdb"]
abllike_files = ["models/4BKJ.pdb","models/3ZOS.pdb","models/2OIQ.pdb","models/1OPJ.pdb"]
srclike_files = ["models/2SRC.pdb","models/2G2I.pdb","models/2H8H.pdb","models/3U4W.pdb"]
other_files = ["models/DDR1_Dasatinib.pdb","models/DDR1_VX680.pdb","models/4GT5.pdb","models/1IRK.pdb"]

active = md.load(active_files)
abllike = md.load(abllike_files)
srclike = md.load(srclike_files)
other = md.load(other_files)

# the coordinates for the 'full ensembler' models and the old trajectories
#  are shorter so there are two versions of the numbering.

# Load 2SRC reference structures

SRC2 = md.load("models/2SRC.pdb")
SRC2_short = md.load("/home/hansons/sims/ddr1/10484/reference-structures/DDR1_2SRC_A.pdb")

# Define hydrogen bond coordinates (0-indexed)

KER = [[51,68],[68,185]]
KER_short = [[45,62],[62,179]]

# Define Activation loop (resid)

Aloop = [181,201]
Aloop_short = [175,195]

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
        Activation_Loop_Src = traj.top.select("backbone and (resid %s to %s)" %(Aloop[0],Aloop[1]))

        SRC2_cut = SRC2.atom_slice(Activation_Loop_SRC2)
        traj_cut = traj.atom_slice(Activation_Loop_Src)

        rmsd.append(10*(md.rmsd(traj_cut,SRC2_cut,frame=0))) # 10x because mdtraj is naturaly in nm

    # flatten list of arrays

    flattened_difference = np.asarray([val for sublist in difference for val in sublist])
    flattened_rmsd = np.asarray([val for sublist in rmsd for val in sublist])

    return [flattened_rmsd, flattened_difference]

# generate data

[rmsd_act,difference_act] = shukla_coords(active,KER,Aloop,SRC2)
[rmsd_abl,difference_abl] = shukla_coords(abllike,KER,Aloop,SRC2)
[rmsd_src,difference_src] = shukla_coords(srclike,KER,Aloop,SRC2)
[rmsd_oth,difference_oth] = shukla_coords(other,KER,Aloop,SRC2)

[rmsd,difference] = shukla_coords(t_models,KER_short,Aloop_short,SRC2_short)

#plot
plt.plot(rmsd, difference[:,0], 'o', markersize=2, label="ensembler models", color='grey')
sns.kdeplot(rmsd,difference[:,0],shade=True,log=True, cut=10)
plt.plot(rmsd_act, difference_act[:,0], '*',markersize=20, label="active",color='yellow')
plt.plot(rmsd_abl, difference_abl[:,0], '*',markersize=20, label="abl-like inactive",color='green')
plt.plot(rmsd_src, difference_src[:,0], '*',markersize=20, label="src-like inactive",color='blue')
plt.plot(rmsd_oth, difference_oth[:,0], '*',markersize=20, label="other inactive",color='red')


plt.xlabel('RMSD Activation Loop ($\AA$)')
#Equivalent to E310-R409 and K295-E310 in SRC
#Numbering according to DDR1b
plt.ylabel('d(E672-R789) - d(K655-E672) ($\AA$)')
plt.ylim(-20,20)
plt.xlim(0,10)
plt.title('DDR1 new starting structures')
plt.legend()

plt.savefig('plotting_Shukla_fig2_DDR1-colors.png')
