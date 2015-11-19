# import libraries

import matplotlib
matplotlib.use('Agg')

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np

from msmbuilder import dataset

import seaborn as sns
sns.set_style("whitegrid")
sns.set_context("poster")

# load trajectories

trajectories = dataset.MDTrajDataset("../ipynbs/trajectories/*.h5")

WIG = md.load("../original-models/3WIG_model.pdb")
AN2 = md.load("../original-models/4AN2_model.pdb")
EQD = md.load("../original-models/3EQD_model.pdb")
EQI = md.load("../original-models/3EQI_model.pdb")
EQG = md.load("../original-models/3EQG_model.pdb")
ORN = md.load("../original-models/3ORN_model.pdb")

def catkhrd(trajectories):

     # define empty lists

     D218 = []
     D222 = []

     for traj in trajectories:

          #append h188s218 difference

          h188s218 = md.compute_contacts(traj, [[120,151]],scheme='ca')
          D218.append(h188s218[0])

          #append k97s222 difference

          k97s222 = md.compute_contacts(traj, [[29,155]],scheme='ca')
          D222.append(k97s222[0])

     #flatten these lists of arrays

     flattened_h188s218 = np.asarray([val for sublist in D218 for val in sublist])
     flattened_k97s222 = np.asarray([val for sublist in D222 for val in sublist])

     return [flattened_h188s218, flattened_k97s222]
      
#plot

[s1_traj,s2_traj] = catkhrd(trajectories)

[s1_3WIG,s2_3WIG] = catkhrd(WIG)
[s1_4AN2,s2_4AN2] = catkhrd(AN2)
[s1_3EQD,s2_3EQD] = catkhrd(EQD)
[s1_3EQI,s2_3EQI] = catkhrd(EQI)
[s1_3EQG,s2_3EQG] = catkhrd(EQG)
[s1_3ORN,s2_3ORN] = catkhrd(ORN)

sns.kdeplot(s1_traj[:,0],s2_traj[:,0],shade=True,log=True)

plt.plot(s1_3WIG[:,0],s2_3WIG[:,0], '*',markersize=20, label="3WIG",color='yellow')
plt.plot(s1_4AN2[:,0],s2_4AN2[:,0], '*',markersize=20, label="4AN2(EE)",color='black')
plt.plot(s1_3EQD[:,0],s2_3EQD[:,0], '*',markersize=20, label="3EQD",color='blue')
plt.plot(s1_3EQI[:,0],s2_3EQI[:,0], '*',markersize=20, label="3EQI",color='red')
plt.plot(s1_3EQG[:,0],s2_3EQG[:,0], '*',markersize=20, label="3EQG",color='magenta')
plt.plot(s1_3ORN[:,0],s2_3ORN[:,0], '*',markersize=20, label="3ORN",color='cyan')

plt.xlabel('d(H188-M219) (nm)')
plt.ylabel('d(K97-F223) (nm)')
plt.legend()

plt.savefig('plot_hrd_151_catk_155_heatmap_stars.png')