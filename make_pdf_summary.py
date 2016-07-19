import matplotlib
matplotlib.use('Agg')

import mdtraj as md
import seaborn as sns
import numpy as np

import sys
import os

from msmbuilder import dataset 
####
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")
sns.set_context("poster")

import kinalysis
import json

from glob import glob
import argparse

# Define argparse stuff

parser = argparse.ArgumentParser(description="""Make initial plots of your kinase simulations:  > python make_summary_pdf.py --proje
ct 11401""")
parser.add_argument("--project", help="the FAH project you want to analyze, e.g. 11401", action="store", default=False, type=int)
parser.add_argument("--byruns", help="use this if you want to add the plots by run", action="store_true", default=False)
args = parser.parse_args()

with open('projects.json', 'r') as fp:
    projects = json.load(fp)

if args.project:
    myproject = args.project
    protein =  projects.keys()[projects.values().index(myproject)]
    files = "/cbio/jclab/projects/fah/fah-data/munged3/no-solvent/%s/*.h5" %args.project
    
    # make new results folder with project name
    newpath = "./results/%s/" %args.project
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    os.chdir(newpath)

    if args.byruns:
        print "*** kinalysis: analyzing project %s (%s) BY RUNS ***" % (args.project,protein)
    else:
        print "*** kinalysis: analyzing project %s (%s) ***" % (args.project,protein)
else:
    myproject = 'no project'
    protein = 'SRC'
    files = "trajectories/*.h5"

# Define our trajectories

trajectories = dataset.MDTrajDataset(files)

### LETS FIND OUT SOME THINGS ABOUT ALL OF OUR TRAJECTORIES.
print "This script is analyzing %s simulations." % len(trajectories)
sim_num = len(trajectories)

lens = []
max_length = 0
for i,traj in enumerate(trajectories):
       if len(traj) > max_length:
          max_length = len(traj)
       if len(traj) == max_length:
          longest_traj=traj
       print i
       lens.append(len(traj))

frame = np.arange(len(longest_traj))[:, np.newaxis]

time = frame * longest_traj.timestep*1e-3

all_frames = sum(lens)

total_time = all_frames * longest_traj.timestep*1e-3

sim_time = time[-1] * 1e-3
total_sim_time = total_time* 1e-3

print "Totalling %s us." %  total_sim_time
print "The longest simulation is %s us." % ''.join(map(str, sim_time))

plt.hist(lens,bins=50);
plt.xlabel('Trajectory length')
plt.ylabel('Occurrences')
plt.title('Distribution of trajectory lengths (frames)')

plt.tight_layout()
plt.savefig('traj_distribution.png')

### THE FIRST COUPLE PLOTS ARE ONLY CONDUCTED ON THE LONGEST SIM.

### RMSD PLOT

plt.clf()

rmsd = md.rmsd(longest_traj,longest_traj,frame=0)
plt.plot(time, rmsd)
plt.xlabel('time (ns)')
plt.ylabel('RMSD(nm)')
plt.title('RMSD');

plt.tight_layout()
plt.savefig('rmsd.png')

### SECONDARY STRUCTURE PLOT

plt.clf()

dssp = md.compute_dssp(longest_traj)
dssp_counts = []
for d in dssp:
    unique, counts = np.unique(d, return_counts=True)
    dssp_counts.append(dict(zip(unique, counts)))

total_vals = sum(dssp_counts[0].values())
C_values = []
for d in dssp_counts:
    C_values.append(d['C']/float(total_vals))
E_values = []
for d in dssp_counts:
    E_values.append(d['E']/float(total_vals))
H_values = []
for d in dssp_counts:
    H_values.append(d['H']/float(total_vals))

plt.plot(time,E_values,label='beta');
plt.plot(time,H_values,label='helix');
plt.plot(time,C_values,label='neither');
plt.xlabel('time (ns)')
plt.ylabel('Percent Occupancy')
plt.ylim(0,1)
plt.legend();

plt.tight_layout()
plt.savefig('dssp.png')

### Now lets analyze all the trajectories.
### DFG PLOT

plt.clf()

with open('DFG.json', 'r') as fp:
    DFG = json.load(fp)

[Src_dihedral] = kinalysis.DFG_dihedral(trajectories, DFG[protein])

#Save dihedral
np.save('dihedral_%s.npy'%protein,Src_dihedral)

sns.distplot(Src_dihedral, color="b",label=protein)
plt.xlabel('Dihedral (radians)')
plt.ylabel('Occupancy')
plt.legend()

plt.tight_layout()
plt.savefig('dihedral.png')

### MAKE SHUKLA PLOT

plt.clf()

with open('KER_hbond.json', 'r') as fp:
    KER_hbond = json.load(fp)
with open('Aloop_def.json', 'r') as fp:
    Aloop_def = json.load(fp)

SRC2 = md.load("reference_pdbs/%s_2SRC_A.pdb" %protein)

[rmsd,difference] = kinalysis.shukla_coords(trajectories,KER_hbond[protein],Aloop_def[protein],SRC2)

#save rmsd and difference data
np.save('rmsd_%s.npy'%protein,rmsd)
np.save('difference_%s.npy'%protein,difference)

sns.kdeplot(rmsd,difference[:,0],shade=True,log=True,cmap="Reds",shade_lowest=False)

plt.xlabel('RMSD Activation Loop ($\AA$)')
plt.ylabel('d(E310-R409) - d(K295-E310) ($\AA$)')
plt.ylim(-20,20)
plt.xlim(0,10)
plt.title('%s sims' % (protein) )

plt.savefig('shukla.png' )

#make PDF summary

kinalysis.pdf_summary(protein,sim_num,total_sim_time,sim_time,myproject)

