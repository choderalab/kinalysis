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
    newpath = "./results/%s" %args.project
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    if args.byruns:
        print "*** kinalysis: analyzing project %s (%s) BY RUNS ***" % (args.project,protein)
    else:
        print "*** kinalysis: analyzing project %s (%s) ***" % (args.project,protein)
else:
    myproject = 'no project'
    protein = 'SRC'
    files = "trajectories/*.h5"
    newpath = "./results/%s" %protein
    if not os.path.exists(newpath):
        os.makedirs(newpath)

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

# DFG BYRUNS PLOT
# DFG dihedral by runs, finds out whether your sim started in DFG-in or DFG-out and then plots these two DFG histograms on the same plot.

plt.clf()

[dihedral] = kinalysis.DFG_dihedral_byrun(files, DFG[protein])

#Save dihedral
np.save('dihedral_separate_%s.npy'%protein,dihedral)

# Load dihedral: Abl_dihedral = np.load('dihedral_separate_Abl.npy')

# Calculate dihedral at the first frame of clone0 of each run

path_base = files.split('*')[0]
clone0_trajectories = dataset.MDTrajDataset("%s/*clone0.h5" % path_base )
firstframes = [traj[0] for traj in clone0_trajectories]
[lines] = kinalysis.DFG_dihedral(firstframes, DFG[protein])

# Rotate dihedral so histogram doesn't get cut in figure
import math

dihedral_rotate = [ np.array( [A-(2*math.pi) if A >= 1.9 else A for A in run ] ) for run in dihedral[0] ]
line_rotate =  [A-(2*math.pi) if A >= 1.9 else A for A in lines]

# Define which sims start in DFG-in vs DFG-out conformation
line_rotate = np.asarray(line_rotate)

DFG_in = np.where(line_rotate > -0.5)
DFG_out = np.where(line_rotate < -0.5)

# Accumulate these in a loop
dihedral_in = [ dihedral_rotate[index] for index in DFG_in[0] ]
dihedral_out = [ dihedral_rotate[index] for index in DFG_out[0] ]

dihedral_in_flattened = [val for sublist in dihedral_in for val in sublist]
dihedral_out_flattened = [val for sublist in dihedral_out for val in sublist]

#Plot just vertical lines at dihedral for first frames

try:
    line_rotate_out = np.vstack([ line_rotate[index] for index in DFG_out[0]])
except:
    line_rotate_out = []

try:
    line_rotate_in = np.vstack([ line_rotate[index] for index in DFG_in[0]])
except:
    line_rotate_in = []

for i in range(len(line_rotate_in)):
        plt.axvline(line_rotate_in[i], color="m", ymin=0.95)
for i in range(len(line_rotate_out)):
        plt.axvline(line_rotate_out[i], color="r", ymin=0.95)

# Plot histogram with special seaborn sauce
try:
    sns.distplot(dihedral_in_flattened, color="m",label="DFG in (%s) " %len(DFG_in[0]) )
except:
    pass

try:
    sns.distplot(dihedral_out_flattened, color="r",label="DFG out (%s) " %len(DFG_out[0]) )
except:
    pass

plt.xlabel('Dihedral (radians)')
plt.ylabel('Occupancy')
plt.legend()

plt.savefig('%s_DFG_dihedral_inout_lines.png'%protein)


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

plt.savefig('shukla.png')

### MAKE SHUKLA PLOT BY RUNS

[rmsd_separate,difference_separate] = kinalysis.shukla_coords_byrun(files,KER_hbond[protein],Aloop_def[protein],SRC2)

#save rmsd and difference data
np.save('rmsd_separate_%s.npy'%protein,rmsd_separate)
np.save('difference_separate_%s.npy'%protein,difference_separate)

import matplotlib
colors = matplotlib.cm.hsv(np.linspace(0, 1, len(rmsd_separate)))

for i in range(len(rmsd_separate)):
    plt.plot(rmsd_separate[i],difference_separate[i],'.',color=colors[i],alpha=0.15,ms=2)

for i in range(len(rmsd_separate)):
    try:
        plt.plot(rmsd_separate[i][0],difference_separate[i][0][0],'*',color=colors[i])
    except:
        pass

plt.xlabel('RMSD Activation Loop ($\AA$)')
plt.ylabel('d(E310-R409) - d(K295-E310) ($\AA$)')
plt.ylim(-20,20)
plt.xlim(0,10)
plt.title('%s sims' % (protein) )

plt.savefig('shukla_byrun_%s.png' %protein,dpi=700)


#make PDF summary

kinalysis.pdf_summary(protein,sim_num,total_sim_time,sim_time,myproject)

# Move all files to results folder.

print '***Finished! Now moving...***'

import shutil

for file in glob('*.png'):
    print "...%s"%file
    shutil.move(file,newpath)

for file in glob('*.npy'):
    print "...%s"%file
    shutil.move(file,newpath)

for file in glob('pdf_summary.*'):
    print "...%s"%file
    shutil.move(file,newpath)

print '***to %s!!!***' %newpath

