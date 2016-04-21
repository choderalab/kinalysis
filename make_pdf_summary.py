
import mdtraj as md
import seaborn as sns
import numpy as np

import sys

from msmbuilder import dataset 
####
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")
sns.set_context("poster")

import kinalysis
import json

from glob import glob

protein = 'SRC'

#trajectory = md.load("run24-clone15.xtc",top="run24-clone15.pdb")

#files = sys.argv[1]
files = "trajectories/*.h5"

trajectories = dataset.MDTrajDataset(files)

### THE FIRST COUPLE PLOTS ARE ONLY CONDUCTED ON THE LONGEST SIM.
[max_length] = [max(len(traj) for traj in trajectories)]

for traj in trajectories:
    if len(traj) == max_length:
        longest_traj=traj

frame = np.arange(len(longest_traj))[:, np.newaxis]

# Using 0.25 so that units are in ns.
time = frame * longest_traj.timestep*1e-3

sim_time = time[-1] * 1e-3

print "This script is analyzing %s simulations." % len(trajectories)   
print "The longest simulation is %s us." % ''.join(map(str, sim_time)) 

### RMSD PLOT

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


### Inspired by http://stackoverflow.com/questions/8085520/generating-pdf-latex-with-python-script
### Make PDF Summary
#####
import argparse
import os
import subprocess

content = r'''\documentclass{article}
\usepackage[pdftex]{graphicx}
\usepackage[margin=0.65in]{geometry}
\usepackage{caption}
\usepackage{subcaption}
\renewcommand{\familydefault}{\sfdefault}
\begin{document}
\textbf{\huge PDF Summary of %(protein)s kinase \\}
\\

\begin{figure}[!hbp]
 \begin{subfigure}[b]{0.4\textwidth}
  \centering
  \includegraphics[width=\textwidth]{rmsd.png}
  \caption{RMSD from initial}
  \label{fig:sub1}
 \end{subfigure}
\hfill
 \begin{subfigure}[b]{0.4\textwidth}
  \centering
  \includegraphics[width=\textwidth]{dssp.png}
  \caption{Secondary Structure}
  \label{fig:sub2}
 \end{subfigure}
\caption{Length of the longest sim is %(time)s us.}
\label{fig:longest_sim}
\end{figure}

\noindent\rule{16cm}{0.4pt}
\begin{figure}[h]
 \centering
 \includegraphics[width=7cm]{dihedral.png}
 \caption{DFG dihedral histogram of all sims}
\end{figure}
\begin{figure}[h]
 \centering
 \includegraphics[width=7cm]{shukla.png}
 \caption{Shukla plot: KDE of all sims}
\end{figure}
\end{document}
'''

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--time', default='%s' % ''.join(map(str, sim_time)) )
parser.add_argument('-p', '--protein', default='%s' %protein) 
parser.add_argument('-s', '--school', default='My U')

args = parser.parse_args()

with open('pdf_summary.tex','w') as f:
    f.write(content%args.__dict__)

cmd = ['pdflatex', '-interaction', 'nonstopmode', 'pdf_summary.tex']
proc = subprocess.Popen(cmd)
proc.communicate()

retcode = proc.returncode
if not retcode == 0:
    os.unlink('pdf_summary.pdf')
    raise ValueError('Error {} executing command: {}'.format(retcode, ' '.join(cmd))) 

os.unlink('pdf_summary.tex')
os.unlink('pdf_summary.log')
