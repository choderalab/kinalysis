
import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

import sys

from msmbuilder import dataset 
####

import seaborn as sns
sns.set_style("whitegrid")
sns.set_context("poster")

import kinalysis
import json

protein = 'SRC'

trajectory = md.load("run24-clone15.xtc",top="run24-clone15.pdb")

#files = sys.argv[1]
#trajectories = dataset.MDTrajDataset(files)

### Length sim.
frame = np.arange(len(trajectory))[:, np.newaxis]

# Using 0.25 so that units are in ns.
time = frame * trajectory.timestep*1e-3

sim_time = time[-1] * 1e-3

print "Length of this sim is %s us." % ''.join(map(str, sim_time)) 

### RMSD PLOT

rmsd = md.rmsd(trajectory,trajectory,frame=0)
plt.plot(time, rmsd)
plt.xlabel('time (ns)')
plt.ylabel('RMSD(nm)')
plt.title('RMSD')

plt.tight_layout()
plt.savefig('rmsd.png')

### DFG PLOT

plt.clf()

with open('DFG.json', 'r') as fp:
    DFG = json.load(fp)

[Src_dihedral] = kinalysis.DFG_dihedral(trajectory, DFG[protein])

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

SRC2 = md.load("Shukla_plot/%s_2SRC_A.pdb" %protein)

[rmsd,difference] = kinalysis.shukla_coords(trajectory,KER_hbond[protein],Aloop_def[protein],SRC2)

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
\usepackage[margin=0.75in]{geometry}
\renewcommand{\familydefault}{\sfdefault}
\begin{document}
\textbf{\huge PDF Summary of %(protein)s kinase \\}
\begin{figure}[h]
\includegraphics[width=7cm]{rmsd.png}
\end{figure}
Length of this sim is %(time)s us.
\begin{figure}[h]
\includegraphics[width=7cm]{dihedral.png}
\end{figure}
\begin{figure}[h]
\includegraphics[width=7cm]{shukla.png}
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
