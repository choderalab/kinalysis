
"""
Tools for assisting in rough initial analysis of kinase simulations..
"""

#=============================================================================================
# Imports
#=============================================================================================

import numpy as np
import mdtraj as md
from msmbuilder import dataset
from glob import glob
import re

#=============================================================================================
# Getting the DFG dihedral from simulations
#=============================================================================================


def DFG_dihedral(trajectories,def_DFG):

    dihedral = []

    for traj in trajectories:
        try:
            dihedral.append(md.compute_dihedrals(traj,[def_DFG]))
        except:
            continue

    flattened_dihedral = np.asarray([val for sublist in dihedral for val in sublist])

    return [flattened_dihedral]

def DFG_dihedral_byrun(files,def_DFG):

    # Since we are going to sort files by where they are in first frame of clone0
    #   we can only analyze trajectories with a clone0 present.

    path_base = files.split('*')[0]

    clone0_files = "%s/*clone0.h5" % path_base

    globfiles = glob(clone0_files)

    runs_list = []

    for filename in globfiles:
        run_string = re.search('run([^-]+)',filename).group(1)
        run = int(run_string)
        if run not in runs_list:
            runs_list.append(run)
        runs_list.sort()

    dihedral = []
    dihedral_combinetrajs = []

    for run in runs_list:

        trajectories = dataset.MDTrajDataset("%s/run%d-*.h5" % (path_base,run))
        print "Run %s has %s trajectories." % (run,len(trajectories))

        for traj in trajectories:

            dihedral_combinetrajs.append(md.compute_dihedrals(traj,[def_DFG]))
        # flatten
        dihedral_combinetrajs = [val for sublist in dihedral_combinetrajs for val in sublist]

        dihedral.append(dihedral_combinetrajs)
        dihedral_combinetrajs = []

    dihedral = np.asarray([dihedral])

    return [dihedral]


#=============================================================================================
# Getting the Shukla Coordinates (salt bridge and Aloop RMSD)  from simulations
#=============================================================================================

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

    # flatten list of arrays
    flattened_difference = np.asarray([val for sublist in difference for val in sublist])
    flattened_rmsd = np.asarray([val for sublist in rmsd for val in sublist])

    return [flattened_rmsd, flattened_difference]

# Shukla by runs plots a star at frame 0 of clone 0 and then dots on the shukla plot for the trajectories of these runs


def shukla_coords_byrun(files,KER,Aloop,SRC2):

    difference = []
    rmsd = []

    difference_combinetrajs = []
    rmsd_combinetrajs = []

    path_base = files.split('*')[0]
    clone0_files = "%s/*clone0.h5" % path_base
    globfiles = glob(clone0_files)

    runs_list = []

    for filename in globfiles:
        run_string = re.search('run([^-]+)',filename).group(1)
        run = int(run_string)
        if run not in runs_list:
            runs_list.append(run)
        runs_list.sort()


    for run in runs_list:

        trajectories = dataset.MDTrajDataset("%s/run%d-clone*1.h5" % (path_base,run))
        print "Run %s has %s trajectories." % (run,len(trajectories))

        for traj in trajectories:

            # append difference
            k295e310 = md.compute_contacts(traj, [KER[0]])
            e310r409 = md.compute_contacts(traj, [KER[1]])
            difference_combinetrajs.append(10*(e310r409[0] - k295e310[0])) # 10x because mdtraj is naturally in nm

            # append rmsd
            Activation_Loop_SRC2 = SRC2.top.select("backbone and (resid %s to %s)" %(Aloop[0],Aloop[1]))
            Activation_Loop_kinase = traj.top.select("backbone and (resid %s to %s)" %(Aloop[0],Aloop[1]))

            SRC2_cut = SRC2.atom_slice(Activation_Loop_SRC2)
            traj_cut = traj.atom_slice(Activation_Loop_kinase)

            rmsd_combinetrajs.append(10*(md.rmsd(traj_cut,SRC2_cut,frame=0))) # 10x because mdtraj is naturaly in nm

        # flatten list of arrays
        difference_combinetrajs = np.asarray([val for sublist in difference_combinetrajs for val in sublist])
        rmsd_combinetrajs = np.asarray([val for sublist in rmsd_combinetrajs for val in sublist])

        difference.append(difference_combinetrajs)
        difference_combinetrajs = []

        rmsd.append(rmsd_combinetrajs)
        rmsd_combinetrajs = []

    return [rmsd, difference]



#=============================================================================================
# Making a PDF Summary of the results
#=============================================================================================

### Inspired by http://stackoverflow.com/questions/8085520/generating-pdf-latex-with-python-script

def pdf_summary(protein,sim_num,sim_tot,sim_time,project):

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
   \textbf{{\huge PDF Summary of %(protein)s kinase} Project: %(project)s \\}
   \\

   \begin{figure}[!hbp]
    \begin{subfigure}[b]{0.3\textwidth}
     \centering
     \includegraphics[width=\textwidth]{traj_distribution.png}
     \caption{Trajectory Distribution}
     \label{fig:traj_dist}
    \end{subfigure}
   \hfill
    \begin{subfigure}[b]{0.3\textwidth}
     \centering
     \includegraphics[width=\textwidth]{rmsd.png}
     \caption{RMSD from initial (longest sim)}
     \label{fig:rmsd}
    \end{subfigure}
   \hfill
    \begin{subfigure}[b]{0.3\textwidth}
     \centering
     \includegraphics[width=\textwidth]{dssp.png}
     \caption{Secondary Structure (longest sim)}
     \label{fig:sec_struct}
    \end{subfigure}
   \caption{There are %(sim_num)s simulations, totalling %(sim_tot)s us. The longest sim is %(time)s us.}
   \label{fig:initial_anals}
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
   parser.add_argument('-j', '--project', default='%s' %project )
   parser.add_argument('-a', '--sim_num', default='%s' %sim_num )
   parser.add_argument('-b', '--sim_tot', default='%s' %sim_tot )
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
