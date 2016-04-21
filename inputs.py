
# This dictionary defines which project you are running by your protein name.
# This will have to be hardcorded.

projects = dict()

projects['ABL'] = 11400
projects['SRC'] = 11401
projects['DDR1'] = 11402
projects['DDR1-pro'] = 11403
projects['CK2'] = 11406
projects['SYK'] = 11407

projects['SRC-ss-NVT'] = 10467
projects['ABL-ss-NVT'] = 10468
projects['EGFR-ss-NVT'] = 10469

projects['SRC-ss'] = 10471
projects['ABL-ss'] = 10472
projects['EGFR-ss'] = 10473

projects['DDR1-1s'] = 10484
projects['SYK-3s'] = 10485
projects['CK2-1s'] = 10486
projects['HER2-2s'] = 10487
projects['MEK-4s'] = 10488

#mTOR not included
#AURKA not included

# This dictionary defines which dihedrals define the DFG flip for each protein. Hopefully this will eventually
# be defined programmatically rather than by hand.

DFG = dict()

DFG['ABL'] = [2257,2255,2265,2270]
DFG['SRC'] = [2190,2188,2198,2203]
DFG['CK2'] = [2895,2893,2912,2917]
DFG['SYK'] = [2414,2412,2423,2428]

print 'Projects:'
print projects
print 'DFG dihedral definitions:'
print DFG

# These two dictionary defines our Shukla coordinates

# Define hydrogen bond coordinates (0-indexed)
KER_hbond = { 'SRC' : [[28,43],[43,142]],
              'ABL' : [[29,44],[44,144]],
              'DDR1': [[51,68],[68,185]]}

# Define Activation loop (resid)
Aloop_def = { 'SRC': [138,158],
          'ABL': [140,160],
          'DDR1': [181,201]}

print 'KER hydrogen bond definitions:'
print KER_hbond
print 'Activation loop definitions :'
print Aloop_def

# Here we store these dictionaries as JSON files.

import json
with open('projects.json', 'w') as fp:
    json.dump(projects, fp)
with open('DFG.json', 'w') as fp:
    json.dump(DFG, fp)
with open('KER_hbond.json', 'w') as fp:
    json.dump(KER_hbond, fp)
with open('Aloop_def.json', 'w') as fp:
    json.dump(Aloop_def, fp)
