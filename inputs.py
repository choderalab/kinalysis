
# This dictionary defines which project you are running by your protein name.
# This will have to be hardcorded.

projects = dict()

projects['ABL'] = 11400
projects['SRC'] = 11401
projects['CK2'] = 11406
projects['SYK'] = 11407

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

# Here we store these dictionaries as JSON files.

import json
with open('projects.json', 'w') as fp:
    json.dump(projects, fp)
with open('DFG.json', 'w') as fp:
    json.dump(KER_hbond , fp)
with open('KER_hbond .json', 'w') as fp:
    json.dump(projects, fp)
with open('Aloop_def.json', 'w') as fp:
    json.dump(Aloop_def, fp)
