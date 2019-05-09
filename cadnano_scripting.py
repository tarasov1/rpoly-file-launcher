import sys
import math
import re 
import cadnano
from cadnano.document import Document
from cadnano.part.nucleicacidpart import NucleicAcidPart
from cadnano.part.createvhelixcmd import CreateVirtualHelixCommand
from cadnano.strand import Strand
from cadnano.strandset import StrandSet
from cadnano.strandset import CreateStrandCommand
from cadnano.cntypes import (
    StrandSetT
)
from cadnano.proxies.cnenum import (
	GridEnum
)
from pyquaternion import Quaternion
import numpy as np
#import matplotlib as mpl
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt

#Read File 
if len(sys.argv) == 1:
    print('Need filename')
    sys.exit(-1)

#'data' stores helix coordinates 
data = []
# fwd_helix_connections and rev_helix_connections: number rows is ammount of helices, 
# every row stores an oligo connection 
fwd_helix_connections = []
rev_helix_connections = []
count = 0

polyFile = open(sys.argv[1], 'r')

try: 
	for line in polyFile: 
		if line.startswith('hb'):
			data.insert(count, line.split(' '))
			count += 1 
		elif line.startswith('c'):
			if 'f3' not in line:  
				#rev_helix_connections.append(line.split(' '))
				rev_helix_connections.append([int(re.search('c helix_(.+?) ',
					line).group(1)),int(re.search('\' helix_(.+?) ',line).group(1))])
			else: 
				fwd_helix_connections.append([int(re.search('c helix_(.+?) ',
					line).group(1)),int(re.search('\' helix_(.+?) ',line).group(1))])
except Exception: 
	print('Failed to read the file')

print('The total number of helices is', len(data)) 
 
#Cadnano procedure 
app = cadnano.app()
doc = app.document = Document() 
listOfHelices = [] 

# Should be grid_type=GridEnum.NONE but cadnano crashes 
part = doc.createNucleicAcidPart(grid_type=GridEnum.SQUARE, is_lattice=False)

#part = doc.createNucleicAcidPart()
oligo_pos = []
vhCounter = 0 
length = 0

# Reads orientation from the "data" and produces rotations from the Quaternian coordinates 
for n,i in enumerate(data): 	
	q = Quaternion(i[6:10]) 
	vec = float(i[2])*q.rotate(np.array([0.0,0.0,1.0]))/6.0 
	part.createVirtualHelix(float(i[3]), float(i[4]), float(i[5]), 
		length=int(i[2]), direction=(vec[0],vec[1],vec[2]), id_num=n)
	#listOfHelices.append(part)
	vhCounter += 1 
	#Determine position for oligo strands (first oligo: from 1 --> length/2 second oligo: next (length/2)+1 --> length )
	oligo_pos.append([1,int(int(i[2])/2),int(int(i[2])/2)+1,int(i[2])])
	
# Max length of the Strand Set
MAX_LENGTH = 64
offset_idx = 0

# Make the container for the StrandSets 
ss_fwd = [None]*vhCounter
ss_rev = [None]*vhCounter

# Create the StrandsSets and fill them with Strands 
# The assumption here is that Scaffold is forward strand and that the reverse strands is Staple Strands
try: 
	for i in range(vhCounter):	
		ss_fwd[i] = StrandSet(is_fwd=True, id_num=i, part=part, initial_size=MAX_LENGTH)
		ss_rev[i] = StrandSet(is_fwd=False, id_num=i, part=part, initial_size=MAX_LENGTH)   
		ss_fwd[i].createStrand(base_idx_low=offset_idx, base_idx_high=offset_idx+oligo_pos[i][3], color='#0066cc')
		ss_rev[i].createStrand(base_idx_low=offset_idx, base_idx_high=offset_idx+oligo_pos[i][1], color='#f44242')
		ss_rev[i].createStrand(base_idx_low=offset_idx+oligo_pos[i][2], base_idx_high=offset_idx+oligo_pos[i][3], color='#007200')
		part.fwd_strandsets[i] = ss_fwd[i]
		part.rev_strandsets[i] = ss_rev[i] 
except Exception: 
	print('Failed to create StrandSet/Strands')

"""
print('\t', fwd_ss, '\t', [s.idxs() for s in fwd_ss.strands()], '\n\t\t\t\t',
          [s.getColor() for s in fwd_ss.strands()])
print('\t', rev_ss, '\t', [s.idxs() for s in rev_ss.strands()], '\n\t\t\t\t',
          [s.getColor() for s in rev_ss.strands()])

for oligo in oligos: 
	print("{0}\t{1}\t\'{2}\'\t{3}".format(oligo,oligo.length(),oligo.getCOlor(),oligo.sequence()))
"""

# Creation of specified crossovers from the rpoly-file for the Scaffold  
try:  
	for i in range(vhCounter): 
		if i==vhCounter-1:
			part.createXover(strand5p=part.fwd_strandsets[i].getStrand(base_idx=offset_idx+oligo_pos[i][3]), idx5p=offset_idx+oligo_pos[i][3],
				 strand3p=part.fwd_strandsets[0].getStrand(base_idx=offset_idx), idx3p=offset_idx,
				 update_oligo=True, allow_reordering=True, use_undostack=True
				)
		else: 
			part.createXover(strand5p=part.fwd_strandsets[i].getStrand(base_idx=offset_idx+oligo_pos[i][3]), idx5p=offset_idx+oligo_pos[i][3],
				 strand3p=part.fwd_strandsets[i+1].getStrand(base_idx=offset_idx), idx3p=offset_idx,
				 update_oligo=True, allow_reordering=True, use_undostack=True
				)	
except Exception: 
	print('Failed to connect Forward Oligos')

# Creation of the crossovers for the Staple Strands 
# rev_helix_connections has information about the connections 
try: 
	for i in range(vhCounter): 
		part.createXover(strand5p=part.rev_strandsets[rev_helix_connections[i][0]-1].getStrand(base_idx=offset_idx), idx5p=offset_idx,
				 strand3p=part.rev_strandsets[rev_helix_connections[i][1]-1].getStrand(base_idx=offset_idx+oligo_pos[i][3]), idx3p=offset_idx+oligo_pos[i][3],
				 update_oligo=True, allow_reordering=True, use_undostack=True
				)
except Exception: 
	print('Failed to connect Reverse Oligos')

#Encoding
try: 
	doc.writeToFile(sys.argv[1]+'.json')
except Exception: 
	print('Failed to write to file')