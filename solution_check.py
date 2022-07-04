import numpy as np
import sys

inFile = sys.argv[1] 'chembl_emerging_closed_8_5.csv'
NUM_TRANS = int(sys.argv[2])

# Pattern size check:
print('Pattern size check:')
f = open('results/'+inFile, 'r')

line = f.readline()

lenPattern = []
for line in f:
    if line !='\n':
        data = line.strip().split(',')
        pattern = data[1].strip('"').split(' ')
        lenPattern.append(len(pattern))
        
print('min size: '+str(np.min(lenPattern)))
print('max size: '+str(np.max(lenPattern)))
print('mean size: '+str(np.round(np.mean(lenPattern),1)))
print('median size: '+str(np.round(np.median(lenPattern),1)))
#%%
# Coverage check:
print('\nCoverage check:')

f = open('results/'+inFile, 'r')
    
line = f.readline()

coveredMolecules = []
for line in f:
    if line !='\n':
        data = line.strip().split(',')
        molAct = data[2].strip('"').split(' ')
        if molAct != ['']:
            for mol in molAct:
                if int(mol) not in coveredMolecules:
                    coveredMolecules.append(int(mol))
        molInact = data[3].strip('"').split(' ')
        if molInact != ['']:
            for mol in molInact:
                if int(mol) not in coveredMolecules:
                    coveredMolecules.append(int(mol))
f.close()

if len(coveredMolecules) == NUM_TRANS:
    print('All molecules are covered by the solution')
else:
    print('Only '+str(np.round(len(coveredMolecules)*100/float(NUM_TRANS),2))+'% molecules are covered by the solution\n')
#%%
# Solution coverage checker:
f = open('results/'+inFile, 'r')

line = f.readline()

lenSolution = []
for line in f:
    if line !='\n':
        data = line.strip().split(',')
        molTest = (data[2].strip('"')+' '+data[3].strip('"')).strip(' ')
        solution = molTest.split(' ')
        lenSolution.append(len(solution))
        
print('min support: '+str(np.min(lenSolution)))
print('max support: '+str(np.max(lenSolution)))
print('mean support: '+str(np.round(np.mean(lenSolution),1)))
print('median support: '+str(np.round(np.median(lenSolution),1)))
#%%
