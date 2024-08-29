"""
MIT License

Copyright (c) 2022 Ulf Norinder

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

#!/usr/bin/env python

from rdkit import Chem
from rdkit.RDLogger import logger
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator
from rdkit.Chem import AllChem as Chem
import sys
logger=logger()
#import cPickle


try:
    sys.argv[1]
except IndexError:
    print ("You need to specify an input file")
    sys.exit(1)

sep = 's'
sep2 = '\ '
try:
    sys.argv[2]
    sep = sys.argv[2]
except IndexError:
    sep = 't'
    sep2 = '\t'


calc = MolecularDescriptorCalculator(['Chi0' , 'Chi0n' , 'Chi0v' , 'Chi1' , 'Chi1n' , 'Chi1v' , 'Chi2n' , 'Chi2v' , 'Chi3n' , 'Chi3v' , 'Chi4n' , 'Chi4v' , 'EState_VSA1' , 'EState_VSA10' , 'EState_VSA11' , 'EState_VSA2' , 'EState_VSA3' , 'EState_VSA4' , 'EState_VSA5' , 'EState_VSA6' , 'EState_VSA7' , 'EState_VSA8' , 'EState_VSA9' , 'FractionCSP3' , 'HallKierAlpha' , 'HeavyAtomCount' , 'Kappa1' , 'Kappa2' , 'Kappa3' , 'LabuteASA' , 'MolLogP' , 'MolMR' , 'MolWt' , 'NHOHCount' , 'NOCount' , 'NumAliphaticCarbocycles' , 'NumAliphaticHeterocycles' , 'NumAliphaticRings' , 'NumAromaticCarbocycles' , 'NumAromaticHeterocycles' , 'NumAromaticRings' , 'NumHAcceptors' , 'NumHDonors' , 'NumHeteroatoms' , 'NumRotatableBonds' , 'NumSaturatedCarbocycles' , 'NumSaturatedHeterocycles' , 'NumSaturatedRings' , 'PEOE_VSA1' , 'PEOE_VSA10' , 'PEOE_VSA11' , 'PEOE_VSA12' , 'PEOE_VSA13' , 'PEOE_VSA14' , 'PEOE_VSA2' , 'PEOE_VSA3' , 'PEOE_VSA4' , 'PEOE_VSA5' , 'PEOE_VSA6' , 'PEOE_VSA7' , 'PEOE_VSA8' , 'PEOE_VSA9' , 'RingCount' , 'SMR_VSA1' , 'SMR_VSA10' , 'SMR_VSA2' , 'SMR_VSA3' , 'SMR_VSA4' , 'SMR_VSA5' , 'SMR_VSA6' , 'SMR_VSA7' , 'SMR_VSA8' , 'SMR_VSA9' , 'SlogP_VSA1' , 'SlogP_VSA10' , 'SlogP_VSA11' , 'SlogP_VSA12' , 'SlogP_VSA2' , 'SlogP_VSA3' , 'SlogP_VSA4' , 'SlogP_VSA5' , 'SlogP_VSA6' , 'SlogP_VSA7' , 'SlogP_VSA8' , 'SlogP_VSA9' , 'TPSA' , 'VSA_EState1' , 'VSA_EState10' , 'VSA_EState2' , 'VSA_EState3' , 'VSA_EState4' , 'VSA_EState5' , 'VSA_EState6' , 'VSA_EState7' , 'VSA_EState8' , 'VSA_EState9' ])

if sep == 't':
    suppl = Chem.SmilesMolSupplier(sys.argv[1],delimiter="\t",smilesColumn = 2,nameColumn=0,titleLine=False)
if sep == 's':
    suppl = Chem.SmilesMolSupplier(sys.argv[1],delimiter="\ ",smilesColumn = 2,nameColumn=0,titleLine=False)

#bb =sys.argv[1]+".rdkit"
bb =sys.argv[1]+".rdkit.txt"
#w = Chem.SmilesWriter(bb)
#w.SetProps(nms)
f2 = open(bb,'w')

nms = list(calc.GetDescriptorNames())
nms2 = "\t".join(str(x) for x in nms)
nms2 = 'name\ttarget\t' + nms2 + '\n'
f2.write(nms2)

nDone=0
f = open(sys.argv[1],'r')
for mol in suppl:
    line = f.readline()
    line = line.strip()
    ID, label, smi  = line.split(sep2)
    nDone += 1
    line = line.strip()
    if not nDone%1000: logger.info("Done %d"%nDone)
    try:
        if mol.GetNumAtoms():
            descrs = calc.CalcDescriptors(mol)
            descrs2 = "\t".join(str(x) for x in descrs)
            descrs2 = str(ID) + '\t' + str(label) + '\t' + descrs2 + '\n'
            f2.write(descrs2)
    except AttributeError:
        print ("Error",mol)

