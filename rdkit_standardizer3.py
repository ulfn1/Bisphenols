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

import os,sys
import rdkit
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Descriptors

rdkit.RDLogger.DisableLog('rdApp.info')

try:
    sys.argv[1]
except IndexError:
    print ("You need to specify an input tab or space sep smiles file (smiles and id)")
    sys.exit(1)

try:
    sys.argv[2]
except IndexError:
    print ("You need to specify tab (t) or space (s) separated")
    sys.exit(1)

try:
    sys.argv[3]
except IndexError:
    print ("You need to specify 'smiles name' (s) , 'name smiles' (n) , 'name data smiles' (d) order in the file or smiles id target (d2)")
    sys.exit(1)

try:
    sys.argv[4]
except IndexError:
    print ("You need to specify 'remomve salts cmpds' (r) or not (n)")
    sys.exit(1)

infile_path = sys.argv[1]

outfile_path = sys.argv[1]+".rdkit_std.smi"
if os.path.isfile(outfile_path):
    os.remove(outfile_path)

outfile = open(outfile_path, 'w')

nDone=0
mol_weightcutoff = 2000.0

# Process the input
with open(infile_path, 'r') as f:
    if  sys.argv[3] == 'd':
        line = f.readline()
        outfile.write(line.replace(' ', ''))
    for line in f:
        line = line.strip()
        if  sys.argv[2] == 't':
            try:
                if  sys.argv[3] == 's':
                    smiles, nameid = line.split("\t")
                if  sys.argv[3] == 'n':
                    nameid, smiles = line.split("\t")
                if  sys.argv[3] == 'd':
                    try:
                        nameid, data, smiles = line.split("\t")
                    except ValueError:
                        try:
                            nameid, data, smiles, xtra = line.split("\t")
                        except ValueError:
                            nameid, data, smiles, xtra, xtr2, xtra3 = line.split("\t")
                if  sys.argv[3] == 'd2':
                    try:
                        smiles, nameid, data  = line.split("\t")
                    except ValueError:
                        try:
                            smiles, nameid, data, xtra = line.split("\t")
                        except ValueError:
                            smiles, nameid, data, xtra, xtr2, xtra3 = line.split("\t")
            except ValueError:
                print ("ValueError smiles, nameid")
                smiles = 'error_smiles'
        if  sys.argv[2] == 's':
            try:
                if  sys.argv[3] == 's':
                    smiles, nameid = line.split(" ")
                if  sys.argv[3] == 'n':
                    nameid, smiles = line.split(" ")
                if  sys.argv[3] == 'd':
                    try:
                        nameid, data, smiles = line.split(" ")
                    except ValueError:
                        try:
                            nameid, data, smiles, xtra = line.split(" ")
                        except ValueError:
                            nameid, data, smiles, xtra, xtr2, xtra3 = line.split(" ")
                if  sys.argv[3] == 'd2':
                    try:
                        smiles, nameid, data  = line.split(" ")
                    except ValueError:
                        try:
                            smiles, nameid, data, xtra = line.split(" ")
                        except ValueError:
                            smiles, nameid, data, xtra, xtr2, xtra3 = line.split(" ")
            except ValueError:
                print ("ValueError smiles, nameid")
                smiles = 'error_smiles'

        nopass = 0
        if  sys.argv[4] == 'r':
            if "." in smiles:
                print ("//////////////////////////////// removing salt cmpd",smiles,"////////////////////////////////")
                nopass = 1
            
        if  nopass == 0:
            try:
                if "." in smiles:
#                    nameid = nameid + '_removed_salts'
                    print ("/////// salt:", nameid)
                smiles = smiles.replace("@", "")
                smiles = smiles.replace("/", "")
                smiles = smiles.replace("\\", "")
                m = Chem.MolFromSmiles(smiles)
                mol_weight = Descriptors.MolWt(m)
                natoms = m.GetNumHeavyAtoms()
                if natoms > 1:
                    if mol_weight < mol_weightcutoff:
                        new_sml2 = Chem.MolToSmiles(m)
                        print(new_sml2)
                        std = rdMolStandardize.FragmentParent(m)
                        u = rdMolStandardize.Uncharger()
                        mol = u.uncharge(std)
                        mol = Chem.RemoveHs(mol)
                        new_sml  = Chem.MolToSmiles(mol)
                        new_sml = str(new_sml)
                        if  sys.argv[3] == 's':
                            new_line = "%s\t%s\n" % (new_sml, nameid)
                        if  sys.argv[3] == 'n':
                            new_line = "%s\t%s\n" % (nameid, new_sml)
                        if  sys.argv[3] == 'd':
                            new_line = "%s\t%s\t%s\n" % (nameid, data, new_sml)
                        if  sys.argv[3] == 'd2':
                            new_line = "%s\t%s\t%s\n" % (new_sml, nameid, data)
                        outfile.write(new_line)
                        outfile.flush()
                    else:
                        print ("********************", nameid, mol_weight, ">= MW",mol_weightcutoff, "skipping ********************")
                else:
                    print ("********************", nameid, natoms, "< than 2 heavy atoms skipping ********************")
		
                nDone += 1
                if not nDone%1000: print("Smiles done %d"%nDone)

            except:
                print ("error 2", smiles)

outfile.close()

