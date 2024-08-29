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

#Script for the p-value calculation

# imports
import numpy as np
import sys
import pandas as pd


################## Main ##################

    
try:
    sys.argv[1]
except IndexError:
    print ("You need to specify an input tab sep file with p-values")
    sys.exit(1)

try:
    sys.argv[2]
except IndexError:
    print ("You need to specify a number of models to be used (<0 for all)")
    sys.exit(1)

inpath = sys.argv[1]

compound_dict_id = {}
compound_dict_0 = {}
compound_dict_1 = {}
notificationstatus = 0
notificationstatus1 = 0

with open(inpath) as f:
    header = next(f)
    idx =  0
    for line in f:
        line = line.strip()
        readstatus = 0
        try:
            ID, pos0, pos1, label, model  = line.split("\t")
            pos0 = float(pos0)
            if pos0 == 1 or pos0 <= 0:
                readstatus = 1
        except:
            readstatus = 1
        if readstatus == 1:
            readstatus = 0
            if notificationstatus < 5:
                print ("Trying median p-value format instead")
                notificationstatus = notificationstatus + 1
            try:
                ID, label, pos0, pos1, dummy  = line.split("\t")
                model = 0
            except:
                readstatus = 1
        if readstatus == 1:
            if notificationstatus1 < 5:
                print ("Trying median p-value format instead with additional column")
                notificationstatus1 = notificationstatus1 + 1
            try:
                ID, label, pos0, pos1, predcl, dummy  = line.split("\t")
                model = 0
            except:
                sys.exit(1)

        try:
            pos0 = float(pos0)
        except:
            print("Error reading line:", line)
            line = next(f)
            line = line.strip()
            ID, pos0, pos1, label, model  = line.split("\t")
            
        repeat = 1
        idx =  idx + 1
        if idx % 1000000 == 0 and idx > 0: print ("done:",idx)

        try:
            compound_dict_id[ID]
        except:
            compound_dict_id[ID] = idx

        label = float(label)
        pos0 = float(pos0)
        pos1 = float(pos1)
        model = int(model)

        if model > int(sys.argv[2]) and int(sys.argv[2]) > 0:
            model = int(sys.argv[2])
            print ("Stopping after model",model)
            break

        if int(label) <= 0:
            try:
                compound_dict_0[ID].append([pos0, pos1])
            except:
                compound_dict_0[ID] = [[pos0, pos1]]
        else:
            try:
                compound_dict_1[ID].append([pos0, pos1])
            except:
                compound_dict_1[ID] = [[pos0, pos1]]

f.close()

# Calculate the performance at the set sig_lvls

inpath = sys.argv[1] + '.conf_pred_results_median_class_models_'  + str(model)

compound_tot = []

for key in compound_dict_1:    
    aa =  key  + '\t' + str(np.median(compound_dict_1[key], axis=0)[0]) + '\t' + str(np.median(compound_dict_1[key], axis=0)[1]) + '\t' + str(compound_dict_id[key]) + '\n'
    compound_tot.append(aa)

for key in compound_dict_0:
    aa =  key + '\t' + str(np.median(compound_dict_0[key], axis=0)[0]) + '\t' + str(np.median(compound_dict_0[key], axis=0)[1]) + '\t' + str(compound_dict_id[key]) # + '\n'
    compound_tot.append(aa)

df = pd.DataFrame({"Col": compound_tot})
df = df.pop('Col').str.split('\t',expand=True)
cols = ['id','cp_p-value_class_0','cp_p-value_class_1','order']
df.columns = cols
df = df.astype({"order": int})
df = df.sort_values(by=['order'], ascending=True)
df = df.drop(['order'], axis=1, errors='ignore')
df.to_csv(inpath, sep = '\t', header=True, index=False) 
