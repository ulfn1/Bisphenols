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
 
def data_loader(infile,sep):
    if sep == 't':
        data = pd.read_csv(infile, sep='\t', header = 0, index_col = None)
    if sep == 'c':
        data = pd.read_csv(infile, sep=',', header = 0, index_col = None)
    if 'Unnamed: 0' or 'name' or 'Name' or 'Molecule name' or 'CASRN' in data.columns:
        data.rename(columns={'Unnamed: 0': 'id', 'name': 'id', 'Name': 'id', 'Molecule name': 'id', 'CASRN': 'id'}, inplace=True)
    if ' Label' or 'class' or 'very_toxic' or 'nontoxic' in data.columns:
        data.rename(columns={' Label': 'target', 'class': 'target', 'very_toxic': 'target', 'nontoxic': 'target'}, inplace=True)
    return data

def writeOutListsApp(filePath, dataLists, delimiter = '\t'):
    with open(filePath, 'a') as fd:
        numOfLists = len(dataLists)
        lenOfLists = len(dataLists[0])
        for lineNum in range(lenOfLists):
            tmpString = ""
            for column in range(numOfLists):
                tmpString += str(dataLists[column][lineNum]) + delimiter
            fd.write(tmpString.strip() + "\n")
        fd.flush()
    fd.close()

import os,sys
from sklearn.ensemble import RandomForestClassifier
from nonconformist.icp import IcpClassifier
from nonconformist.nc import ProbEstClassifierNc, margin
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from sklearn.utils import shuffle
import cloudpickle
from sklearn.impute import SimpleImputer
import random
import gzip

parser = ArgumentParser()
parser.add_argument('-i','--infile', help='input training file')
parser.add_argument('-n','--nmodels', type=str, help='number of models (default 20 models)')
parser.add_argument('-m','--mode', type=str, choices=['t','p', 'b'], help='mode: build models, predict new data from models, both build and predict')
parser.add_argument('-s','--sep', type=str, choices=['t','c'], help='file separator: tab or comma')
parser.add_argument('-p','--predfile', help='input prediction file if mode == p')
args = parser.parse_args()

infile = args.infile
nmodels = args.nmodels
mode = args.mode
sep = args.sep
predfile = args.predfile

if infile == None or mode == None or sep == None:
    parser.print_help(sys.stderr)
    sys.exit(1)

if mode == 'p' and predfile == None or mode == 'b' and predfile == None:
    parser.print_help(sys.stderr)
    sys.exit(1)

if nmodels == None:
    nmodels = 20
nmodels = int(nmodels)

modelfile = infile +"_nonconf"+".model.gz"
if mode != 'p':
    if os.path.isfile(modelfile):
        os.remove(modelfile)
#    f= open(modelfile, mode='ab')
    with gzip.open(modelfile, mode="ab") as f:
        cloudpickle.dump(nmodels, f)
    data = data_loader(infile,sep)
    print (data)

    try:
        labels = data['id'].values
    except:
        data.rename(columns={ data.columns[0]: "id" }, inplace = True)
        labels = data['id'].values
        print("/////////////////////// renaming column 1 to 'id' ///////////////////////")

    try:
        data.loc[data['target'] < 0, 'target'] = 0
    except:
        data.rename(columns={ data.columns[1]: "target" }, inplace = True)
        data.loc[data['target'] < 0, 'target'] = 0
        print("/////////////////////// renaming column 2 to 'target' ///////////////////////")
        print (data)

    labels = data['id']
    labels2 = pd.DataFrame(data.loc[:, 'id'])
    data['idsplit'] = data['id']
    print (labels2)
    idlist = labels2['id'].unique()
    print (idlist)
    data = data.drop(['dataset'], axis=1, errors='ignore')
    data.replace([np.inf, -np.inf], np.nan, inplace=True)
    nancount = data.isnull().sum().sum()
    if nancount > 0:
        print ("imputing missing values:",nancount)
        imp = SimpleImputer(missing_values=np.nan, strategy='mean')
        impmodel = imp.fit(data)
        data = impmodel.transform(data)


    part1 = int(0.7*len(idlist))
    for xx in range(1, nmodels+1):
        modelfile2 = infile +"_nonconf"+"_"+str(xx)+".model.gz"
        print ("Working on model", xx)
 
        SEED = 1234 * xx
        random.seed(SEED)
        random.shuffle(idlist)
        idx = idlist

        trainset = idx[:part1]
        calset = idx[part1:]
        print(trainset)
        print(calset)
        dftrain = data[data['idsplit'].isin(trainset)]
        dfcal = data[data['idsplit'].isin(calset)]
        targettrain = dftrain['target'].values
        train = dftrain.drop(columns=['id','idsplit', 'target'], axis=1, errors='ignore').values
        targetcal = dfcal['target'].values
        calibr = dfcal.drop(columns=['id','idsplit', 'target'], axis=1, errors='ignore').values
        del dftrain, dfcal

        nc = ProbEstClassifierNc(RandomForestClassifier,margin,model_params={'n_estimators': 100})
        icp_norm = IcpClassifier(nc, condition=lambda instance: instance[1])

        icp_norm.fit(train, targettrain)
        icp_norm.calibrate(calibr, targetcal)
        with gzip.open(modelfile, mode="ab") as f:
            cloudpickle.dump(icp_norm, f)
    f.close()

if mode != 't':
    print("Reading in data ...")
    data = data_loader(predfile,sep)

    try:
        labels = data['id'].values
    except:
        data.rename(columns={ data.columns[0]: "id" }, inplace = True)
        labels = data['id'].values
        print("/////////////////////// renaming column 1 to 'id' ///////////////////////")

    try:
        data.loc[data['target'] < 0, 'target'] = 0
    except:
        data.rename(columns={ data.columns[1]: "target" }, inplace = True)
        data.loc[data['target'] < 0, 'target'] = 0
        print("/////////////////////// renaming column 2 to 'target' ///////////////////////")
        print (data)
    labels = data['id'].values
    target = data['target'].values
    cls_uniq = data['target'].nunique()
    print("Number of classes:",cls_uniq)
    test = data.drop(['id'], axis=1, errors='ignore')
    test = test.drop(['dataset'], axis=1, errors='ignore')
    for col in test.columns:
        test[col] = pd.to_numeric(test[col], errors='coerce')
    test.replace([np.inf, -np.inf], np.nan, inplace=True)
    nancount = test.isnull().sum().sum()
    test = test.drop(['target'], axis=1, errors='ignore').values
    del data

    if nancount > 0:
        print ("imputing missing values:",nancount)
        imp = SimpleImputer(missing_values=np.nan, strategy='mean')
        impmodel = imp.fit(test)
        test = impmodel.transform(test)

#    f= open(modelfile, mode='rb')
    with gzip.open(modelfile, mode="rb") as f:
        nmodels_built = cloudpickle.load(f)
        print ("Models built:",  nmodels_built)
        outfile = predfile + "_nonconf_pred" + str(nmodels_built)
        outfile = outfile + "sum.csv"

        f2 = open(outfile,'w')
        if cls_uniq > 2:
            f2.write('id\tp-value_low_class\tp-value_high_class\tp-value_highest_class\tclass\tmodel\n')
        else:
            f2.write('id\tp-value_low_class\tp-value_high_class\tclass\tmodel\n')
        f2.close()
        if nmodels > nmodels_built:
            print ("More models ordered (",nmodels,") than the file contains. Setting the number of models to built models.")
            nmodels = nmodels_built

        for xx in range(1, nmodels+1):
#            print ("\r Predicting from model", xx, end='')
            print ("\r Predicting from model", xx)
            modelfile2 = infile +"_nonconf"+"_"+str(xx)+".model.gz"
            icp_norm = cloudpickle.load(f)
        
            chunk_size = int(50000)
#            chunk_size = int(500)
            chunk_tot = float(len(test)) / chunk_size
            i = 0
            for start in range(0, len(test), chunk_size):
                i = i+ 1
                print ("\r using chunk (",chunk_tot,")", i , end='')
#                print ("\r using chunk (",chunk_tot,")", i)
                test_subset = test[start:start + chunk_size]
                target_subset = target[start:start + chunk_size]
                labels_subset = labels[start:start + chunk_size]
                ll = len(labels_subset)
                num = [xx]*ll
                predicted = icp_norm.predict(test_subset)
                print ("predicted raw output")
                print (predicted)
                predicted0 = [x[0] for x in predicted]
                predicted1 = [x[1] for x in predicted]
                if cls_uniq > 2:
                    predicted2 = [x[2] for x in predicted]
                    writeOutListsApp(outfile, [labels_subset, predicted0,  predicted1,  predicted2, target_subset, num])
                else:
                    writeOutListsApp(outfile, [labels_subset, predicted0,  predicted1, target_subset, num])
        f.close()

print (" - finished\n")
 
