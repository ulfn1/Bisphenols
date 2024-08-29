# Installing:
conda install python=3.6.6 numpy=1.19.2 scikit-learn=0.23.2 cloudpickle=1.6.0 pandas=1.1.5 rdkit=2020.09.1.0 scikit-learn=0.23.2
pip install nonconformist==1.2.5

# Data download
Download assaydata from Pubchem, e.g. AID_xxx_datatable.csv
Create a tab separated file with cid or sid identifiers followed by SMILES (<smiles_file>)

# SMILES standardization
python rdkit_standardizer3.py <smiles_file> t s n

# Merging assay data with SMILES
./pubchem_read_assay_csv_cid_cutoff.pl AID_xxx_datatable.csv <smiles_file>.rdkit_std.smi 0.6667 c (or s)

# Calculating RdKit descriptors
python rdkit_pc_smi_col1_3.py AID_xxx_datatable.csv.data

# Running model building:
python conformal_prediction.py -n 10 -s t -m t -i AID_xxx_datatable.csv.data.rdkit.txt
# Running model prediction:
python conformal_prediction.py -n 10 -s t -m p -i AID_xxx_datatable.csv.data.rdkit.txt -p <predfile>

# or
# Running model building and prediction:
python conformal_prediction.py -n 10 -s t -m b -i AID_xxx_datatable.csv.data.rdkit.txt -p <predfile>

# Median merging conformal prediction p-values from several models for each class
python cp_median.py <predfile>_nonconf_pred10sum.csv -1



format of AID_xxx_datatable.csv.data.rdkit.txt, <predfile>:
header row: id<sep>target<sep>feature1<sep>feture2 ... featureX
followed by rows for each example (compound)
