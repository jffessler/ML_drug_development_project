import sys
import pandas as pd
from IPython.display import display
import numpy as np
from rdkit import Chem
from rdkit.Chem import Lipinski, Descriptors

#import data for analysis
df = pd.read_csv("bioactivity_preprocessed_data.csv")
# display(df)

## Lipinski descriptors
## molecular weight < 500 Dalton
## Octanol-water partition energy (LogP) < 5
## Hydrogen bond donors < 5
## Hydrogen bond acceptors <10

def lipinski(smiles, verbose = False):
    moldata = []

    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1,1)
    i = 0

    for mol in moldata:
        desc_MolWT = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)

        row = np.array([desc_MolWT, desc_MolLogP, desc_NumHDonors, desc_NumHAcceptors])

        if (i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData,row])
        i += 1
    
    columnNames = ["MW", "LogP", "NumHDonors", "NumHAcceptors"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)
    return descriptors

df_lipinski = lipinski(df.canonical_smiles)
display(df_lipinski)