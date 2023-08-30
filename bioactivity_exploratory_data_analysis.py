import sys
import pandas as pd
from IPython.display import display
import numpy as np
from rdkit import Chem
from rdkit.Chem import Lipinski, Descriptors
import seaborn as sns
import matplotlib.pyplot as plt
import shutil

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
# display(df_lipinski)
# display(df)
df_combined = pd.concat([df,df_lipinski], axis=1)
# display(df_combined)

def norm(input):
    norm = []

    for i in input["standard_value"]:
        if i > 100000000:
            i = 100000000
        norm.append(i)
    
    input["standard_value_norm"] = norm
    x = input.drop(columns = "standard_value")
    return x

def pIC50(input):
    pIC50 = []

    for i in input["standard_value_norm"]:
        molar = i*(10**-9) #convert nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop(columns='standard_value_norm')

    return x

# print(df_combined.standard_value.describe())

df_norm = norm(df_combined)
# display(df_norm)
# print(df_combined.standard_value_norm.describe())

df_final = pIC50(df_norm)
# display(df_final)
# print(df_final.pIC50.describe())

df_2class = df_final[df_final.bioactivity_class != "intermediate"]
# print(df_2class.pIC50.describe())
display(df_2class)

#### Exploratory Data Analysis ####
#### Chemical Space Analysis #### <- correct terminology for drug testing
sns.set(style='ticks')

###plot 1: bar plot of bioactivity
# plt.figure(figsize=(5.5,5.5))
# sns.countplot(x="bioactivity_class", data=df_2class, edgecolor="black")

# plt.xlabel("Bioactivity Class", fontsize=14,fontweight="bold")
# plt.ylabel("Frequency", fontsize=14,fontweight="bold")

# plt.savefig('plot_bioactivity_class.pdf')
# shutil.move("/Users/johnfessler/Desktop/Coding Practice/ML_drug_development_project/plot_bioactivity_class.pdf", "/Users/johnfessler/Desktop/Coding Practice/ML_drug_development_project/plots/plot_bioactivity_class.pdf")


### plot 2: scatter plot of molecular weight and logP
# plt.figure(figsize=(5.5, 5.5))
# ax = plt.subplot(111)

# sns.scatterplot(x='MW', y='LogP', data=df_2class, hue='bioactivity_class', size='pIC50', edgecolor='black', alpha=0.7)

# plt.xlabel('MW', fontsize=14, fontweight='bold')
# plt.ylabel('LogP', fontsize=14, fontweight='bold')
# # plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0)

# # Shrink current axis by 20%
# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# # Put a legend to the right of the current axis
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# plt.savefig('plot_MW_vs_LogP.pdf')
# shutil.move("/Users/johnfessler/Desktop/Coding Practice/ML_drug_development_project/plot_MW_vs_LogP.pdf", "/Users/johnfessler/Desktop/Coding Practice/ML_drug_development_project/plots/plot_MW_vs_LogP.pdf")

###plot 3: box plot of pIC50 (bioactivity)
# plt.figure(figsize=(5.5,5.5))
# sns.boxplot(x="bioactivity_class", y="pIC50", data=df_2class)

# plt.xlabel("Bioactivity Class", fontsize=14,fontweight="bold")
# plt.ylabel("pIC50", fontsize=14,fontweight="bold")

# plt.savefig('plot_pIC50.pdf')
# shutil.move("/Users/johnfessler/Desktop/Coding Practice/ML_drug_development_project/plot_pIC50.pdf", "/Users/johnfessler/Desktop/Coding Practice/ML_drug_development_project/plots/plot_pIC50.pdf")

###  Mann Whitney Test
def mannwhitney(descriptor,df_2class):
  # https://machinelearningmastery.com/nonparametric-statistical-significance-tests-in-python/
  from numpy.random import seed
  from numpy.random import randn
  from scipy.stats import mannwhitneyu

# seed the random number generator
  seed(1)

# actives and inactives
  selection = [descriptor, 'bioactivity_class']
  df = df_2class[selection]
  active = df[df['bioactivity_class'] == 'active']
  active = active[descriptor]

  selection = [descriptor, 'bioactivity_class']
  df = df_2class[selection]
  inactive = df[df['bioactivity_class'] == 'inactive']
  inactive = inactive[descriptor]

# compare samples
  stat, p = mannwhitneyu(active, inactive)
  #print('Statistics=%.3f, p=%.3f' % (stat, p))

# interpret
  alpha = 0.05
  if p > alpha:
    interpretation = 'Same distribution (fail to reject H0)'
  else:
    interpretation = 'Different distribution (reject H0)'
  
  results = pd.DataFrame({'Descriptor':descriptor,
                          'Statistics':stat,
                          'p':p,
                          'alpha':alpha,
                          'Interpretation':interpretation}, index=[0])
  filename = 'mannwhitneyu_' + descriptor + '.csv'
  results.to_csv(filename)

  return results

print(mannwhitney('pIC50',df_2class))