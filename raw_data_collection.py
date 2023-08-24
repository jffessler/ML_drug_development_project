#import necessary libraries 
import pandas as pd
from chembl_webresource_client.new_client import new_client
from IPython.display import display

#form data frame of collected proteins
target = new_client.target
target_query = target.search("coronavirus")
targets = pd.DataFrame.from_dict(target_query)

#select target protein
selected_target = targets.target_chembl_id[6]

#build data frame of target proteins, filtered for IC50
activity = new_client.activity
res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
df = pd.DataFrame.from_dict(res)
# display(df.head(3))
# print(df.standard_type.unique())

####view data set
## import webbrowser
## from tempfile import NamedTemporaryFile
## f = open("dataframeview.html",mode="w+")
## browser_table = df.to_html()
## f.write(browser_table)
## webbrowser.open_new("f.html")

##note: "Standard Value" is the potency of the drug, higher number lower potency, we are looking for a lower standard value

#create .csv file of target bioactivity
# df.to_csv("bioactivity_data.csv", index=False)
data = pd.read_csv("bioactivity_data.csv")

#drop data with missing Standard Value
df2 = df[df.standard_value.notna()]

#define active compounds with IC50 < 1000nM
bioactivity_class = []
for i in df2.standard_value:
    if float(i) >= 10000:
        bioactivity_class.append("inactive")
    elif float(i) <= 1000:
        bioactivity_class.append("active")
    else:
        bioactivity_class.append("intermediate")
# print(bioactivity_class)

#identify need columns
selection = ['molecule_chembl_id','canonical_smiles', 'standard_value']
df3 = df2[selection]

#combine data sets
df4 = pd.concat([df3,pd.Series(bioactivity_class,name='bioactivity_class')], axis=1)

#save new refined set of data to .csv
# df4.to_csv("bioactivity_preprocessed_data.csv", index=False)


