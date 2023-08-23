#import necessary libraries 
import pandas as pd
from chembl_webresource_client.new_client import new_client
from IPython.display import display

#form data frame of collected proteins
target = new_client.target
target_query = target.search("coronavirus")
targets = pd.DataFrame.from_dict(target_query)
# print(targets.columns)
# print(targets)

#select target protein
selected_target = targets.target_chembl_id[6]
print(selected_target)

#build data frame of target proteins, filtered for IC50
activity = new_client.activity
res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
df = pd.DataFrame.from_dict(res)
# print(df)
display(df)

####visualization of data set
import webbrowser
from tempfile import NamedTemporaryFile
f = open("dataframeview.html",mode="w+")
browser_table = df.to_html()
f.write(browser_table)
webbrowser.open_new("f.html")
