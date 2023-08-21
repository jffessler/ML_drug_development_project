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

#
activity = new_client.activity
res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
df = pd.DataFrame.from_dict(res)
# print(df)
display(df)

#visualization
# import webbrowser
# from tempfile import NamedTemporaryFile
# with NamedTemporaryFile(mode="w+b", delete=False, suffix='.html') as f:
#     print("1")
#     df.to_html(f)
#     print("2")
#     webbrowser.open(f.name)
#     print("3")
