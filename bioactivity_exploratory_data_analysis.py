import sys
import pandas as pd
from IPython.display import display

#import data for analysis
df = pd.read_csv("bioactivity_preprocessed_data.csv")
display(df)