import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import train_test_split

""" 
Pandas here is a bit like Xarray, i.e., it let's us read in data
and maninipulte it aithe statements like data['XYZ'] etc.
"""

file='../Datasets/housing/housing.csv'

housing = pd.read_csv( file )

print("got housing data")

""" 
Set up new income category 'income_cat' as a new key of 'housing'.
See, just like xarray dataset.
Note: median_income is reported as a numerical value [0,14] i.e. 
units of 10K$. Hence the odd bin edges for 'bins'.
"""
housing["income_cat"] = pd.cut( housing["median_income"],
                                bins=[0., 1.5, 3.0, 4.5, 6.0, np.inf],
                                labels=[1, 2, 3, 4, 5] )



""" Naive random split """
train_set_1, test_set_1 =  train_test_split( housing, test_size=0.2, random_state=42 )


""" better splitting that recognizes need to preserve Income stats """
split = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42 )

for train_index, test_index in split.split( housing, housing["income_cat"]):
    strat_train_set = housing.loc[train_index]
    strat_test_set  = housing.loc[test_index]

""" go on """

bins=np.arange(6)*1.0+.5

d0=housing["income_cat"]
hd0=np.histogram(d0,bins=bins)
print(hd0[0]/len(d0))

d1=test_set_1["income_cat"]
hd1=np.histogram(d1,bins=bins)
print(hd1[0]/len(d1))

d2=strat_test_set["income_cat"]
hd2=np.histogram(d2,bins=bins)
print(hd2[0]/len(d2))
