import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import train_test_split
from sklearn.impute import SimpleImputer

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

""" better splitting that recognizes need to preserve Income stats """
split = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42 )

for train_index, test_index in split.split( housing, housing["income_cat"]):
    strat_train_set = housing.loc[train_index]
    strat_test_set  = housing.loc[test_index]

"""               Start to clean data            """

""" 
First 'drop' learning target from training set.
Not really sure what these lines do
"""    
housing = strat_train_set.drop("median_house_value", axis=1 )
housing_labels = strat_train_set["median_house_value"].copy()


housing_num = housing.drop("ocean_proximity", axis=1 )

imputer = SimpleImputer( strategy = "median" )

imputer.fit( housing_num )

X=imputer.transform( housing_num )  # X is an np ndarray

""" 
Check out effect of 'imputer' on total_bedrooms.  
Pandas isnull() will be True when data is NaN inf etc
"""
brs=housing['total_bedrooms']
noo=brs.isnull() # noo here is a 'Pandas series'. Could be converted to np ndarray. Doesn't simplify though.
oo=np.where( noo==True )[0] # The [0] makes this a comprehensible vector of indices instaed of some weird Python tuple 

print("\nFirst 11 values of original training data in total_bedrooms. Selected with isnull() ")
print(brs.values[oo[0:10]])

print("\nFirst 11 values in relevant column of np ndarray X")
print(X [ oo[0:10], 4 ] )

print( "\n shows that values were replaced with \n\n    brs.median()=",  brs.median() ,"\n\n by scikit imputer" )




