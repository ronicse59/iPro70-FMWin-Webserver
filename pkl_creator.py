from __future__ import division
import methods
import random

# By pass warnings
#====================================================================
import warnings
def warn(*args, **kwargs): pass
warnings.warn = warn

# Define Important Library
#====================================================================
import pandas as pd

# scikit-learn library import
#====================================================================
from sklearn.linear_model import LogisticRegression
from sklearn.externals import joblib

#====================================================================
# Common error to set file path
# http://help.pythonanywhere.com/pages/NoSuchFileOrDirectory/
# There are two approach to setting csv file path
# One is to use absolute path: in our case: /home/ipro70/mysite
'''
import urllib2
url = 'https://raw.githubusercontent.com/ronicse59/iPromoter70/master/selected_feature.csv'
file_path = urllib2.urlopen(url)
#D = pd.read_csv(file_path, header=None)
'''

#====================================================================
# First time add belows comment out code for make pkl file by your model fitted data

file_path = '27_features.csv'
D = pd.read_csv(file_path, header=None)

# Divide features (X) and classes (y) :
#====================================================================
X = D.iloc[:, :-1].values
y = D.iloc[:, -1].values

# Encoding y :
#====================================================================
from sklearn.preprocessing import LabelEncoder
y = LabelEncoder().fit_transform(y)


# Define classifiers within a list
#====================================================================

model = LogisticRegression(n_jobs=1000)

# Save to file in the current working directory
joblib_file = "model.pkl"
joblib.dump(model.fit(X, y), joblib_file)