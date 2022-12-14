import os
import tarfile
import urllib.request



def fetch_housing_data(housing_url,housing_path):
    os.makedirs( housing_path , exist_ok=True)
    tgz_path = os.path.join( housing_path, "housing.tgz" )
    urllib.request.urlretrieve( housing_url, tgz_path )
    housing_tgz =  tarfile.open( tgz_path )
    housing_tgz.extractall(path=housing_path)
    housing_tgz.close()


DOWNLOAD_ROOT="https://raw.githubusercontent.com/ageron/handson-ml2/master/"
HOUSING_PATH=os.path.join( "../Datasets_2","housing")
HOUSING_URL=DOWNLOAD_ROOT+"datasets/housing/housing.tgz"


fetch_housing_data(housing_url=HOUSING_URL  ,  housing_path=HOUSING_PATH)



