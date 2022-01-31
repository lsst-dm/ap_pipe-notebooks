#! /usr/bin/env python

from reg_tools import list_datasets
import lsst.daf.butler as dafButler
import sys
import pandas as pd

repo = sys.argv[1]
collection = sys.argv[2]

butler = dafButler.Butler(repo)

rr = list_datasets(butler,
                   collections=collection,
                   dataID='',#{'instrument': 'LSSTCam-imSim','skymap':'DC2'},
                   datasetType='*',
                   print_type=True,
                   max_print=0,
                   #where=''
                   );

print(f'Total items in {collection}: {len(rr)}')
print()

#-- show some summary stats of the contents
# note that we may have a hetrogeneous dataframe composed of various types of rows!!
ids = []
for r in rr:
    ids.append(r.dataId.byName())
    ids[-1]['type']=r.datasetType.name
    
df = pd.DataFrame(ids)
#df = df.astype({"visit": int, "detector": int},errors='ignore')


#- counts of items of each type
print(df.groupby('type')['instrument'].count())
print()

#- ranges of visits/detectors or tracts/patches
try:
    print(df.groupby('type').
          agg({'visit':['min','max'],
               'detector':['min','max'],
               }))
except:
    pass

try:
    print(df.groupby('type').
          agg({'tract':['min','max'],
               'patch':['min','max'],
               }))
except:
    pass
print()

#- a simple describe() of the whole content
print(df.describe().apply(lambda s: s.apply('{0:.5f}'.format)))
print()
