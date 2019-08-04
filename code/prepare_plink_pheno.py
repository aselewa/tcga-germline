#!/bin/python

import numpy as np
import pandas as pd

# somatic mutation file
ssms = pd.read_csv('SSM_MAFs/stddata__2016_01_28/ALL_SSMs_unfiltered.txt',sep='\t',header=None)

#list of driver genes and their tumor type(s)
drivers = pd.read_csv('driver_genes_ttype.txt',sep='\t',header=None)

#plink family file
fam = pd.read_csv('PLINK/TCGA_European.fam',sep=' ', header=None)
fam = fam.drop([2,3,4,5],axis=1)
fam.index = fam[1]
fam.columns = ['FID','IID']

n = drivers.shape[0]

for i in range(n):
    g = drivers.iloc[i,0]
    ttype = ['TCGA-'+s for s in drivers.iloc[i,1].split(',')]

    fam[g] = 1
    ids = ssms.loc[ssms[0]==g,5]
    samples = np.unique(['-'.join(i.split('-')[0:3]) for i in ids])
    in_samples = np.intersect1d(samples, fam.index.values)
    fam.loc[in_samples, g] = 2

    BOL = np.isin(fam['FID'].values, ttype)    
    fam.loc[~BOL,g] = -9

numInd = (fam==2).sum(axis=0)
numInd[0:2] = np.inf
filtered_fam = fam.loc[:,(numInd > 20)]
filtered_fam.to_csv('PLINK/drivers_pooled_Fisher/SSM_phenotypes_for_plink.txt',sep=' ', index=False, header=True)
print(filtered_fam.shape)