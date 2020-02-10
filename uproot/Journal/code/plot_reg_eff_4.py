
import uproot 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

def reg_eff(trees, lum, label, soft=True):
    f = uproot.pandas.iterate(trees+"/*.root", "t", reportpath = True,
                          branches=['beam_energy', 'mass', 'trigger', 'procedure'])
    dataRAW = pd.concat([j.assign(name=float(re.findall(r"\d+\.\d+", path)[0])) for path,j in f])
    dataRAW = dataRAW.query('procedure!=1').drop(['procedure', 'trigger'], axis=1).rename({'name':'label'}, axis=1)
    dataRAW = dataRAW.groupby('label').agg({'mass': 'count'}).reset_index()
    if soft:
        lum = pd.read_csv(lum).\
                drop_duplicates(subset=['label']).sort_values(by='label')
        dataRAW = dataRAW.merge(lum, on=['label'])
    else:
        dataRAW['lum'] = 20000
    dataRAW['eff'] = dataRAW.mass/dataRAW.lum
    dataRAW['err'] = np.sqrt(dataRAW.mass)/dataRAW.lum
    plt.errorbar(data=dataRAW, x='label', y='eff', yerr='err', fmt='o', label=label)
