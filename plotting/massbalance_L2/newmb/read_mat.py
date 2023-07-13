import numpy as np
import h5py
import pandas as pd

#outpath = '/data/project3/kesf/sharing/budget_ww/outputs_N_PNDN_only_realistic_0-200m/'
#outpath = './budget_ww/outputs_C_PNDN_only_realistic_0-200m/'
outpath = './budget_ww/outputs_O2_PNDN_only_realistic_0-200m/'
exp = ['L2SCB',
       'L2SCB_AP',
       'PNDN_only_realistic',
       'FNDN_only_realistic',
       'pndn50_realistic',
       'pndn90_realistic']


matb = h5py.File(outpath+'MATBGCF.mat','r')
matc = h5py.File(outpath+'MATVARC.mat','r')

# data is now dictionary like
# data.keys() to see all keys
datab = matb.get('MATBGCF')
datac = matc.get('MATVARC')

#nh4up = np.squeeze(data['PHOTO_NH4'])

