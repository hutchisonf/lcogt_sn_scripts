import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from pathlib import Path
from astropy.time import Time
from astropy.table import QTable
import astropy.units as u
import requests

field_observations_path = ''
supernova_observation_path = ''
color_term_path = None
pan_search_radius = 2

field_observations_path = Path(field_observations_path)
supernova_observations_path = Path(supernova_observation_path)

# field_spec = np.dtype({'names': [
#     'Name', 'RA', 'Dec', 'Time', 'Filter', 'Mag', 'Mag err'
#     ], 'formats': [
#     'U32',float, float, float, 'U1', float, float
# ]})

# supernova_spec = np.dtype({'names': [
#     'Name', 'RA', 'Dec', 'Time', 'Filter', 'Mag', 'Mag err'
#     ], 'formats': [
#     'U32',float, float, float, 'U1', float, float
# ]})

field_observations = np.load(field_observations_path)
supernova_observations = np.load(supernova_observations_path)

fields = QTable([
    field_observations['Name'],
    field_observations['RA']*u.deg,
    field_observations['Dec']*u.deg,
    Time(field_observations['Time'],format='mjd'),
    field_observations['Filter'],
    field_observations['Mag']*u.mag,
    field_observations['Mag err']*u.mag
])

supernovae = QTable([
    supernova_observations['Name'],
    supernova_observations['RA']*u.deg,
    supernova_observations['Dec']*u.deg,
    Time(supernova_observations['Time'],format='mjd'),
    supernova_observations['Filter'],
    supernova_observations['Mag']*u.mag,
    supernova_observations['Mag err']*u.mag
])

fields = fields[fields['Dec'] >= -30*u.deg]
supernovae = supernovae[supernovae['Dec'] >= -30*u.deg]

field_names, inds = np.unique(fields['Name'],return_index=True)
unique_fields = fields[inds]

# Perform PanSTARRS Search
url = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv"
data = {}
data['radius'] = pan_search_radius/3600
data['nDetections.gt'] = 1
columns = """objID,raMean,decMean,nDetections,ng,nr,ni,nz,ny,
    gMeanPSFMag,rMeanPSFMag,iMeanPSFMag,zMeanPSFMag,yMeanPSFMag,
    gMeanPSFMagErr,rMeanPSFMagErr,iMeanPSFMagErr,zMeanPSFMagErr,yMeanPSFMagErr""".split(',')
columns = [x.strip() for x in columns]
columns = [x for x in columns if x and not x.startswith('#')]
data['columns'] = '[{}]'.format(','.join(columns))
results = []
for f in unique_fields:
    data['ra'] = f['RA']
    data['dec'] = f['Dec']
    r = requests.get(url, params=data)
    r.raise_for_status()
    if r.text == '':
        tab = []
    else:
        tab = ascii.read(r.text)
        for f in 'grizy':
            col = f + 'MeanPSFMag'
            tab[col][tab[col] == -999.0] = np.nan
            col = f + 'MeanPSFMagErr'
            tab[col][tab[col] == -999.0] = np.nan
    results += [tab]
# Every item of the results list is either an empty list or an astropy table with the results