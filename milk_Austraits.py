import pandas as pd
import numpy as np

df_traits = pd.read_csv('traits.csv',low_memory=False)

#### Select variable of interest
df_SLA = df_traits[df_traits['trait_name'] == 'specific_leaf_area']

### Grab info on phenology, growth form, and photosynthetic pathway
df_PHEN = df_traits[df_traits['trait_name'] == 'leaf_phenology']
df_PGF = df_traits[df_traits['trait_name'] == 'plant_growth_form']
df_PHOT = df_traits[df_traits['trait_name'] == 'photosynthetic_pathway']

### All taxon names into list
TN = df_SLA.taxon_name.to_list()

### Loop through all names and assign phenology, growth form,
### and photosyntetic pathway. Loops need to be separate to avoid double 
### counting of nan

PHEN = []
PGF = []
PHOT = []

for i in TN:
    try:
        PHEN.append(df_PHEN[df_PHEN.taxon_name == i].loc[:,'value'].mode()[0])
    except KeyError:
        PHEN.append(np.nan)

for i in TN:
    try:
        PGF.append(df_PGF[df_PGF.taxon_name == i].loc[:,'value'].mode()[0])
    except KeyError:
        PGF.append(np.nan)

for i in TN:
    try:
        PHOT.append(df_PHOT[df_PHOT.taxon_name == i].loc[:,'value'].mode()[0])
    except KeyError:
        PHOT.append(np.nan)
