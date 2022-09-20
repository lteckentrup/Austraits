import pandas as pd
import numpy as np
from sklearn import linear_model

df_AusTraits = pd.read_csv('traits.csv',low_memory=False)

def grab_trait(trait_name):
    #### Select variable of interest
    df_trait = df_AusTraits[df_AusTraits['trait_name'] == trait_name]
    df_trait = df_trait.reset_index()

    ### Grab info on phenology, growth form, and photosynthetic pathway
    df_PHEN = df_AusTraits[df_AusTraits['trait_name'] == 'leaf_phenology']
    df_PGF = df_AusTraits[df_AusTraits['trait_name'] == 'plant_growth_form']
    df_PHOT = df_AusTraits[df_AusTraits['trait_name'] == 'photosynthetic_pathway']
    df_LeafType = df_AusTraits[df_AusTraits['trait_name'] == 'leaf_type']

    ### All taxon names into list
    TN = df_trait.taxon_name.to_list()

    '''
    Loop through all names and assign phenology, growth form,
    and photosyntetic pathway. Loops need to be separate to avoid double
    counting of nan
    '''

    PHEN = []
    PGF = []
    PHOT = []
    LeafType = []

    ### Phenology
    for i in TN:
        try:
            PHEN.append(df_PHEN[df_PHEN.taxon_name == i].loc[:,'value'].mode()[0])
        except KeyError:
            PHEN.append('not_defined')

    ### Plant growth form
    for i in TN:
        try:
            PGF.append(df_PGF[df_PGF.taxon_name == i].loc[:,'value'].mode()[0])
        except KeyError:
            PGF.append('not_defined')

    ### Photosynthetic pathway
    for i in TN:
        try:
            PHOT.append(df_PHOT[df_PHOT.taxon_name == i].loc[:,'value'].mode()[0])
        except KeyError:
            PHOT.append('not_defined')

    ### Leaf type
    for i in TN:
        try:
            LeafType.append(df_LeafType[df_LeafType.taxon_name == i].loc[:,'value'].mode()[0])
        except KeyError:
            LeafType.append('not_defined')

    ### Add info to dataframe
    df_trait['PHEN'] = PHEN
    df_trait['PGF'] = PGF
    df_trait['PHOT'] = PHOT
    df_trait['LeafType'] = LeafType

    '''
    Now grab site information: latitude and longitude so traits can be grouped into
    temperate and tropical vegetation (LPJ-GUESS has temperatue and tropical PFTs)
    '''

    ### Read in site info
    df_sites = pd.read_csv('sites.csv')
    df_lats = df_sites[df_sites['site_property']=='latitude (deg)']
    df_lons = df_sites[df_sites['site_property']=='longitude (deg)']

    lats = []
    lons = []

    SN = df_trait.site_name.to_list()
    DI = df_trait.dataset_id.to_list()

    ### Get latitude
    for sn,di in zip(SN,DI):
        try:
            lat = df_lats[(df_lats.site_name == sn) &
                          (df_lats.dataset_id == di)].loc[:,'value'].iloc[0]
            lats.append(float(lat))
        except (IndexError,ValueError):
            lats.append(np.nan)

    ### Get longitude
    for sn,di in zip(SN,DI):
        try:
            lon = df_lons[(df_lons.site_name == sn) &
                          (df_lons.dataset_id == di)].loc[:,'value'].iloc[0]
            lons.append(float(lon))
        except (IndexError,ValueError):
            lons.append(np.nan)

    df_trait['lat'] = lats
    df_trait['lon'] = lons

    return (df_trait)

def get_wooddens(veggroup):
    ### AusTraits reports wood density in kg/m3 but LPJ GUESS takes kgC/m3
    df_wooddens = grab_trait('wood_density')
    df_woodC = grab_trait('wood_C_per_dry_mass')

    ### Filter dataframes for nan
    df_wooddens_veggroup = df_wooddens[df_wooddens['PGF'].str.contains(veggroup)]
    df_woodC_veggroup = df_woodC[df_woodC['PGF'].str.contains(veggroup)]

    obs_ID = df_woodC_veggroup.observation_id.to_list()

    wooddens = []
    lat = []
    for i in obs_ID:
        try:
            woodC = df_woodC_veggroup[df_woodC_veggroup.observation_id == i].loc[:,'value'].iloc[0]
            wooddens_kg = df_wooddens_veggroup[df_wooddens_veggroup.observation_id == i].loc[:,'value'].iloc[0]
            LAT = df_wooddens_veggroup[df_wooddens_veggroup.observation_id == i].loc[:,'lat'].iloc[0]
            wooddens.append(float(woodC)*float(wooddens_kg))
            lat.append(LAT)
        except IndexError:
            pass

    print(veggroup+' median')
    print(np.median(wooddens))
    print(veggroup+' 25 percentile')
    print(np.quantile(wooddens,0.25))
    print(veggroup+' 75 percentile')
    print(np.quantile(wooddens,0.75))   
    
get_wooddens('tree')
get_wooddens('shrub')
 
#-------------------------------------------------------------------------------
### Update calcsla function
   
def isla_params(leaftype):
    ### Define equation based on Reich et al., 1992 - calculate SLA from leaflong
    df_SLA = grab_trait('specific_leaf_area')
    df_lifespan = grab_trait('leaf_lifespan')

    ### Select broadleafed vegetation
    df_SLA_leaf = df_SLA[df_SLA['LeafType']==leaftype]
    df_lifespan_leaf = df_lifespan[df_lifespan['LeafType']==leaftype]

    SLA = []
    LIFESPAN = []
    obs_ID = df_SLA_leaf.observation_id.to_list()

    ### Grab data that match
    for i in obs_ID:
        try:
            sla = df_SLA_leaf[df_SLA_leaf.observation_id == i].loc[:,'value'].iloc[0]
            lifespan = df_lifespan_leaf[df_lifespan_leaf.observation_id == i].loc[:,'value'].iloc[0]

            SLA.append(float(sla)*10)
            LIFESPAN.append(float(lifespan))

        except (IndexError,KeyError):
            pass

    SLA = np.array(SLA)
    LIFESPAN = np.array(LIFESPAN)
    
    ### Export to csv and calculate regression separately in plotting script
    df = pd.DataFrame()
    df['leaflong'] = LIFESPAN
    df['SLA'] = SLA
    df.to_csv(leaftype+'_calcsla.csv')

isla_params('broadleaf')
isla_params('needle')
isla_params('broadleaf')
isla_params('needle')

#-------------------------------------------------------------------------------
### Get leaf lifespan
df_lifespan = grab_trait('leaf_lifespan')

### Grasses
print(df_lifespan[(df_lifespan['PGF']=='herb')|(df_lifespan['PGF']=='graminoid')])
### only six entries, either 4 or 24?? not enough, also no C4. Use SLA equation

### Tropical evergreen trees: no data
print('TrBE')
print(df_lifespan[(df_lifespan['PGF']=='tree')&
                  (df_lifespan['PHEN']=='evergreen')&
                  (df_lifespan['LeafType']=='broadleaf')&
                  (df_lifespan['lat']>-23.5)].value.dropna().astype(float).median()/12)

### Tropical raingreen trees: no data
print('TrBR')
print(df_lifespan[(df_lifespan['PGF']=='tree')&
                  (df_lifespan['PHEN']=='drought_deciduous')&
                  (df_lifespan['LeafType']=='broadleaf')&
                  (df_lifespan['lat']>-23.5)].value.dropna().astype(float).median()/12)


### Temperate broadleaf evergreen trees, 20 values
print('Median TeBE')
print(df_lifespan[(df_lifespan['PGF']=='tree')&
                  (df_lifespan['PHEN']=='evergreen')&
                  (df_lifespan['LeafType']=='broadleaf')&
                  (df_lifespan['lat']<-23.5)].value.dropna().astype(float).median()/12)

print('25 percentile TeBE')
print(df_lifespan[(df_lifespan['PGF']=='tree')&
                  (df_lifespan['PHEN']=='evergreen')&
                  (df_lifespan['LeafType']=='broadleaf')&
                  (df_lifespan['lat']<-23.5)].value.dropna().astype(float).quantile(0.25)/12)

print('75 percentile TeBE')
print(df_lifespan[(df_lifespan['PGF']=='tree')&
                  (df_lifespan['PHEN']=='evergreen')&
                  (df_lifespan['LeafType']=='broadleaf')&
                  (df_lifespan['lat']<-23.5)].value.dropna().astype(float).quantile(0.75)/12)

### Temperate broadleaf deciduous trees: no data
print(df_lifespan[(df_lifespan['PGF']=='tree')&
                  (df_lifespan['PHEN']=='deciduous')&
                  (df_lifespan['LeafType']=='broadleaf')&
                  (df_lifespan['lat']<-23.5)].value.dropna().astype(float).median()/12)

### Temperate needleleaf evergreen trees: no data
print(df_lifespan[(df_lifespan['PGF']=='tree')&
                  (df_lifespan['PHEN']=='evergreen')&
                  (df_lifespan['LeafType']=='needle')&
                  (df_lifespan['lat']<-30)].value.dropna().astype(float).median()/12)

### Grass: only six values, all c3, some 4 months, some 24?
print(df_lifespan[(df_lifespan['PGF']=='herb')|
                  (df_lifespan['PGF']=='graminoid')].value.dropna().astype(float).median()/12)

#-------------------------------------------------------------------------------
### For the missing data calculate leaf life span from specific leaf area.
### Lifespan is also needed for ctonmin so keep lifespan instead of updating
### sla following updated equation

### Get specific leaf area
df_SLA = grab_trait('specific_leaf_area')

### Grasses
df_grass = df_SLA[(df_SLA['PGF']=='herb')|(df_SLA['PGF']=='graminoid')]
print('Median C3 grass')
print(df_grass[df_grass['PHOT'] == 'c3'].value.dropna().astype(float).median())
print('25 percentile C3 grass')
print(df_grass[df_grass['PHOT'] == 'c3'].value.dropna().astype(float).quantile(0.25))
print('75 percentile C3 grass')
print(df_grass[df_grass['PHOT'] == 'c3'].value.dropna().astype(float).quantile(0.75))

print('Median C4 grass')
print(df_grass[df_grass['PHOT'] == 'c4'].value.dropna().astype(float).median())
print('25 percentile C4 grass')
print(df_grass[df_grass['PHOT'] == 'c4'].value.dropna().astype(float).quantile(0.25))
print('75 percentile 4 grass')
print(df_grass[df_grass['PHOT'] == 'c4'].value.dropna().astype(float).quantile(0.75))

### Tropical evergreen trees
print('Median TrBE')
print(df_SLA[(df_SLA['PGF']=='tree')&
             (df_SLA['PHEN']=='evergreen')&
             (df_SLA['LeafType']=='broadleaf')&
             (df_SLA['lat']>-23.5)].value.dropna().astype(float).median())
print('25 percentile TrBE')
print(df_SLA[(df_SLA['PGF']=='tree')&
             (df_SLA['PHEN']=='evergreen')&
             (df_SLA['LeafType']=='broadleaf')&
             (df_SLA['lat']>-23.5)].value.dropna().astype(float).quantile(0.25))
print('75 percentile TrBE')
print(df_SLA[(df_SLA['PGF']=='tree')&
             (df_SLA['PHEN']=='evergreen')&
             (df_SLA['LeafType']=='broadleaf')&
             (df_SLA['lat']>-23.5)].value.dropna().astype(float).quantile(0.75))

### Tropical evergreen trees. Leaftype undefined; checked and they are BDL
print('Median TrBR')
print(df_SLA[(df_SLA['PGF']=='tree')&
             (df_SLA['PHEN']=='drought_deciduous')&
             (df_SLA['lat']>-23.5)].value.dropna().astype(float).median())
### only six :(
print('25 percentile TrBR')
print(df_SLA[(df_SLA['PGF']=='tree')&
             (df_SLA['PHEN']=='drought_deciduous')&
             (df_SLA['lat']>-23.5)].value.dropna().astype(float).quantile(0.25))
print('75 percentile TrBR')
print(df_SLA[(df_SLA['PGF']=='tree')&
             (df_SLA['PHEN']=='drought_deciduous')&
             (df_SLA['lat']>-23.5)].value.dropna().astype(float).quantile(0.75))

### Temperate evergreen trees
print('Median TeBE')
print(df_SLA[(df_SLA['PGF']=='tree')&
             (df_SLA['PHEN']=='evergreen')&
             (df_SLA['LeafType']=='broadleaf')&
             (df_SLA['lat']<-23.5)].value.dropna().astype(float).median())
print('25 percentile TeBE')
print(df_SLA[(df_SLA['PGF']=='tree')&
             (df_SLA['PHEN']=='evergreen')&
             (df_SLA['LeafType']=='broadleaf')&
             (df_SLA['lat']<-23.5)].value.dropna().astype(float).quantile(0.25))
print('75 percentile TeBE')
print(df_SLA[(df_SLA['PGF']=='tree')&
             (df_SLA['PHEN']=='evergreen')&
             (df_SLA['LeafType']=='broadleaf')&
             (df_SLA['lat']<-23.5)].value.dropna().astype(float).quantile(0.75))

### Temperate summergreen trees
print(df_SLA[(df_SLA['PGF']=='tree')&
             (df_SLA['PHEN']=='deciduous')&
             (df_SLA['LeafType']=='broadleaf')&
             (df_SLA['lat']<-23.5)].value.dropna().astype(float).median())

### Temperate evergreen NDL trees
print('Median TeNE')
print(df_SLA[(df_SLA['PGF']=='tree')&
             (df_SLA['PHEN']=='evergreen')&
             (df_SLA['LeafType']=='needle')&
             (df_SLA['lat']<-23.5)].value.dropna().astype(float).median())
print('25 percentile TeNE')
print(df_SLA[(df_SLA['PGF']=='tree')&
             (df_SLA['PHEN']=='evergreen')&
             (df_SLA['LeafType']=='needle')&
             (df_SLA['lat']<-23.5)].value.dropna().astype(float).quantile(0.25))
print('75 percentile TeNE')
print(df_SLA[(df_SLA['PGF']=='tree')&
             (df_SLA['PHEN']=='evergreen')&
             (df_SLA['LeafType']=='needle')&
             (df_SLA['lat']<-23.5)].value.dropna().astype(float).quantile(0.75))

#-------------------------------------------------------------------------------

### Leaf to sapwood area ratio. Inverse of huber_value. 
### Drop drought affected vegetation
df_klatosa = grab_trait('huber_value')

### All trees
print('Median tree')
print(1/df_klatosa[(df_klatosa['context_name']!='Drought')&
                   (df_klatosa['PGF']=='tree')].value.dropna().astype(float).median())
print('25 percentile tree')
print(1/df_klatosa[(df_klatosa['context_name']!='Drought')&
                   (df_klatosa['PGF']=='tree')].value.dropna().astype(float).quantile(0.25))
print('75 percentile tree')
print(1/df_klatosa[(df_klatosa['context_name']!='Drought')&
                   (df_klatosa['PGF']=='tree')].value.dropna().astype(float).quantile(0.75))
### All shrubs
print('Median shrub')
print(1/df_klatosa[(df_klatosa['context_name']!='Drought')&
                   (df_klatosa['PGF']=='shrub')].value.dropna().astype(float).median())
print('25 percentile shrub')
print(1/df_klatosa[(df_klatosa['context_name']!='Drought')&
                   (df_klatosa['PGF']=='shrub')].value.dropna().astype(float).quantile(0.25))
print('75 percentile shrub')
print(1/df_klatosa[(df_klatosa['context_name']!='Drought')&
                   (df_klatosa['PGF']=='shrub')].value.dropna().astype(float).quantile(0.75))

#-------------------------------------------------------------------------------
'''
Reich, P.B., Walters, M.B. and Ellsworth, D.S. (1992), Leaf Life-Span in Relation to Leaf, Plant, 
and Stand Characteristics among Diverse Ecosystems. Ecological Monographs, 62: 365-392. 
https://doi.org/10.2307/2937116
'''
