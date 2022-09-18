'''
Uses AustTraits to update LPJ GUESS PFTs parameters and the calculation that translates
leaf lifespan to specific leaf area following

Reich, P. B., Walters, M. B., & Ellsworth, D. S. (1992). Leaf Life-Span in Relation to Leaf, 
Plant, and Stand Characteristics among Diverse Ecosystems. Ecological Monographs, 62(3), 
365–392. https://doi.org/10.2307/2937116
'''

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
            LeafType.append(np.nan)

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
    df_wooddens = grab_trait('wood_density')
    df_woodC = grab_trait('wood_C_per_dry_mass')

    ### Filter dataframes for nan
    df_wooddens_filter = df_wooddens[['observation_id','value','PGF']].dropna()
    df_woodC_filter = df_woodC[['observation_id','value','PGF']].dropna()

    df_wooddens_veggroup = df_wooddens_filter[df_wooddens_filter['PGF'].str.contains(veggroup)]
    df_woodC_veggroup = df_woodC_filter[df_woodC_filter['PGF'].str.contains(veggroup)]

    obs_ID = df_woodC_veggroup.observation_id.to_list()

    wooddens = []

    for i in obs_ID:
        try:
            woodC = df_woodC_veggroup[df_woodC_veggroup.observation_id == i].loc[:,'value'].iloc[0]
            wooddens_kg = df_wooddens_veggroup[df_wooddens_veggroup.observation_id == i].loc[:,'value'].iloc[0]
            wooddens.append(float(woodC)*float(wooddens_kg))
        except IndexError:
            pass

    print(veggroup)
    print(np.median(wooddens))

get_wooddens('tree')
get_wooddens('shrub')

def isla_params():
    ### Define equation based on Reich et al., 1992 - calculate SLA from leaflong
    df_SLA = grab_trait('specific_leaf_area')
    df_lifespan = grab_trait('leaf_lifespan')

    ### Select broadleafed vegetation
    df_SLA_BDL = df_SLA[df_SLA['PGF']=='broadleaf']
    df_lifespan_BDL = df_lifespan[df_lifespan['PGF']=='broadleaf']

    SLA = []
    LIFESPAN = []
    obs_ID = df_SLA_BDL.observation_id.to_list()

    ### Grab data that match
    for i in obs_ID:
        try:
            sla = df_SLA_BDL[df_SLA_BDL.observation_id == i].loc[:,'value'].iloc[0]
            lifespan = df_lifespan_BDL[df_lifespan_BDL.observation_id == i].loc[:,'value'].iloc[0]

            SLA.append(float(sla)*10)
            LIFESPAN.append(float(lifespan)*10)

        except (IndexError,KeyError):
            pass

    lifespan = np.array(LIFESPAN)
    ransac.fit(np.log10(lifespan.reshape((-1, 1))),np.log10(df.sla))
    
    print(ransac.estimator_.coef_)
    print(ransac.estimator_.intercept_)
