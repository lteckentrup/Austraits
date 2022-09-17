import pandas as pd
import numpy as np

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
            PHEN.append(np.nan)

    ### Plant growth form
    for i in TN:
        try:
            PGF.append(df_PGF[df_PGF.taxon_name == i].loc[:,'value'].mode()[0])
        except KeyError:
            PGF.append(np.nan)

    ### Photosynthetic pathway
    for i in TN:
        try:
            PHOT.append(df_PHOT[df_PHOT.taxon_name == i].loc[:,'value'].mode()[0])
        except KeyError:
            PHOT.append(np.nan)

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
        except IndexError:
            lats.append(np.nan)

    ### Get longitude
    for sn,di in zip(SN,DI):
        try:
            lon = df_lons[(df_lons.site_name == sn) &
                          (df_lons.dataset_id == di)].loc[:,'value'].iloc[0]
            lons.append(float(lon))
        except IndexError:
            lons.append(np.nan)

    df_trait['lat'] = lats
    df_trait['lon'] = lons

    return (df_trait)

df_SLA = grab_trait('specific_leaf_area')
df_SLA = grab_trait('leaf_C_per_dry_mass')

df_SLA_filter = df_SLA[['value','PHEN','PGF','PHOT','LeafType','lat','lon']].dropna()

df_SLA_herb = df_SLA_filter[df_SLA_filter['PGF'].str.contains('herb|graminoid_not_tussock')]
df_SLA_tree = df_SLA_filter[df_SLA_filter['PGF'].str.contains('tree')]
