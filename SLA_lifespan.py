import pandas as pd
import numpy as np
from sklearn import linear_model
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error

'''
Derive leaf lifespan from SLA based on Reich et al. 1992. 
Most PFTs in LPJ do not have matching leaf lifespan but SLA is 
available for all PFTs. Leaf lifespan is also needed to derive
cton_min.
'''

fig = plt.figure(figsize=(9.5,6))

fig.subplots_adjust(hspace=0.22)
fig.subplots_adjust(wspace=0.2)
fig.subplots_adjust(right=0.95)
fig.subplots_adjust(left=0.12)
fig.subplots_adjust(bottom=0.35)
fig.subplots_adjust(top=0.95)

plt.rcParams['text.usetex'] = False
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['font.size'] = 11
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11

ax1=fig.add_subplot(1,2,1)
ax2=fig.add_subplot(1,2,2)

def get_TheilSen_regression(leaftype):
  
    ### Get SLA and leaf lifespan data
    df = pd.read_csv(leaftype+'_calcsla.csv')
    leaflong = np.array(df.leaflong)
    SLA = np.array(df.SLA)

    ### Fit Theil Sen model: robust against outliers
    TheilSen = linear_model.TheilSenRegressor(random_state=42)
    TheilSen.fit(np.log10(leaflong.reshape((-1, 1))),np.log10(SLA))
    return(TheilSen.coef_,TheilSen.intercept_,
           df.leaflong,df.SLA)

def plot_regression(leaftype,axis):
    slope, intercept, leaflong, SLA = get_TheilSen_regression(leaftype)
    print(slope)
    print(intercept)
    axis.scatter(leaflong,SLA,color='#00678a')
    fake_SLA = (slope*np.log10(leaflong))+intercept

    ### Default values LPJ
    if leaftype == 'broadleaf':
        SLA_LPJ = (-0.38*np.log10(leaflong))+2.41
    elif leaftype == 'needle':
        SLA_LPJ = (-0.4*np.log10(leaflong))+2.29
    
    ### Plot LPJ regression
    axis.plot(leaflong,10**SLA_LPJ,color='#984464')
    
    ### Plot new regression
    axis.plot(leaflong,10**fake_SLA,color='#00678a')
     
    ### Use logscales
    axis.set_xscale('log')
    axis.set_yscale('log')

plot_regression('broadleaf',ax1)
plot_regression('needle',ax2)

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)

plt.show()
