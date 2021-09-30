import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib as mpl
import os
import copy

#os.environ['PROJ_LIB']='C:/Users/Julianne Quinn/Anaconda3/pkgs/proj4-5.2.0-ha925a31_1/Library/share' # Home Desktop
#os.environ['PROJ_LIB']='C:/Users/jdq21/Anaconda3/pkgs/proj4-4.9.3-vc14_5/Library/share' # Laptop
os.environ['PROJ_LIB']='C:/ProgramData/Anaconda3/pkgs/proj4-4.9.3-vc14_5/Library/share' # Work Desktop

from mpl_toolkits.basemap import Basemap

# load simulation data
simData = pd.read_csv("Sensitivity_cfrobust.csv")
data = pd.read_csv('Simdata_cellid.csv')

# input variables to model
variables = ['temp_p1','temp_p2','temp_p3','precip_p1','precip_p2','precip_p3',
             'soilcec','soilph','acidity','soilom','soiln','claypct','siltpct','bulkdens',
             'elevm','mp','up','int']

# outputs of interest (yield different and investment return)
outputs = ['yp_dif','irr']

# create data frames to store delta changes in irr or yp_dif for each input variable
irrDeltas = pd.DataFrame(columns = ['cellid','iter','irr_baseline'] + variables)
yp_difDeltas = pd.DataFrame(columns = ['cellid','iter','yp_dif_baseline'] + variables[0:-3])

# copy cell id from simData
irrDeltas['cellid'] = simData['cellid']
yp_difDeltas['cellid'] = simData['cellid']

# copy iter from simData
irrDeltas['iter'] = simData['iter']
yp_difDeltas['iter'] = simData['iter']

# copy baseline irr and yp_dif from simData
irrDeltas['irr_baseline'] = simData['irr_baseline']
yp_difDeltas['yp_dif_baseline'] = simData['yp_dif_baseline']

# calculate absolute value of output difference for +-5% change in each input
for var in variables:
    irrDeltas[var] = np.abs(simData['irr_' + var + '_hi'] - simData['irr_' + var + '_lo'])
    if var != 'mp' and var != 'up' and var != 'int':
        yp_difDeltas[var] = np.abs(simData['yp_dif_' + var + '_hi'] - simData['yp_dif_' + var + '_lo'])
        
# remove cells with nan for irr because we don't have price data
irrDeltas.dropna(subset=['irr_baseline'], inplace=True)
        
# sort simulated irr and yp_dif values for each cell
irrDeltas = irrDeltas.sort_values(['cellid','irr_baseline'], ascending=[True, True])
yp_difDeltas = yp_difDeltas.sort_values(['cellid','yp_dif_baseline'], ascending=[True, True])

# convert variables to a numpy array and use it to find to which variable each model run in each cell was most sensitive
variables = np.array(variables)
irrDeltas['MostSensitive'] = variables[np.argmax(np.array(irrDeltas.iloc[:,3::]),axis=1)]
yp_difDeltas['MostSensitive'] = variables[np.argmax(np.array(yp_difDeltas.iloc[:,3::]),axis=1)]


# find variable to which each cell is most frequently most sensitive
irrMostSensitive = irrDeltas.groupby(['cellid'])['MostSensitive'].agg(pd.Series.mode)
yp_difMostSensitive = yp_difDeltas.groupby(['cellid'])['MostSensitive'].agg(pd.Series.mode)

irrMostSensitive.to_csv('irr_most_sensitive.csv')
yp_difMostSensitive.to_csv('yp_dif_most_sensitive.csv')


# find average delta for each variable across 1000 simulations for each cell
irrAvgDeltas = irrDeltas.groupby(['cellid'])[variables].mean()
yp_difAvgDeltas = yp_difDeltas.groupby(['cellid'])[variables[0:-3]].mean()

# rank sensitivities based on average delta
irrAvgDeltas['MostSensitive'] = variables[np.argmax(np.array(irrAvgDeltas),axis=1)]
yp_difAvgDeltas['MostSensitive'] = variables[np.argmax(np.array(yp_difAvgDeltas),axis=1)]

irrAvgDeltas['MostSensitive'].to_csv('irr_avg_most_sensitive.csv')
yp_difAvgDeltas['MostSensitive'].to_csv('yp_dif_avg_most_sensitive.csv')


# make sensitivity map for each percentile of yp_dif/irr
# simulated IRR distributions at different sites
shortlist_data = pd.read_csv('TrialSites_Sim_Shortlist.csv')
selectSites = np.unique(shortlist_data['sitecode'])
site_names = selectSites.astype('str')
index = range(len(np.where(shortlist_data['sitecode']==16)[0]))
CDFs = pd.DataFrame(index=index, columns=site_names)
for i in range(len(selectSites)):
    CDFs[site_names[i]] = shortlist_data['irr_sim'].iloc[np.where(shortlist_data['sitecode']==selectSites[i])[0]].to_numpy()

# locations of sites with simulated IRRs
lat = np.arange(-39.875,40.125,0.25)
lon = np.arange(-19.875,55.125,0.25)
select_data = pd.DataFrame(index=range(len(site_names)), columns=data.columns)
for i in range(len(selectSites)):
    row = i*len(np.where(shortlist_data['sitecode']==shortlist_data['sitecode'][0])[0]) # row of shortlist_data where new sitecode is located, every n simulations
    index = np.intersect1d(np.where(data['lat']==shortlist_data['lat'][row])[0], np.where(data['lon']==shortlist_data['lon'][row])[0])
    select_data.iloc[i,:] = data.iloc[index,:].to_numpy()
    
# color of line for each site and where to label it on the map
siteColors = ['#fdbf6f','#1f78b4','#b2df8a','#a6cee3','#cab2d6','#ff7f00','#fb9a99','#6a3d9a','#e31a1c','#33a02c']
label_ys = np.array(select_data['lat']-1)
label_xs = np.array(select_data['lon']-1.5)
label_xs[0] = label_xs[0] - 2.5
label_xs[3] = label_xs[3] - 5
label_xs[5] = label_xs[5] - 6
label_xs[[1,2,8]] = label_xs[[1,2,8]] + 3
label_ys[[2,4,6,7]] = label_ys[[2,4,6,7]] + 2
label_ys[9] = label_ys[9] - 3

def getGroup(sensitivity):
    # classify most sensitive variables into groups
    sensitivity['group'] = np.empty(len(sensitivity['MostSensitive']))
    sensitivity['group'][sensitivity.index[np.hstack([np.where(sensitivity['MostSensitive']=='temp_p1')[0],
                                    np.where(sensitivity['MostSensitive']=='temp_p2')[0],
                                    np.where(sensitivity['MostSensitive']=='temp_p3')[0]])]] = 0
    sensitivity['group'][sensitivity.index[np.hstack([np.where(sensitivity['MostSensitive']=='precip_p1')[0],
                                    np.where(sensitivity['MostSensitive']=='precip_p2')[0],
                                    np.where(sensitivity['MostSensitive']=='precip_p3')[0]])]] = 1
    sensitivity['group'][sensitivity.index[np.where(sensitivity['MostSensitive']=='soilph')[0]]] = 2
    sensitivity['group'][sensitivity.index[np.hstack([np.where(sensitivity['MostSensitive']=='soilcec')[0],
                                    np.where(sensitivity['MostSensitive']=='acidity')[0],
                                    np.where(sensitivity['MostSensitive']=='soilom')[0],
                                    np.where(sensitivity['MostSensitive']=='soiln')[0],
                                    np.where(sensitivity['MostSensitive']=='claypct')[0],
                                    np.where(sensitivity['MostSensitive']=='siltpct')[0],
                                    np.where(sensitivity['MostSensitive']=='bulkdens')[0]])]] = 3
    sensitivity['group'][sensitivity.index[np.where(sensitivity['MostSensitive']=='elevm')[0]]] = 4
    sensitivity['group'][sensitivity.index[np.hstack([np.where(sensitivity['MostSensitive']=='mp')[0],
                                    np.where(sensitivity['MostSensitive']=='up')[0],
                                    np.where(sensitivity['MostSensitive']=='int')[0]])]] = 5
    
    return sensitivity

def makeSensitivityMap(yp_difData, irrData, lat, lon, selectSites, select_data, label_xs, label_ys, percentile, figname):    
    C_yp_dif = np.array([[228,26,28],[55,126,184],[255,127,0],[166,86,40],[255,255,51]])/255.0
    classes_yp_dif = ['Temperature','Precipitation','Soil pH','Other Soil\nVariables','Elevation']
    bounds_yp_dif = np.array([-0.5,0.5,1.5,2.5,3.5,4.5])
    
    C_irr = np.array([[228,26,28],[55,126,184],[255,127,0],[166,86,40],[255,255,51],[77,175,74]])/255.0
    classes_irr = ['Temperature','Precipitation','Soil pH','Other Soil\nVariables','Elevation','Prices']
    bounds_irr = np.array([-0.5,0.5,1.5,2.5,3.5,4.5,5.5])
    
    class_array_yp_dif, bounds_yp_dif = get_array(yp_difData, lat, lon, 'group', bounds_yp_dif)
    class_array_irr, bounds_irr = get_array(irrData, lat, lon, 'group', bounds_irr)
    
    sns.set()
    fig = plt.figure()
    ax = fig.add_subplot(1,2,1)
    makeSubplot(ax, 'group', class_array_yp_dif, bounds_yp_dif, selectSites, select_data, label_xs, label_ys, \
                str(percentile) + '%-ile yield difference', C_yp_dif, '', classes_yp_dif, True, 1.0)
    
    ax = fig.add_subplot(1,2,2)
    makeSubplot(ax, 'group', class_array_irr, bounds_irr, selectSites, select_data, label_xs, label_ys, \
                str(percentile) + '%-ile IRR', C_irr, '', classes_irr, True, 1.0)
    
    fig.set_size_inches([15,5])
    fig.savefig(figname)
    fig.clf()
    
    return None

def get_array(data, lat, lon, variable, bounds):
    array = np.empty((80*4,75*4))
    array[:] = np.NAN;
    
    for i in range(np.shape(data)[0]):
        row = int((data['lat'][i]-np.min(lat))/0.25)
        col = int((data['lon'][i]-np.min(lon))/0.25)
        if variable != 'class' and variable != 'classnrt':
            if np.isnan(data[variable][i]) == False:
                array[row,col] = data[variable][i]
        elif variable == 'class':
            if data['class_yesyes_probT'][i] == 1:
                array[row,col] = 1 # always
            elif data['class_nono_probT'][i] == 1:
                array[row,col] = 0 # never
            elif data['class_type1_probT'][i] == 1:
                array[row,col] = 2 # robust
            elif data['class_type2_probT'][i] == 1:
                array[row,col] = 3 # naive
        elif variable == 'classnrt':
            if data['classnrt_yesyes_probT'][i] == 1:
                array[row,col] = 1 # always
            elif data['classnrt_nono_probT'][i] == 1:
                array[row,col] = 0 # never
            elif data['classnrt_type1_probT'][i] == 1:
                array[row,col] = 2 # robust
            elif data['classnrt_type2_probT'][i] == 1:
                array[row,col] = 3 # naive
                
    # update bounds on lower and upper end
    bounds[0] = np.round(np.min([bounds[0], np.nanmin(array)]),1)
    bounds[-1] = np.round(np.max([bounds[-1], np.nanmax(array)]),1)
    
    return array, bounds

def makeSubplot(ax, variable, array, bounds, selectSites, select_data, label_xs, label_ys, title, \
                C, label, classes, colorbar, shrink, points=False):
    
    # plot basemap and countries
    m = Basemap(llcrnrlat=-40,urcrnrlat=40,llcrnrlon=-20,urcrnrlon=55,resolution='l')
    m.drawlsmask()
    m.drawcountries(color='0.2', linewidth=0.5)
    m.drawparallels(np.arange(-40,41,20), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.8', fontsize=16)
    m.drawmeridians(np.arange(-20,56,15), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.8', fontsize=16)
        
    array_mask = np.ma.masked_invalid(array)
    x, y = np.meshgrid(lon, lat)
    x, y = m(x,y)
    if type(C) != str:
        cmap = mpl.colors.ListedColormap(C)
    else:
        cmap = C
        
    if variable == 'irr_sim_probT_cfrobust' or variable == 'irr_sim_probT_naive' or variable == 'cr_sim_probT_cfrobust':
        norm = mpl.colors.BoundaryNorm(boundaries=bounds,ncolors=(len(bounds)-1))
        m.pcolormesh(x, y, array_mask, cmap=cmap, norm=norm, rasterized=False, edgecolor='0.6', linewidth=0)
        cbar = plt.colorbar(shrink=shrink)
    else:
        m.pcolormesh(x, y, array_mask, cmap=cmap, rasterized=False, edgecolor='0.6', linewidth=0)
        if colorbar == True and variable != 'class' and variable != 'classnrt' and variable != 'group':
            cbar = plt.colorbar(shrink=shrink)
        if points == True:
            ax.scatter(select_data['lon'],select_data['lat'],c='k')
            for i in range(len(label_xs)):
                ax.text(label_xs[i], label_ys[i], selectSites[i], fontsize=16)
        
    if colorbar == True:
        if variable == 'cr_sim_probT_cfrobust':
            ticks = np.zeros([len(bounds)-1])
            for i in range(len(ticks)):
                ticks[i] = bounds[i] + (bounds[i+1] - bounds[i])/2
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(classes)
        elif variable == 'class' or variable == 'classnrt' or variable == 'group':
            formatter = plt.FuncFormatter(lambda val, loc: classes[val])
            cbar = plt.colorbar(ticks=np.arange(len(C)), format=formatter, shrink=shrink)
            plt.clim(-0.5,len(C)-0.5)
            cbar.solids.set_edgecolor("face")
        
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.set_ylabel(label,fontsize=18)
        
    ax.set_title(title, fontsize=18)
    
    return ax

for i in range(100):
    irrSubset = irrDeltas.iloc[(10*(i+1)-1)::1000]
    yp_difSubset = yp_difDeltas[(10*(i+1)-1)::1000]
    
    irrSensitivity = getGroup(copy.deepcopy(irrSubset))
    yp_difSensitivity = getGroup(copy.deepcopy(yp_difSubset))

    irrData = pd.merge(data,irrSensitivity)
    yp_difData = pd.merge(data,yp_difSensitivity)
    
    makeSensitivityMap(yp_difData, irrData, lat, lon, selectSites, select_data, label_xs, label_ys, i+1,
                       'gifFigures/SensitivityMap' + str(i) + '.png')
    
# make 3x2 panel figure for supplement of 10th, 50th and 90th percentile
C_yp_dif = np.array([[228,26,28],[55,126,184],[255,127,0],[166,86,40],[255,255,51]])/255.0
classes_yp_dif = ['Temperature','Precipitation','Soil pH','Other Soil\nVariables','Elevation']
bounds_yp_dif = np.array([-0.5,0.5,1.5,2.5,3.5,4.5])

C_irr = np.array([[228,26,28],[55,126,184],[255,127,0],[166,86,40],[255,255,51],[77,175,74]])/255.0
classes_irr = ['Temperature','Precipitation','Soil pH','Other Soil\nVariables','Elevation','Prices']
bounds_irr = np.array([-0.5,0.5,1.5,2.5,3.5,4.5,5.5])

sns.set()
fig = plt.figure()
percentile = [5, 50, 95]
for i in range(3):
    irrSubset = irrDeltas.iloc[(10*percentile[i]-1)::1000]
    yp_difSubset = yp_difDeltas[(10*percentile[i]-1)::1000]
    
    irrSensitivity = getGroup(copy.deepcopy(irrSubset))
    yp_difSensitivity = getGroup(copy.deepcopy(yp_difSubset))

    irrData = pd.merge(data,irrSensitivity)
    yp_difData = pd.merge(data,yp_difSensitivity)
    
    class_array_yp_dif, bounds_yp_dif = get_array(yp_difData, lat, lon, 'group', bounds_yp_dif)
    class_array_irr, bounds_irr = get_array(irrData, lat, lon, 'group', bounds_irr)

    ax = fig.add_subplot(2,3,i+1)
    if i == 2:
        makeSubplot(ax, 'group', class_array_yp_dif, bounds_yp_dif, selectSites, select_data, label_xs, label_ys, \
                str(int(percentile[i])) + 'th percentile yield difference', C_yp_dif, '', classes_yp_dif, True, 1.0)
    else:
        makeSubplot(ax, 'group', class_array_yp_dif, bounds_yp_dif, selectSites, select_data, label_xs, label_ys, \
                str(int(percentile[i])) + 'th percentile yield difference', C_yp_dif, '', classes_yp_dif, False, 1.0)
    
    ax = fig.add_subplot(2,3,i+4)
    if i == 2:
        makeSubplot(ax, 'group', class_array_irr, bounds_irr, selectSites, select_data, label_xs, label_ys, \
                str(int(percentile[i])) + 'th percentile IRR', C_irr, '', classes_irr, True, 1.0)
    else:
        makeSubplot(ax, 'group', class_array_irr, bounds_irr, selectSites, select_data, label_xs, label_ys, \
                str(int(percentile[i])) + 'th percentile IRR', C_irr, '', classes_irr, False, 1.0)

fig.set_size_inches([18,9])
fig.savefig('SensitivityMap_percentile.pdf')
fig.clf()

