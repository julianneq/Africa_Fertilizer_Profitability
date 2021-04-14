import numpy as np
from netCDF4 import Dataset
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib as mpl
import os

os.environ['PROJ_LIB']='C:/ProgramData/Anaconda3/pkgs/proj4-4.9.3-vc14_5/Library/share' # Work Desktop

from mpl_toolkits.basemap import Basemap

# read in grid points
gridPts = pd.read_csv('../AfricaDatabase.csv',delimiter=',')
gridPts['ELEV'] = np.zeros(np.shape(gridPts)[0])-9999

# read in elevation data
elev_data = Dataset('elev.0.25-deg.nc')

elev = elev_data.variables['data'][:]
lat = elev_data.variables['lat'][:]
lon = elev_data.variables['lon'][:]

# convert lon from (.125, 359.875) to (-179.875, 179.875)
lon[np.where(lon>180)[0]] = -180 + (lon[np.where(lon>180)[0]] - 180.0)

for j in range(np.shape(gridPts)[0]):
    # find which 0.25 degree elevation grid cell each Africa database cell is in
    lat_row = np.where(lat==gridPts['LAT'][j])[0][0]
    lon_row = np.where(lon==gridPts['LON'][j])[0][0]
    gridPts['ELEV'][j] = elev[0,lat_row,lon_row]

# write database to file
gridPts.to_csv('AfricaDatabase.csv')

# plot elevation on a map
def get_array(data, lat, lon, variable, bounds):
    array = np.empty((80*4,75*4))
    array[:] = np.NAN;
    
    for i in range(np.shape(data)[0]):
        row = int((data['LAT'][i]-np.min(lat))/0.25)
        col = int((data['LON'][i]-np.min(lon))/0.25)
        if np.isnan(data[variable][i]) == False:
                array[row,col] = data[variable][i]
                
    # update bounds on lower and upper end
    bounds[0] = np.min([bounds[0], np.nanmin(array)])
    bounds[-1] = np.max([bounds[-1], np.nanmax(array)])
    
    return array, bounds

m = Basemap(llcrnrlat=-40,urcrnrlat=40,llcrnrlon=-20,urcrnrlon=55,resolution='l')
m.drawlsmask()
m.drawcountries(color='0.2', linewidth=0.5)
m.drawparallels(np.arange(-40,41,20), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.8', fontsize=16)
m.drawmeridians(np.arange(-20,56,15), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.8', fontsize=16)

elev_array, bounds = get_array(gridPts, np.arange(-39.875,40.125,0.25), np.arange(-19.875,55.125,0.25), 'ELEV', np.array([0.0,1.0]))

array_mask = np.ma.masked_invalid(elev_array)
x, y = np.meshgrid(np.arange(-19.875,55.125,0.25), np.arange(-39.875,40.125,0.25)) # lons, lats
x, y = m(x,y)
m.pcolormesh(x, y, array_mask, cmap='terrain', rasterized=False, edgecolor='0.6', linewidth=0)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Elevation (m)')
plt.savefig('ElevationMap.png')
plt.clf()
