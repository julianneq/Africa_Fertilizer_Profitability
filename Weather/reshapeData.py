import numpy as np
from netCDF4 import Dataset

# variables:
# 'prate': 365/6 x 720 x 1440 (Precipitation Rate in mm/day, 0.25 degree x 0.25 degree)
# 'tmax': 365/6 x 720 x 1440 (Daily Maximum Temperature in degrees C, 0.5 degree x 0.5 degree)
# 'tmin': 365/6 x 720 x 1440 (Daily Minimum Temperature in degrees C, 0.5 degree x 0.5 degree)
# 'tavg': 365/6 x 720 x 1440 (Daily Average Temperature in degrees C, 0.5 degree x 0.5 degree)
# 'srad': 365/6 x 720 x 1440 (Solar Radiation in MJ/m^2/day, 1 degree x 1 degree)
# 'rhstmax': 365/6 x 720 x 1440 (Relative Humidity at Time of Maximum daily Maximum Temperature in %, 0.25 degree x 0.25 degree)
# 'wndspd': 365/6 x 720 x 1440 (Daily Average Wind Speed in m/s, 0.25 degree x 0.25 degree)
# 'time': 365/6 (Day of the Year)
# 'latitude': 720 (coordinates of center of grid cell: 89.875, 89.625, ... -89.625, -89.875)
# 'longitude': 1440 (coordinates of center of grid cell: 0.125, 0.375, ..., 359.625, 359.875)

# output format:
#*WEATHER DATA : 152_434.psims.nc
#@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT
#    CI   14.250   36.750   -99  29.2   6.3   -99    10
#@DATE  SRAD  TMAX  TMIN  TAVG  RAIN  WIND
#80001  19.8  33.4  18.1 25.8  0.0 146.9

# read in grid points
gridPts = np.loadtxt('AfricaGridCells.csv',delimiter=',',skiprows=1,usecols=[3,4,5,6])
tClimatology = np.zeros(np.shape(gridPts)[0])
ndays = np.zeros(np.shape(gridPts)[0])

years = range(1980,2011)
for i in range(len(years)):
    # read in data for year i
    precipData = Dataset('agmerra/netCDF/AgMERRA_' + str(years[i]) + '_prate.nc4')
    rain = precipData.variables['prate'][:]
    lat = precipData.variables['latitude'][:]
    lon = precipData.variables['longitude'][:]
    doy = precipData.variables['time'][:]
    tempMaxData = Dataset('agmerra/netCDF/AgMERRA_' + str(years[i]) + '_tmax.nc4')
    tmax = tempMaxData.variables['tmax'][:]
    tempMinData = Dataset('agmerra/netCDF/AgMERRA_' + str(years[i]) + '_tmin.nc4')
    tmin = tempMinData.variables['tmin'][:]
    tempAvgData = Dataset('agmerra/netCDF/AgMERRA_' + str(years[i]) + '_tavg.nc4')
    tavg = tempAvgData.variables['tavg'][:]
    sradData = Dataset('agmerra/netCDF/AgMERRA_' + str(years[i]) + '_srad.nc4')
    srad = sradData.variables['srad'][:]
    windData = Dataset('agmerra/netCDF/AgMERRA_' + str(years[i]) + '_wndspd.nc4')
    wind = windData.variables['wndspd'][:]
    
    for j in range(np.shape(gridPts)[0]):
        if i == 0:
            # find intial estimate of temperature climatology
            tClimatology[j] = np.mean(tavg[:,int(gridPts[j,2]),int(gridPts[j,3])])
            ndays[j] = len(doy)
        
            # start new file with reformatted data
            f = open('agmerra/txt/lat' + str(gridPts[j,0]) + '_long' + str(gridPts[j,1]) + '.txt','w')
            f.write('*WEATHER DATA:\n@ INSI LAT LONG ELEV TAV AMP REFHT WNDHT\n')
            f.write('CI   ' + str(gridPts[j,0]) + ' ' + str(gridPts[j,1]) + ' -99 ' + str(np.round(tClimatology[j],1)) + \
                ' -99 -99 10\n')
            f.write('@DATE SRAD TMAX TMIN TAVG RAIN WIND\n')
            f.close()
        else:
            # update estimate of temperature climatology
            tClimatology[j] = (tClimatology[j]*ndays[j] + np.mean(tavg[:,int(gridPts[j,2]),int(gridPts[j,3])])*len(doy)) / (ndays[j] + len(doy))
            ndays[j] += len(doy)
            
        # append reformatted data to existing file
        f = open('agmerra/txt/lat' + str(gridPts[j,0]) + '_long' + str(gridPts[j,1]) + '.txt','a')
        for k in range(len(doy)):
            if doy[k] < 10:
                str_doy = '00' + str(int(doy[k]))
            elif doy[k] >= 10 and doy[k] < 100:
                str_doy = '0' + str(int(doy[k]))
            else:
                str_doy = str(int(doy[k]))
                
            f.write(str(years[i])[2:4] + str_doy + ' ' + str(srad[k,int(gridPts[j,2]),int(gridPts[j,3])]) + '  ' + \
                str(tmax[k,int(gridPts[j,2]),int(gridPts[j,3])]) + ' ' + str(tmin[k,int(gridPts[j,2]),int(gridPts[j,3])]) + ' ' + \
                str(tavg[k,int(gridPts[j,2]),int(gridPts[j,3])]) + ' ' + \
                str(rain[k,int(gridPts[j,2]),int(gridPts[j,3])]) + ' ' + str(wind[k,int(gridPts[j,2]),int(gridPts[j,3])]) + '\n')
            
        f.close()
                
        #if i == len(years)-1:
            # re-write estimate of temperature climatology
        #    f = open('agmerra/txt/lat' + str(gridPts[j,0]) + '_long' + str(gridPts[j,1]) + '.txt','r+b')
        #    f_content = f.readlines()
        #    f_content[2] = 'CI   ' + str(gridPts[j,0]) + ' ' + str(gridPts[j,1]) + ' -99 ' + str(np.round(tClimatology[j],1)) + \
        #        ' -99 -99 10\n'
        #    f.seek(0)
        #    f.write(''.join(f_content))
        #    f.close()
