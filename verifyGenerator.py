import numpy as np
import pandas as pd
from datetime import datetime
from matplotlib import pyplot as plt
import seaborn as sns

sns.set()

# range of historical data
dates = pd.date_range(start=datetime(1980,1,1),end=datetime(2010,12,31))
    
# load database with lat/long coordinates and associated plant dates
database = np.loadtxt('AfricaDatabase2.csv',delimiter=',',skiprows=1)
    
# read in sites with data and verify generator there
fileList = np.loadtxt('ValidationCells.csv',dtype=str)

def getWeather(data, database, row):
        
    # find historical monthly mean temperatures and total precipitation amounts in 3 consecutive 30-day periods post planting
    # first find how many of the 3 30-day periods post plant date are complete and size monthly data matrices accordingly
    # if plant date is Feb 29th, change to Feb 28th since there is no Feb 29th in 2010
    if database[row,3] == 2 and database[row,4] == 29:
        database[row,4] = 28

    lastPlantDate = np.where(dates==datetime(2010,int(database[row,3]),int(database[row,4])))[0][0] # planting date in last year
    remainingDays = len(dates) - lastPlantDate + 1 # number of days left in the calendar year
    # initialize all months with 30 years; add an extra year if there are enough remaining days in 2010 for a full month
    month1 = np.zeros([30,2])
    month2 = np.zeros([30,2])
    month3 = np.zeros([30,2])
    month4 = np.zeros([30,2])
    month5 = np.zeros([30,2])
    if remainingDays >= 30: # 1st month has an extra year of data
        month1 = np.zeros([31,2])
    elif remainingDays >= 60: # 1st 2 months have an extra year of data
        month2 = np.zeros([31,2])
    elif remainingDays >= 90: # 1st 3 months have an extra year of data
        month3 = np.zeros([31,2])
    elif remainingDays >= 120: # 1st 4 months have an extra year of data
        month4 = np.zeros([31,2])
    elif remainingDays >= 150: #all months have an extra year of data
        month5 = np.zeros([31,2])
        
    # populate monthly data matrices
    months = [month1, month2, month3, month4, month5]
    nMonths = len(months)
    for k in range(nMonths): # loop through 5 months
        for i in range(np.shape(months[k])[0]): # loop through 30 or 31 years
            day1 = np.where(dates==datetime(1980+i,int(database[row,3]),int(database[row,4])))[0][0] # planting date in year i
            months[k][i,0] = np.mean(data[(day1+30*k):(day1+30*(k+1)),0]) # mean temperature
            months[k][i,1] = np.sum(data[(day1+30*k):(day1+30*(k+1)),1]) # total precipitation
            
    return months

def makeTempSubplot(ax, months, tempData, monthNo):
    ax.hist(months[monthNo-1][:,0], color='b', alpha=0.5, density=True)
    ax.hist(tempData[:,monthNo-1], color='g', alpha=0.5, density=True)
    ax.set_title('Temperature Month ' + str(monthNo))
    
    return None
    
def makePrecipSubplot(ax, months, precipData, monthNo, label=False):
    if label == False:
        ax.hist(months[monthNo-1][:,1], color='b', alpha=0.5, density=True)
        ax.hist(precipData[:,monthNo-1], color='g', alpha=0.5, density=True)
    else:
        ax.hist(months[monthNo-1][:,1], color='b', alpha=0.5, density=True, label='Historical')
        ax.hist(precipData[:,monthNo-1], color='g', alpha=0.5, density=True, label='Synthetic')        
    ax.set_title('Precipitation Month ' + str(monthNo))
    
    return None
        

for file in fileList:
    try:
        obsWeather = np.loadtxt('Weather/agmerra/txt/' + file, skiprows=4, usecols=[4,5])
        
        latLong = file.split("_")
        lat = float(latLong[0][3::])
        long = float(latLong[1][4:-4])
        
        row = np.intersect1d(np.where(database[:,1] == lat)[0], np.where(database[:,2] == long)[0])
        
        if len(row) > 0:
            obsMonths = getWeather(obsWeather, database, row)
            simPrecip = np.loadtxt('Weather/synthetic/' + file[0:-4] + '_precip.txt')
            simTemp = np.loadtxt('Weather/synthetic/' + file[0:-4] + '_temp.txt')
            
            fig = plt.figure()
    
            for i in range(len(obsMonths)):
                ax = fig.add_subplot(2,5,i+1) # temp on top, precip on bottom; columns = months
                makeTempSubplot(ax, obsMonths, simTemp, i+1)
                
                ax = fig.add_subplot(2,5,i+6)
                if i != len(obsMonths)-1:
                    makePrecipSubplot(ax, obsMonths, simPrecip, i+1)
                else:
                    makePrecipSubplot(ax, obsMonths, simPrecip, i+1, True)
            
            handles, labels = plt.gca().get_legend_handles_labels()
            fig.subplots_adjust(bottom=0.1)
            fig.legend(handles, labels, loc='lower center',ncol=2,fontsize=16,frameon=True)
            fig.set_size_inches([19.2,9.6])
            fig.savefig('Weather/validation/' + file[0:-4] + '.png')
            fig.clf()
        
    except:
        pass
