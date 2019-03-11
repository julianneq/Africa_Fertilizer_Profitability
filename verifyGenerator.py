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

    lastPlantDate = np.where(dates==datetime(2010,int(database[row,3]),int(database[row,4])))[0][0]
    remainingDays = len(dates) - lastPlantDate + 1
    if remainingDays < 30:
        month1 = np.zeros([30,2])
        month2 = np.zeros([30,2])
        month3 = np.zeros([30,2])
    elif remainingDays >= 30 and remainingDays < 60:
        month1 = np.zeros([31,2])
        month2 = np.zeros([30,2])
        month3 = np.zeros([30,2])
    elif remainingDays >= 60 and remainingDays < 90:
        month1 = np.zeros([31,2])
        month2 = np.zeros([31,2])
        month3 = np.zeros([30,2])
    else:
        month1 = np.zeros([31,2])
        month2 = np.zeros([31,2])
        month3 = np.zeros([31,2])
        
    # populate monthly data matrices
    for i in range(np.shape(month1)[0]):
        day1 = np.where(dates==datetime(1980+i,int(database[row,3]),int(database[row,4])))[0][0]
        month1[i,0] = np.mean(data[day1:(day1+30),0])
        month1[i,1] = np.sum(data[day1:(day1+30),1])
        
    for i in range(np.shape(month2)[0]):
        day1 = np.where(dates==datetime(1980+i,int(database[row,3]),int(database[row,4])))[0][0]
        month2[i,0] = np.mean(data[(day1+30):(day1+60),0])
        month2[i,1] = np.sum(data[(day1+30):(day1+60),1])
        
    for i in range(np.shape(month3)[0]):
        day1 = np.where(dates==datetime(1980+i,int(database[row,3]),int(database[row,4])))[0][0]
        month3[i,0] = np.mean(data[(day1+60):(day1+90),0])
        month3[i,1] = np.sum(data[(day1+60):(day1+90),1])
            
    return month1, month2, month3
        

for file in fileList:
    try:
        obsWeather = np.loadtxt('Weather/agmerra/txt/' + file, skiprows=4, usecols=[4,5])
        
        latLong = file.split("_")
        lat = float(latLong[0][3::])
        long = float(latLong[1][4:-4])
        
        row = np.intersect1d(np.where(database[:,1] == lat)[0], np.where(database[:,2] == long)[0])
        
        if len(row) > 0:
            obsMonth1, obsMonth2, obsMonth3 = getWeather(obsWeather, database, row)
            simPrecip = np.loadtxt('Weather/synthetic/' + file[0:-4] + '_precip.txt')
            simTemp = np.loadtxt('Weather/synthetic/' + file[0:-4] + '_temp.txt')
            
            fig = plt.figure()
            ax = fig.add_subplot(3,2,1) # temp on left, precip on right; rows = months
            ax.hist(obsMonth1[:,0], color='b', alpha=0.5, density=True)
            ax.hist(simTemp[:,0], color='g', alpha=0.5, density=True)
            ax.set_title('Temperature Month 1')
            
            ax = fig.add_subplot(3,2,2) # temp on left, precip on right; rows = months
            ax.hist(obsMonth1[:,1], color='b', alpha=0.5, density=True)
            ax.hist(simPrecip[:,0], color='g', alpha=0.5, density=True)
            ax.set_title('Precipitation Month 1')
            
            ax = fig.add_subplot(3,2,3) # temp on left, precip on right; rows = months
            ax.hist(obsMonth2[:,0], color='b', alpha=0.5, density=True)
            ax.hist(simTemp[:,1], color='g', alpha=0.5, density=True)
            ax.set_title('Temperature Month 2')
            
            ax = fig.add_subplot(3,2,4) # temp on left, precip on right; rows = months
            ax.hist(obsMonth2[:,1], color='b', alpha=0.5, density=True)
            ax.hist(simPrecip[:,1], color='g', alpha=0.5, density=True)
            ax.set_title('Precipitation Month 2')
            
            ax = fig.add_subplot(3,2,5) # temp on left, precip on right; rows = months
            ax.hist(obsMonth3[:,0], color='b', alpha=0.5, density=True)
            ax.hist(simTemp[:,2], color='g', alpha=0.5, density=True)
            ax.set_title('Temperature Month 3')
            
            ax = fig.add_subplot(3,2,6) # temp on left, precip on right; rows = months
            ax.hist(obsMonth3[:,1], color='b', alpha=0.5, density=True, label='Historical')
            ax.hist(simPrecip[:,2], color='g', alpha=0.5, density=True, label='Synthetic')
            ax.set_title('Precipitation Month 3')
            
            handles, labels = plt.gca().get_legend_handles_labels()
            fig.subplots_adjust(bottom=0.15)
            fig.legend(handles, labels, loc='lower center',ncol=2,fontsize=16,frameon=True)
            fig.set_size_inches([10.5,11.5])
            fig.savefig('Weather/validation/' + file[0:-4] + '.png')
            fig.clf()
        
    except:
        pass