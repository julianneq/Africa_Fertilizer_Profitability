import numpy as np
import pandas as pd
import datetime

gridCells = np.loadtxt('NonZeroMaizeDateCells.csv', skiprows=1, delimiter=',')
maizeDatePts = np.loadtxt('NonZeroMaizeDateJoin.csv', skiprows=1, delimiter=',')
df = pd.DataFrame({'year': [2000]*np.shape(maizeDatePts)[0],
                  'doy': maizeDatePts[:,1]})
df.doy = df.doy.astype(int)
df = pd.to_datetime(df['year']*1000 + df['doy'], format='%Y%j')

gridCellMonths = np.zeros(np.shape(gridCells)[0])
gridCellDays = np.zeros(np.shape(gridCells)[0])
for i in range(len(gridCellDays)):
    indices = [k for k in range(np.shape(maizeDatePts)[0]) if maizeDatePts[k,2] == gridCells[i,1] and \
        maizeDatePts[k,3] == gridCells[i,2]]
    timestamp = np.sort(df[indices])[int(len(df[indices])/2)]
    gridCellMonths[i] = datetime.datetime.utcfromtimestamp(timestamp.tolist()/1e9).month
    gridCellDays[i] = datetime.datetime.utcfromtimestamp(timestamp.tolist()/1e9).day
    
newGridCells = np.zeros([np.shape(gridCells)[0],np.shape(gridCells)[1]+2])
newGridCells[:,0:np.shape(gridCells)[1]] = gridCells
newGridCells[:,-2] = gridCellMonths
newGridCells[:,-1] = gridCellDays

np.savetxt('NonZeroMaizeDayMonthCells.csv',newGridCells,delimiter=',')
