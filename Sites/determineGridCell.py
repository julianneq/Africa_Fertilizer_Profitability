import numpy as np
    
def conditions(i, index):
    direction = []
    if abs(sites[i,index] - int(sites[i,index])) == 0.0 or abs(sites[i,index] - int(sites[i,index])) == 0.25 \
        or abs(sites[i,index] - int(sites[i,index])) == 0.5 or abs(sites[i,index] - int(sites[i,index])) == 0.75:
        direction.append(str(sites[i,index]-0.125))
        direction.append(str(sites[i,index]+0.125))
    elif abs(sites[i,index] - int(sites[i,index])) > 0 and abs(sites[i,index] - int(sites[i,index])) < 0.25:
        if sites[i,index] < 0:
            direction.append(str(int(sites[i,index])-0.125))
        else:
            direction.append(str(int(sites[i,index])+0.125))
    elif abs(sites[i,index] - int(sites[i,index])) > 0.25 and abs(sites[i,index] - int(sites[i,index])) < 0.5:
        if sites[i,index] < 0:
            direction.append(str(int(sites[i,index])-0.375))
        else:
            direction.append(str(int(sites[i,index])+0.375))
    elif abs(sites[i,index] - int(sites[i,index])) > 0.5 and abs(sites[i,index] - int(sites[i,index])) < 0.75:
        if sites[i,index] < 0:
            direction.append(str(int(sites[i,index])-0.625))
        else:
            direction.append(str(int(sites[i,index])+0.625))
    else:
        if sites[i,index] < 0:
            direction.append(str(int(sites[i,index])-0.875))
        else:
            direction.append(str(int(sites[i,index])+0.875))
            
    return direction
    
sites = np.genfromtxt('EIL_site_lat_lon.csv',delimiter=',',skip_header=1,usecols=[4,5])
file1 = []
file2 = []
file3 = []
file4 = []
for i in range(np.shape(sites)[0]):
    lat = conditions(i, 0)
    lon = conditions(i, 1)
            
    file1.append('lat' + lat[0] + '_long' + lon[0] + '.txt')
    if len(lat) > 1:
        file2.append('lat' + lat[1] + '_long' + lon[0] + '.txt')
        if len(lon) > 1:
            file3.append('lat' + lat[0] + '_long' + lon[1] + '.txt')
            file4.append('lat' + lat[1] + '_long' + lon[1] + '.txt')
        else:
            file3.append('')
            file4.append('')
    else:
        if len(lon) > 1:
            file2.append('lat' + lat[0] + '_long' + lon[1] + '.txt')
        else:
            file2.append('')
            
        file3.append('')
        file4.append('')
        
f = open('EIL_site_lat_lon.csv','r')
y = []
for line in f.readlines():
    value = [value for value in line.split(',')]
    y.append(value)
    
f.close()

y[0][-1] = y[0][-1][0:-1] # remove '\n' from end of last column
y[0].extend(['File1','File2','File3','File4\n'])
for i in range(len(y)-1):
    y[i+1][-1] = y[i+1][-1][0:-1]
    y[i+1].extend([file1[i], file2[i], file3[i], file4[i] + '\n'])
    
f = open('EIL_site_lat_lon_files.csv','w')
for i in range(len(y)):
    for j in range(len(y[i])):
        if j != len(y[i]) - 1:
            f.write(y[i][j] + ',')
        else:
            f.write(y[i][j])
            
f.close()
