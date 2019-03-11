import numpy as np

database = np.loadtxt('AfricaDatabase.csv',delimiter=',',skiprows=1)
badrows = [8080, 9519, 9605, 9606, 10091, 10092, 10916, 11135, 11577, 11670, 11671, 11768, 11769, 12065, 12163, 12268, 13068, 13320, 13446, 13573, 14574, 14622, 14805, 15061, 15132]
database = np.delete(database, badrows, 0)
np.savetxt('AfricaDatabase2.csv',database,delimiter=',')