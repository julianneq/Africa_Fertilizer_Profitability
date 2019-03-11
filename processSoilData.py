import os
import arcpy as ap
from convertDBFtoCSV import convertDBFtoCSV

def processSoilData():
    '''Resamples soil data from 250m to 1000m using bilinear interpolation for continuous data, nearest neighbor interpolation
    for discrete data. Then finds the soil values of the cell each site is inside.'''

    #find the working directory
    workingDir = os.getcwd()

    # Overwrite pre-existing files
    #ap.env.overwriteOutput = True

    #Turn on Spatial Statistics package and define points at which soil characteristics will be calculated
    ap.CheckOutExtension("Spatial")
    in_point_features = workingDir + "\\Sites\\Sites.shp"

    #Define the projection and change the working directory to the directory with all of the soil data folders
    latLongRef = "Coordinate Systems/Geographic Coordinate Systems/World/WGS 1984.prj"
    os.chdir(workingDir + "\\Soils")
    tiffs250 = [f for f in os.listdir(os.getcwd()) if f[-8::] == '250m.tif']

    cell_size = "1000 1000"
    for i in range(len(tiffs250)):
        #Resample the data from 250m resolution to 1000m resolution
        in_raster = os.getcwd() + "\\" + tiffs250[i]
        out_raster = os.getcwd() + "\\" + tiffs250[i][0:-8] + "1000m.tif"
        if in_raster == os.getcwd() + "\\af_DRAINFAO_T__M_250m.tif":
            resampling_type = "NEAREST"
        else:
            resampling_type = "BILINEAR"

        ap.Resample_management(in_raster, out_raster, cell_size, resampling_type)
            
    tiffs1000 = [f for f in os.listdir(os.getcwd()) if f[-9::] == '1000m.tif']
    for j in range(len(tiffs1000)):
        in_raster = os.getcwd() + "\\" + tiffs1000[j]
        out_point_features = os.getcwd() + "\\" + tiffs1000[j][0:-4] + "_Sites.shp"
        #Calculate Zonal Statistics of soil data at Sites
        ap.sa.ExtractValuesToPoints(in_point_features, in_raster, out_point_features)
 
    #Convert the DBFs with all the soil data to CSVs
    #Change the working directory to the directory with all the soil data folders
    filelist = os.listdir(os.getcwd())
    DBFs = []
    for j in range(len(filelist)):
        name = filelist[j]
        if name[-4::] == '.dbf':
            DBFs.append(name)

    #Convert the DBF to a CSV
    for j in range(len(DBFs)):
        convertDBFtoCSV(os.getcwd() + "\\" + DBFs[j])
            
    return None
