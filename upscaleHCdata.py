import os
import arcpy as ap
from convertDBFtoCSV import convertDBFtoCSV

def upscaleHCdata():
    #set the working directory
    workingDir = os.getcwd()

    # Overwrite pre-existing files
    #ap.env.overwriteOutput = True

    #Turn on Spatial Statistics package and define points at which soil characteristics will be calculated
    ap.CheckOutExtension("Spatial")
    in_zone_data = workingDir + "/Grid/AfricaGrid.shp"
    zone_field = "FID"

    #Define the projection and change the working directory to the directory with all of the soil data folders
    latLongRef = "Coordinate Systems/Geographic Coordinate Systems/World/WGS 1984.prj"
    os.chdir(workingDir + "/HarvestChoice")
    HCdir = os.getcwd()
    directories = [f for f in os.listdir(HCdir) if os.path.isfile(f) == False]
    variables = ['HarvestedArea','Production']

    for i in range(len(directories)):
        for variable in variables:
            #Find all the tiffs with soil data in each soil characteristic folder
            os.chdir(HCdir + "\\" + directories[i] + "\\" + variable)
            filelist = os.listdir(os.getcwd())
            tiffs = [j for j in filelist if j[-4::]=='.tif']
                
            for j in range(len(tiffs)):
                in_value_raster = os.getcwd() + "\\" + tiffs[j]

                #Calculate Zonal Statistics of soil data at AggLevel
                out_table = os.getcwd() + "\\AfricaGrid_" + tiffs[j][0:-4] + ".dbf"
                ap.sa.ZonalStatisticsAsTable(in_zone_data, zone_field, in_value_raster, out_table)

                # convert dbf to csv
                convertDBFtoCSV(out_table)

    return None