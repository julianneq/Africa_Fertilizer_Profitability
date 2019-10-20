import arcpy as ap
import os

# load spatial analysis and geostatistical analysis toolboxes
ap.CheckOutExtension("Spatial")

commodities = ["maiz","urea"] # call maize "maiz" because limits on word length
semivariograms = ["Spherical","Gaussian"] # Krige maize with Spherical and urea with Gaussian
estimates = ["beta","std_error"] # columns of shapefile attribute table to Krige
abbrevs = ["b","se"] # abbreviations of above estimates
regions = ["A","E","W"] # regions over which the Kriging is performed; A=Africa, E=East, W=West
lagsize = [["0.172616","0.086841","0.052430"], # lags for maize Kriged over A, E and W
	["0.144033","0.052656","0.045133"]] # lags for urea Kriged over A, E and W

# specify raster grid cells for interpolation - use aez raster
out_raster_cells = os.getcwd() + "/AEZ/aezraster"

# Krige betas and std_erorrs for maize and urea
# get mean and variance of Kriged predictions
for i, commodity in enumerate(commodities):
	for j, estimate in enumerate(estimates):
		for k, region in enumerate(regions):
			ap.gp.Kriging_sa(os.getcwd() + "/Prices/" + commodity + "_price_pts_" + region + ".shp", # in_point_features
				estimate, # z_field
				os.getcwd() + "/Prices/" + commodity + "_" + abbrevs[j] + "_" + region, # output raster name
				semivariograms[i] + " " + lagsize[i][k], # semivariogram type and lag
				out_raster_cells, # resolution of output raster
				"VARIABLE 12", # variable search radius; 12 points
				os.getcwd() + "/Prices/" + commodity + "_" + abbrevs[j] + "_" + region + "_var") # output variance of prediction raster name


# convert Kriged estimates from raster to point shapefile
for commodity in commodities:
	for j, estimate in enumerate(estimates):
		for region in regions:
			# mean predictions of betas and standard errors
			ap.RasterToPoint_conversion(os.getcwd() + "/Prices/" + commodity + "_" + abbrevs[j] + "_" + region, # input raster
				os.getcwd() + "/Prices/" + commodity + "_" + abbrevs[j] + "_" + region + ".shp", # output shapefile
				"Value") # feature

			# variance of predictions of betas and standard errors
			ap.RasterToPoint_conversion(os.getcwd() + "/Prices/" + commodity + "_" + abbrevs[j] + "_" + region + "_var", # input raster
				os.getcwd() + "/Prices/" + commodity + "_" + abbrevs[j] + "_" + region + "_var.shp", # output shapefile
				"Value") # feature

# perform a spatial join of point shapefiles and Africa grid cells
# mimicked code from here: https://pro.arcgis.com/en/pro-app/tool-reference/analysis/spatial-join.htm
AfricaGridCells = os.getcwd() + "/Grid/AfricaGridJoin"

isvar = ["","_var"] # ending if estimate is mean of variance of prediction

count = 1 # iterate counter each time a field is spatially joined to AfricaGridCells
for commodity in commodities:
	for j, estimate in enumerate(estimates):
		for region in regions:
			for n in range(len(isvar)):
				# specify target (grid cells), join (Kriged estimates) and output (grid w/ Kriged estimates) feature classes
				if count == 1:
					target_features = AfricaGridCells + ".shp"
				else:
					target_features = AfricaGridCells + str(count-1) + ".shp"
				join_features = os.getcwd() + "/Prices/" + commodity + "_" + abbrevs[j] + "_" + region + isvar[n] + ".shp"
				out_feature_class = AfricaGridCells + str(count) + ".shp"

				# create fieldmappings for output (mean of Kriged estimates within Africa grid cell)
				fieldmappings = ap.FieldMappings()
				fieldmappings.addTable(target_features)
				fieldmappings.addTable(join_features)
				Kriged_estimate = fieldmappings.findFieldMapIndex("grid_code")
				fieldmap = fieldmappings.getFieldMap(Kriged_estimate)

				# Get the output field's properties as a field object
				field = fieldmap.outputField

				# Rename the field and pass the updated field object back into the field map
				field.name = commodity + "_" + abbrevs[j] + "_" + region + isvar[n]
				field.aliasName = commodity + "_" + abbrevs[j] + "_" + region + isvar[n]
				fieldmap.outputField = field

				# Set the merge rule to mean and then replace the old fieldmap in the mappings object with the updated one
				fieldmap.mergeRule = "mean"
				fieldmappings.replaceFieldMap(Kriged_estimate, fieldmap)

				# Delete fields that are no longer applicable, as only the first value will be used by default
				x = fieldmappings.findFieldMapIndex("pointid")
				fieldmappings.removeFieldMap(x)

				# join target_features and join_features to create joined out_feature_class
				ap.SpatialJoin_analysis(target_features,
					join_features,
					out_feature_class,
					"#", # {join_operation}, 
					"#", #{join_type}, 
					fieldmappings) # {field_mapping})

				# fix this code - this deletes the data in the shapefile, not the shapefile
				# deleted manually for now
				#if count == 1: # delete former join shapefile
				#	ap.DeleteFeatures_management(AfricaGridCells + ".shp")
				#else:
				#	ap.DeleteFeatures_management(AfricaGridCells + str(count-1) + ".shp")

				count += 1 # iterate count

# rename last joined database - need to fix; deleted manually for now
#ap.CopyFeatures_management(AfricaGridCells + str(count-1) + ".shp", AfricaGridCells + ".shp")
#ap.DeleteFeatures_management(AfricaGridCells + str(count-1) + ".shp")

# need to add following steps to code (did manually for now):

# 1. remove unwanted columns from joined dataset (# of cells averaged in the spatial join and their average ID)
# 2. replace 0s with -9999s
# 3. convert shapefile's .dbf to .csv
