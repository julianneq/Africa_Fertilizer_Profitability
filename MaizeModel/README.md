The files in this folder compile the complete estimation dataset and simulation datasets (in stata), using as inputs the climate and soils data layers. After preparing the data for estimation, the R files estimate the causal forest model (in R). And then the yield response and profitability results are analyzed, and the simulations are completed, in Stata.

To run these files, start by running: 
MaizeModel_master.do
[Note: within that file, the user should set the directory location of the raw climate data.]

Then the user should run the following do files from MaizeModel_master.do:
MaizeModel_1_CreateTrialData.do
MaizeModel_2_CreateSimData.do
MaizeModel_3_PrepDataForR.do
MaizeModel_3.5_RegressionEst.do

Before the remaining do files can be run, these tw R scripts should be run:
MaizeModel_Causal forest.R
MaizeModel_randomforest.R

Then the final two do files can be run from MaizeModel_master.do:
MaizeModel_4_ResultsTrialModel.do
MaizeModel_5_ResultsSimulations.do

