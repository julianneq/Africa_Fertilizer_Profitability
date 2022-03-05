set more off
cap log close

clear mata
clear matrix
mata: mata clear


** Folder locations

gl DROPBOX ""
gl DATA "Data"   /*Data from the Africa fertilizer trial and matched geo data*/
gl TEMP "Temp"
gl RESULTS "Results"
gl FIGS "Figures"
gl CLIMDATA "Y:/Data/FertInvest" 	/* NOTE: Need to set this location - data can be downloaded from the Weather folder here: https://drive.google.com/drive/folders/1KaeIACKpuA1nyuimD1KhtolZ2G1_lszA */


** Trial dataset inclusion
gl toinclude_trial "!(mgmt_lowph==1 | mgmt_strk==1 | mgmt_dr==1)" 


** Price assumptions

gl ureaPrice=260*1.3333 /*USD/MT World Price from Andrew; markup from Chemonics2007 in refs*/
gl maizePrice=170*0.75 /*USD/MT World Price from Andrew; margin from Yamano2011 in refs*/
gl fertilizerTreatmentVolume=0.125 /*MT fert per ha*/
gl int_rate=0.14 /*interest rate % per ag season*/
gl simreps=1000 
gl bootreps=1
set seed 123456
gl TM = 30 /* max temp above which GDD doesn't benefit*/
gl TB = 10 /* base temp*/


** Variables

gl idvar 		 		sitecode management date_plant vargroup site_year_id toinclude country
gl fertvar 				fert 
gl depvarlog 			logyield
gl depvarlev 			levyield
gl xvars_dum 		 	poordrain wortmann var_hy /*omitted open populated varieties*/
gl xvars_cont_clim 		temp_p1 temp_p2 temp_p3 precip_p1 precip_p2 precip_p3 
gl xvars_cont_fixed 	soilcec soilph acidity soilom soiln claypct siltpct bulkdens elevm
gl xvars_cont			${xvars_cont_clim} ${xvars_cont_fixed}
gl clustvar 			site_year_id
gl wvar 				ps_pro


do "MaizeModel_1_CreateTrialData.do"  
do "MaizeModel_2_CreateSimData.do"
do "MaizeModel_3_PrepDataForR.do"
do "MaizeModel_3.5_RegressionEst.do"

** --> Run the "Causal forest.R" file in R before running MaizeModel_4_ResultsTrialModel.do and MaizeModel_5_ResultsSimulations.do

*do "MaizeModel_4_ResultsTrialModel.do"  
*do "MaizeModel_5_ResultsSimulations.do"


