***************************************************************************************
***************************************************************************************
******************************** RESULTS TRIAL MODEL **********************************
***************************************************************************************
***************************************************************************************

** This file takes the Random Forest results (from R) and analyzes them 

set more off
cap log close

clear mata
clear matrix
mata: mata clear

loc critval			1.96 /*1.96 for 95% cis 1.65 for 90% cis*/
loc bwval = 1/30

log using "${TEMP}/MaizeModel_ResultsTrialModel.log", replace


gl allvarslev1 ${xvars_cont} 
gl allvarslev2 ${xvars_dum}


*******************
** Basic results **
*******************

** Load RF data

use "RandomForest/Results_TrialData_rf.dta", clear

rename trialData_id id
rename Y y_act
rename predictions y_pred
rename variance_estimates y_pred_act_var
rename debiased_error y_pred_debiasederror
rename excess_error y_pred_excesserror
rename predictions_1 y_pred_f0 
rename variance_estimates_1 y_pred_f0_var
rename predictions_2 y_pred_f1
rename variance_estimates_2 y_pred_f1_var
foreach type in act f0 f1 {
	g y_pred_`type'_se = y_pred_`type'_var^0.5
	}

merge 1:1 id using "RandomForest/TrialData_act_lev.dta", nogen keep(1 3) keepusing($fertvar $depvarlev $xvars_dum $xvars_cont $clustvar $wvar)

tempfile TrialData_rf_validation
save `TrialData_rf_validation', replace

save "${TEMP}/TrialData_rf_validation.dta", replace


** Load CF data

use "RandomForest/Results_TrialData_cf.dta", clear

rename trialData_id id
rename predictions y_pred_dif
rename debiased_error y_pred_dif_debiasederror
rename excess_error y_pred_dif_excesserror
rename variance_estimates y_pred_dif_var
g y_pred_dif_se = sqrt(y_pred_dif_var)

merge 1:1 id using "RandomForest/TrialData_act_lev.dta", nogen keep(1 3) keepusing($fertvar $depvarlev $xvars_dum $xvars_cont $clustvar $wvar)

tempfile TrialData_cf_results
save `TrialData_cf_results', replace

save "${TEMP}/TrialData_cf_results.dta", replace


** Compare predicted yields to actual yields

use `TrialData_rf_validation', replace

regress y_pred y_act, noconstant 

g y_pred_err = y_act - y_pred
g y_pred_err2 = y_pred_err^2

sktest y_pred_err

bysort $fertvar: su y_pred_err y_pred_err2

su y_pred_err2
loc sse = r(mean)
loc rmse = `sse'^.5
display "RMSE = `rmse'" 
su y_act


** Depict variable importance

* First for RF

use "RandomForest/Varimp_RF", clear

g varid=_n
g varname = ""
replace varname = "$fertvar" in 1
loc i = 2

forvalues j = 1/2 {
	loc K = wordcount("${allvarslev`j'}") 
	forvalues k = 1/`K' {
		loc varname: word `k' of ${allvarslev`j'}
		replace varname = "`varname'" in `i'
		loc i = `i' + 1
		}
	}

gsort + varimp
g sortorder = _n

keep varid varname varimp sortorder

replace sortorder  = sortorder * 12
g min = -.1
g axis = 0
graph twoway (bar varimp sortorder, horizontal barwidth(10) ylabel(none) ytitle("") yscale(lstyle(none)) xlabel(0 (0.05) .3) xtitle("Variable Importance") legend(off)) (scatter sortorder axis,  m(i) mlabel(varname) mlabposition(9) mlabtextstyle(medlarge)) (scatter sortorder min,  m(i)) 
graph display, ysize(11) scale(1.2)

graph export "${FIGS}/VariableImportance_RF.pdf", replace

tempfile varimp_rf
save `varimp_rf', replace


* Then for CF
use  "RandomForest/Varimp_CF", clear

g varid=_n
g varname = ""
loc i = 1

forvalues j = 1/2 {
	loc K = wordcount("${allvarslev`j'}") 
	forvalues k = 1/`K' {
		loc varname: word `k' of ${allvarslev`j'}
		replace varname = "`varname'" in `i'
		loc i = `i' + 1
		}
	}

gsort + varimp
g sortorder = _n

keep varid varname varimp sortorder

replace sortorder  = sortorder * 12
g min = -.06
g axis = 0
graph twoway (bar varimp sortorder, horizontal barwidth(10) ylabel(none) ytitle("") yscale(lstyle(none)) xlabel(0 (0.025) .125) xtitle("Variable Importance") legend(off)) (scatter sortorder axis,  m(i) mlabel(varname) mlabposition(9) mlabtextstyle(medlarge)) (scatter sortorder min,  m(i)) 
graph display, ysize(11) scale(1.2)

graph export "${FIGS}/VariableImportance_CF.pdf", replace

tempfile varimp_cf
save `varimp_cf', replace


** Compare distributions of predicted yields (act fert use)

use `TrialData_rf_validation', clear

* Test distribution diff - predicted yields by actual fert use
ksmirnov y_pred, by($fertvar) /*f1 is larger yield*/
ttest y_pred, by($fertvar) unequal

* Test distribution diff - actual yields by actual fert use
ksmirnov y_act, by($fertvar) /*f1 is larger yield*/
ttest y_act, by($fertvar) unequal

* Test distribution diff - side by side predicted yields simulated fert use
ttest y_pred_f0 = y_pred_f1
reshape long y_pred_f, i(id) j(f2)
ksmirnov y_pred_f, by(f2) /*f1 is larger yields)*/

* Graphically compare predicted yield distributions by actual fert use
use `TrialData_rf_validation', clear
twoway kdensity y_pred if fert==0 [aw=$wvar] || kdensity y_pred if fert==1 [aw=$wvar], xtitle("Predicted yields (t/ha)") ytitle("") legend(row(1) order(1 "No Fertilizer" 2 "Fertilizer" ))
graph export "${FIGS}/Distributions_ypred_act_RF.pdf", replace

* Graphically compare f0 and f1 predicted yield distributions assigned fert use
twoway kdensity y_pred_f0 [aw=$wvar] || kdensity y_pred_f1 [aw=$wvar], xtitle("Predicted yields (t/ha)") legend(row(1) order(1 "No Fertilizer" 2 "Fertilizer" ))
graph export "${FIGS}/Distributions_ypred_f0f1_RF.pdf", replace

* Graphically compare scatter of f1 vs f0 predicted yield assigned fert use
su y_pred
loc min = r(min)
loc max = r(max)
g min = runiform(`min',`max')
sort min 
graph twoway (scatter y_pred_f1 y_pred_f0) (lpoly y_pred_f1 y_pred_f0 [aw=$wvar]) (line min min), xtitle("F=0 (predicted yields, kg/ha)") ytitle("F=1 (predicted yields, kg/ha)") xlabel(2(1)6) ylabel(2(1)6) legend(off)
graph export "${FIGS}/Scatter_ypred_f0f1_RF.pdf", replace

* Graphically compare scatter of predicted fert dif 
g y_pred_dif = y_pred_f1 - y_pred_f0
sort min 
graph twoway (scatter y_pred_dif y_pred_f0) (lpoly y_pred_dif y_pred_f0 [aw=$wvar]), yline(0) xtitle("F=0 (predicted yields, kg/ha)") ytitle("Predicted fert response (kg/ha)") xlabel(2(1)4.5) ylabel(-0.5(0.5)4.5) legend(off)
graph export "${FIGS}/Scatter_ydif_f0_RF.pdf", replace


** Compare predicted yields to actual

*collapse (mean) y_pred y_act, by(sitecode fert)
cap drop min
su y_pred
loc min = r(min)
loc max = r(max)
g min = runiform(`min',`max')
sort min 
graph twoway (scatter y_act y_pred if fert==1) (scatter y_act y_pred if fert==0) (lpoly y_act y_pred if fert==1) (lpoly y_act y_pred if fert==0) (line min min), legend(order(1 "Y pred vs y act (f=0)" 2 "Y pred vs y act (f=1)")) xscale(range(0(5)15)) yscale(range(0(5)20)) xlabel(0(5)15)  ylabel(0(5)20) xtitle("Predicted yields, kg/ha") ytitle("Actual yields, kg/ha")
graph export "${FIGS}/Scatter_ypred_yact_RF.pdf", replace

regress y_act y_pred
regress y_act y_pred, noconstant
regress y_act y_pred [aw=${wvar}]
regress y_act y_pred [aw=${wvar}], noconstant 
regress y_act y_pred [aw=${wvar}], vce(cluster site_year_id)
regress y_act y_pred [aw=${wvar}], noconstant vce(cluster site_year_id)
correlate y_act y_pred
correlate y_act y_pred [aw=${wvar}]


** Compare distributions of predicted yields (sim fert 0 and fert 1)

use `TrialData_rf_validation', clear

* diff between two rvs has a var of var(a) + var(b)
g fertdif_se = (y_pred_f0_var + y_pred_f1_var)^0.5

tempfile TrialData_validation
save `TrialData_rf_validation', replace 
save "${TEMP}/TrialData_rf_validation.dta", replace 


** CF Graphs 
* Graphically compare predicted yield difference distributions by actual fert use
use `TrialData_cf_results', clear


twoway kdensity y_pred_dif [aw=$wvar], xtitle("Predicted yield response (t/ha)") ytitle("") legend(off)
graph export "${FIGS}/Distributions_y_pred_dif_CF.pdf", replace

twoway kdensity y_pred_dif if fert==0 [aw=$wvar] || kdensity y_pred_dif if fert==1 [aw=$wvar], xtitle("Predicted yield response (t/ha)") ytitle("") legend(row(1) order(1 "Fertilizer used" 2 "Fertilizer not used" ))
graph export "${FIGS}/Distributions_y_pred_dif_by_fert_use_CF.pdf", replace


***************************************************************************
** Heterogeneity results - simulate climate variation across trial sites **
***************************************************************************

** Create graphical descriptions for each continous xvar and dummy (lpoly) - for CF

use `TrialData_cf_results', clear
merge 1:1 id using "${TEMP}/TrialData_orig", nogen keep(1 3)

foreach xvar in $xvars_cont {
	qui su `xvar'_orig
	loc `xvar'_orig_mean = r(mean)
	loc `xvar'_orig_sd = r(sd)
	}
	
g temp_ave_orig = (2*temp_p1_orig + temp_p2_orig + 2*temp_p3_orig)/5
la var temp_ave_orig "Mean growing season temperature (degrees C)"
egen precip_tot_orig = rowtotal(precip_p1_orig precip_p2_orig precip_p3_orig)
la var precip_tot_orig "Total growing season precipitation (mm)"

g y_pred_dif_lo = y_pred_dif - `critval'*y_pred_dif_se 
g y_pred_dif_hi = y_pred_dif + `critval'*y_pred_dif_se 

tempfile graphdata
save `graphdata', replace


foreach var in temp_ave precip_tot $xvars_cont {
	use `graphdata', clear
		
	_pctile `var'_orig, p(2.5)
	loc xmin = r(r1)
	_pctile `var'_orig, p(97.5)
	loc xmax = r(r1)
	loc bw = (`xmax'-`xmin')*`bwval'

	kdensity `var'_orig if `var'_orig >= `xmin' & `var'_orig <= `xmax', nograph generate(atx_`var'_dens predd_`var') 

	loc varlabel : var label `var'
	parse "`varlabel'", parse("(")
	loc `var'_title `1'
		
	foreach type in dif {
		lpoly y_pred_`type' `var'_orig [aw=${wvar}], nograph degree(0) bw(`bw') at(atx_`var'_dens) generate(yp_`type')
		lpoly y_pred_`type'_lo `var'_orig [aw=${wvar}], nograph degree(0) bw(`bw') at(atx_`var'_dens) generate(yplo_`type')
		lpoly y_pred_`type'_hi `var'_orig [aw=${wvar}], nograph degree(0) bw(`bw') at(atx_`var'_dens) generate(yphi_`type')
		
		foreach pred in yp yplo yphi {
			replace `pred'_`type' = . if predd_`var' == 0  
			}
		}

	save "${TEMP}/AveFertResponse_CF_`var'.dta", replace
	}

	
loc ysc_yielddif	"-2 3.8"
loc ylab_yielddif	"0(1)3"
loc yaxlab 			"(t/ha)"

foreach var in temp_ave precip_tot $xvars_cont {
	use "${TEMP}/AveFertResponse_CF_`var'.dta", clear
	
	_pctile `var'_orig, p(1)
	loc xmin = r(r1)
	_pctile `var'_orig, p(99)
	loc xmax = r(r1)
	
	qui su predd_`var'
	loc y2scalemax = r(max) * 5

	loc varlabel : var label `var'
	parse "`varlabel'", parse("(")
	loc vartitle `1'

	graph twoway 	(line yp_dif atx_`var'_dens, lpattern(solid) lcolor(black) lwidth(medthick)) ///
					(rarea yplo_dif yphi_dif atx_`var'_dens, lwidth(none) fcolor(gray%20)) ///
					(line predd_`var' atx_`var'_dens, lcolor(gray) lwidth(vthin) lpattern(solid) yaxis(2)) ///
						||, ///
						xsc(r(`xmin' `xmax')) ///
						yscale(axis(1) range(`ysc_yielddif')) ylabel(`ylab_yielddif', nogrid axis(1))  ///
						yscale(axis(2) noline range(0 `y2scalemax')) ylabel(none, nogrid axis(2)) ///
						title(""/*"`vartitle'"*/) ///
						xtitle("`vartitle'"/*"`varlabel'"*/) ///
						ytitle("Fert response (t/ha)"/* "Predicted fertilizer response `yaxlab_`level''"*/, axis(1)) ///
						ytitle(""/*"Density:" "`vartitle'"*/, placement(se) axis(2) justification(left)) ///
						yline(0, lpattern(dash) lcolor(black) axis(1)) ///
						legend(off) ///
						saving("${TEMP}/ydif_`var'", replace)
	graph export "${FIGS}/AveFertResponse_CF_yielddif_`var'.pdf", replace
	}

graph combine "${TEMP}/ydif_temp_p1.gph" "${TEMP}/ydif_temp_p2.gph" "${TEMP}/ydif_temp_p3.gph" "${TEMP}/ydif_precip_p1.gph" "${TEMP}/ydif_precip_p2.gph" "${TEMP}/ydif_precip_p3.gph" /*
*/"${TEMP}/ydif_soilcec.gph" "${TEMP}/ydif_soilph.gph" "${TEMP}/ydif_acidity.gph" "${TEMP}/ydif_soilom.gph" "${TEMP}/ydif_soiln.gph" "${TEMP}/ydif_claypct.gph" "${TEMP}/ydif_siltpct.gph" "${TEMP}/ydif_bulkdens.gph" "${TEMP}/ydif_elevm.gph", cols(3) xsize(6) ysize(6)
graph export "${FIGS}/AveFertResponse_CF_yielddif_allvars.pdf", replace

graph combine "${TEMP}/ydif_temp_ave.gph" "${TEMP}/ydif_precip_tot.gph" /*
*/"${TEMP}/ydif_soilcec.gph" "${TEMP}/ydif_soilph.gph" "${TEMP}/ydif_acidity.gph" "${TEMP}/ydif_soilom.gph" "${TEMP}/ydif_soiln.gph" "${TEMP}/ydif_claypct.gph" "${TEMP}/ydif_siltpct.gph" "${TEMP}/ydif_bulkdens.gph" "${TEMP}/ydif_elevm.gph", cols(3) xsize(6) ysize(6)
graph export "${FIGS}/AveFertResponse_CF_yielddif_shortlistvars.pdf", replace


******************************
** Isolate Climate in Sites **
******************************	

** First establish limits of historic climate record

use "${TEMP}/AfricaSites_bt_tempprecip_all.dta", clear 
su simyr
loc yrmax = r(max)
loc yrmin = r(min)
loc nyears = 1 + `yrmax' - `yrmin'


** RF model - Isolate climate characteristics -- predict fert response with each site's historic climate and look at the distributions

use "RandomForest/Results_TrialSim_rf.dta", clear

rename trialDataSimClimate_sitecode sitecode
rename trialDataSimClimate_simyr simyr
rename trialDataSimClimate_fert $fertvar
rename predictions y_pred_f
rename variance_estimates y_pred_se_f
replace y_pred_se_f = sqrt(y_pred_se_f)

sort sitecode simyr $fertvar
merge 1:1 sitecode simyr $fertvar using "RandomForest/SimData_trialsim_lev", nogen keep(3) keepusing($xvars_cont $xvars_dum)
merge m:1 sitecode using "${TEMP}/psweights_site.dta", nogen keep(1 3) 

reshape wide y_pred_f y_pred_se_f, i(sitecode simyr) j(fert)
rename y_pred_se_f0 y_pred_f0_se
rename y_pred_se_f1 y_pred_f1_se

g y_pred_dif = y_pred_f1 - y_pred_f0
g y_pred_dif_se = (y_pred_f1_se^2 + y_pred_f0_se^2)^0.5

foreach type in dif f0 f1 {
	g y_pred_`type'_lo = y_pred_`type' - `critval'*y_pred_`type'_se
	g y_pred_`type'_hi = y_pred_`type' + `critval'*y_pred_`type'_se
	}
	
foreach var in temp_p1 temp_p2 temp_p3 precip_p1 precip_p2 precip_p3 soilcec soilph claypct siltpct bulkdens acidity soilom soiln elevm {
	display "g `var'_orig = `var' * ``var'_orig_sd' + ``var'_orig_mean'"
	g `var'_orig = `var' * ``var'_orig_sd' + ``var'_orig_mean'
	su `var' `var'_orig 
	}
	
g temp_ave_orig = (2*temp_p1_orig + temp_p2_orig + 2*temp_p3_orig)/5
egen precip_tot_orig = rowtotal(precip_p1_orig precip_p2_orig precip_p3_orig)


tempfile Trialsim_rf
save `Trialsim_rf', replace

g count = 1
collapse (sum) count, by(sitecode)
su count 
loc size = r(mean)
loc size_sd = r(sd)
assert `size_sd' == 0
drop count

expand $simreps
g count = 1
bysort sitecode : g iter = sum(count) 
drop count

generate simyr = floor(`yrmin' + `nyears'*runiform())
su simyr

merge m:1 sitecode simyr using `Trialsim_rf', nogen keep(1 3)

sort sitecode iter
isid sitecode iter

save "${TEMP}/TrialSim_rf_processed.dta", replace


** Isolate climate characteristics -- predict fert response with each site's historic climate and look at the distributions

use "RandomForest/Results_TrialSim_cf.dta", clear

rename trialDataSimClimate_sitecode sitecode
rename trialDataSimClimate_simyr simyr
rename predictions y_pred_dif
rename variance_estimates y_pred_dif_se
replace y_pred_dif_se = sqrt(y_pred_dif_se)

sort sitecode simyr
merge 1:1 sitecode simyr using "RandomForest/SimData_trialsim_lev_cf", nogen keep(3) keepusing($xvars_cont $xvars_dum)
merge m:1 sitecode using "${TEMP}/psweights_site.dta", nogen keep(1 3) 

foreach type in dif {
	g y_pred_`type'_lo = y_pred_`type' - `critval'*y_pred_`type'_se
	g y_pred_`type'_hi = y_pred_`type' + `critval'*y_pred_`type'_se
	}

foreach var in temp_p1 temp_p2 temp_p3 precip_p1 precip_p2 precip_p3 soilcec soilph claypct siltpct bulkdens acidity soilom soiln elevm {
	display "g `var'_orig = `var' * ``var'_orig_sd' + ``var'_orig_mean'"
	g `var'_orig = `var' * ``var'_orig_sd' + ``var'_orig_mean'
	su `var' `var'_orig 
	}
	
g temp_ave_orig = (2*temp_p1_orig + temp_p2_orig + 2*temp_p3_orig)/5
egen precip_tot_orig = rowtotal(precip_p1_orig precip_p2_orig precip_p3_orig)


tempfile Trialsim
save `Trialsim', replace

g count = 1
collapse (sum) count, by(sitecode)
su count 
loc size = r(mean)
loc size_sd = r(sd)
assert `size_sd' == 0
drop count

expand $simreps
g count = 1
bysort sitecode : g iter = sum(count) 
drop count

generate simyr = floor(`yrmin' + `nyears'*runiform())
su simyr

merge m:1 sitecode simyr using `Trialsim', nogen keep(1 3)

sort sitecode iter
isid sitecode iter

save "${TEMP}/TrialSim_cf_processed.dta", replace

su simyr 

loc ysc_yielddif	"-2 3.5"
loc ylab_yielddif	"0(1)3"
loc ysc_yieldcomp	"-1 6"
loc ylab_yieldcomp	"2(1)6"
loc yaxlab 			"(t/ha)"
loc yaxlab 			"(t/ha)"

foreach var in temp_ave temp_p1 temp_p2 temp_p3 precip_tot precip_p1 precip_p2 precip_p3 soilcec soilph claypct siltpct bulkdens acidity soilom soiln elevm {
	use "${TEMP}/TrialSim_cf_processed.dta", clear
	
	_pctile `var'_orig, p(2.5)
	loc xmin = r(r1)
	_pctile `var'_orig, p(97.5)
	loc xmax = r(r1)
	loc bw = (`xmax'-`xmin')*`bwval'


	kdensity `var'_orig if `var'_orig >= `xmin' & `var'_orig <= `xmax', nograph generate(atx_`var'_dens predd_`var') 		
	foreach type in dif {
		lpoly y_pred_`type' `var'_orig [aw=${wvar}], nograph bw(`bw') degree(0) at(atx_`var'_dens) generate(yp_`type')
		lpoly y_pred_`type'_lo `var'_orig [aw=${wvar}], nograph bw(`bw') degree(0) at(atx_`var'_dens) generate(yplo_`type')
		lpoly y_pred_`type'_hi `var'_orig [aw=${wvar}], nograph bw(`bw') degree(0) at(atx_`var'_dens) generate(yphi_`type')
		
		foreach pred in yp yplo yphi {
		    replace `pred'_`type' = . if atx_`var'_dens == . 
			}
		}

	save "${TEMP}/AveFertResponse_TrialDataSimClim_CF_`var'.dta", replace
	}
	

foreach var in temp_ave temp_p1 temp_p2 temp_p3 precip_tot precip_p1 precip_p2 precip_p3 soilcec soilph claypct siltpct bulkdens acidity soilom soiln elevm {
	use "${TEMP}/AveFertResponse_TrialDataSimClim_CF_`var'.dta", clear
	
	_pctile `var'_orig, p(2.5)
	loc xmin = r(r1)
	_pctile `var'_orig, p(97.5)
	loc xmax = r(r1)
	
	qui su predd_`var'
	loc y2scalemax = r(max) * 5

	graph twoway 	(line yp_dif atx_`var'_dens, lpattern(solid) lcolor(black) lwidth(medthick)) ///
					(rarea yplo_dif yphi_dif atx_`var'_dens, lwidth(none) fcolor(gray%20)) ///
					(line predd_`var' atx_`var'_dens, lcolor(gray) lwidth(vthin) lpattern(solid) yaxis(2)) ///
						||, ///
						xsc(r(`xmin' `xmax')) ///
						yscale(axis(1) range(`ysc_yielddif')) ylabel(`ylab_yielddif', nogrid axis(1))  ///
						yscale(axis(2) noline range(0 `y2scalemax')) ylabel(none, nogrid axis(2)) ///
						title(""/*"`vartitle'"*/) ///
						xtitle("``var'_title'") ///
						ytitle("Fert response (t/ha)"/* "Predicted fertilizer response `yaxlab_`level''"*/, axis(1)) ///
						ytitle(""/*"Density:" "`vartitle'"*/, placement(se) axis(2) justification(left)) ///
						yline(0, lpattern(dash) lcolor(black) axis(1)) ///
						legend(off) ///
						saving("${TEMP}/TS_ydif_`var'", replace)
	graph export "${FIGS}/AveFertResponse_TrialDataSimClim_CF_yielddif_`var'.pdf", replace
	}	

graph combine "${TEMP}/TS_ydif_temp_p1.gph" "${TEMP}/TS_ydif_temp_p2.gph" "${TEMP}/TS_ydif_temp_p3.gph" "${TEMP}/TS_ydif_precip_p1.gph" "${TEMP}/TS_ydif_precip_p2.gph" "${TEMP}/TS_ydif_precip_p3.gph" /*
*/ "${TEMP}/TS_ydif_soilcec.gph" "${TEMP}/TS_ydif_soilph.gph" "${TEMP}/TS_ydif_claypct.gph" /*  
*/"${TEMP}/TS_ydif_siltpct.gph" "${TEMP}/TS_ydif_bulkdens.gph" "${TEMP}/TS_ydif_acidity.gph" /*
*/"${TEMP}/TS_ydif_soilom.gph" "${TEMP}/TS_ydif_soiln.gph" "${TEMP}/TS_ydif_elevm.gph" , cols(3) xsize(6) ysize(10)
graph export "${FIGS}/AveFertResponse_TrialDataSimClim_CF_yielddif_allvars.pdf", replace 

graph combine "${TEMP}/TS_ydif_temp_ave.gph" "${TEMP}/TS_ydif_precip_tot.gph" /*
*/ "${TEMP}/TS_ydif_soilcec.gph" "${TEMP}/TS_ydif_soilph.gph" "${TEMP}/TS_ydif_acidity.gph" "${TEMP}/TS_ydif_soilom.gph" "${TEMP}/TS_ydif_soiln.gph" "${TEMP}/TS_ydif_claypct.gph" "${TEMP}/TS_ydif_siltpct.gph" "${TEMP}/TS_ydif_bulkdens.gph" "${TEMP}/TS_ydif_elevm.gph", cols(3) xsize(6) ysize(6)
graph export "${FIGS}/AveFertResponse_TrialDataSimClim_CF_yielddif_shortlistvars.pdf", replace


** Calculate average fertilizer response - trial sites with bootstrapped clim

use "${TEMP}/TrialSim_cf_processed.dta", clear
g err = rnormal(0,y_pred_dif_se)
replace y_pred_dif = y_pred_dif + err
drop err

collapse (mean) y_pred_dif (sd) y_pred_dif_sd=y_pred_dif, by(sitecode simyr)

merge m:1 sitecode using "${TEMP}/psweights_site.dta", nogen keep(3) keepusing(${wvar}) 
su y_pred_dif y_pred_dif_sd [aw=${wvar}]

/* 
site simyr mean (trial sites sim climate)

CF DIFF: 
    Variable |     Obs      Weight        Mean   Std. Dev.       Min        Max
-------------+-----------------------------------------------------------------
  y_pred_dif |   4,260  1520.96215     1.49747   .4461754   .3054003    3.10816
y_pred_dif~d |   4,260  1520.96215    .5954449   .2120033   .1900496   1.721464
*/


use "${TEMP}/TrialSim_cf_processed.dta", clear
g err = rnormal(0,y_pred_dif_se)
replace y_pred_dif = y_pred_dif + err
drop err

collapse (mean) y_pred_dif (sd) y_pred_dif_sd=y_pred_dif, by(sitecode)

merge m:1 sitecode using "${TEMP}/psweights_site.dta", nogen keep(3) keepusing(${wvar}) 
su y_pred_dif y_pred_dif_sd [aw=${wvar}]


/*
site mean (trial sites sim climate)

CF DIFF: 
    Variable |     Obs      Weight        Mean   Std. Dev.       Min        Max
-------------+-----------------------------------------------------------------
  y_pred_dif |     142  50.6987382    1.500136   .3937717    .642675   2.588093
y_pred_dif~d |     142  50.6987382    .6410064   .1529585   .3240473    1.14884

*/


** Calculate average fertilizer response - trial dataset actual climate

use "${TEMP}/TrialData_cf_results.dta", clear
merge 1:1 id using "RandomForest/TrialData_act_lev.dta", keepusing(sitecode site_year_id)

collapse (mean) y_pred_dif (sd) y_pred_dif_sd=y_pred_dif , by(sitecode site_year_id)

merge m:1 sitecode using "${TEMP}/psweights_site.dta", nogen keep(3) keepusing(${wvar}) 
su y_pred_dif y_pred_dif_sd [aw=${wvar}]

/*
site year mean (trial sites actual climate)
  
CF DIFF:
    Variable |     Obs      Weight        Mean   Std. Dev.       Min        Max
-------------+-----------------------------------------------------------------
  y_pred_dif |     320  113.011659    1.615248    .664902   .2607676   3.650363
y_pred_dif~d |     320  113.011659    .0674186   .0905699          0   .4438037
*/

use "${TEMP}/TrialData_rf_validation.dta", clear
merge 1:1 id using "RandomForest/TrialData_act_lev.dta", keepusing(sitecode)

collapse (mean) fertdif y_pred_f0 y_pred_f1 (sd) fertdif_sd=fertdif y_pred_f0_sd=y_pred_f0 y_pred_f1_sd=y_pred_f1, by(sitecode)

merge m:1 sitecode using "${TEMP}/psweights_site.dta", nogen keep(3) keepusing(${wvar}) 
su fertdif fertdif_sd [aw=${wvar}]

su fertdif [aw=${wvar}]
loc rf_dif_mean = r(mean)
su fertdif_sd [aw=${wvar}]
loc rf_dif_sd = r(sd)
su y_pred_f0 [aw=${wvar}]
loc rf_f0_mean = r(mean)
su y_pred_f0_sd [aw=${wvar}]
loc rf_f0_sd = r(sd)
su y_pred_f1 [aw=${wvar}]
loc rf_f1_mean = r(mean)
su y_pred_f1_sd [aw=${wvar}]
loc rf_f1_sd = r(sd)

use "${TEMP}/TrialData_cf_results.dta", clear
merge 1:1 id using "RandomForest/TrialData_act_lev.dta", keepusing(sitecode)

collapse (mean) y_pred_dif (sd) y_pred_dif_sd=y_pred_dif , by(sitecode)

merge m:1 sitecode using "${TEMP}/psweights_site.dta", nogen keep(3) keepusing(${wvar}) 
su y_pred_dif y_pred_dif_sd [aw=${wvar}]

su y_pred_dif [aw=${wvar}]
loc cf_dif_mean = r(mean)
su y_pred_dif_sd [aw=${wvar}]
loc cf_dif_sd = r(sd)

/*
site mean (trial sites actual climate)

CF DIFF:
    Variable |     Obs      Weight        Mean   Std. Dev.       Min        Max
-------------+-----------------------------------------------------------------
 y_pred_dif |     142  50.6987382    1.488977   .5505216   .4464484   3.193114
y_pred_dif~d |     142  50.6987382    .1345267   .1300596          0    .550401
*/


** Now FGLS model
use "${TEMP}\Vselect_fgls_lev_fertresponse.dta", clear
merge 1:1 id using "RandomForest/TrialData_act_lev", nogen keep(1 3) keepusing(sitecode site_year_id)

su rmse_oos 
loc fgls_rmse = r(mean)

collapse (mean) y_pred_dif y_pred_f0 y_pred_f1 (sd) y_pred_dif_sd=y_pred_dif y_pred_f0_sd=y_pred_f0 y_pred_f1_sd=y_pred_f1, by(sitecode site_year_id)
merge m:1 sitecode using "${TEMP}/psweights_site.dta", nogen keep(3) keepusing(${wvar}) 

su y_pred_dif y_pred_dif_sd [aw=${wvar}]


use "${TEMP}\Vselect_fgls_lev_fertresponse.dta", clear
merge 1:1 id using "RandomForest/TrialData_act_lev", nogen keep(1 3) keepusing(sitecode site_year_id)

collapse (mean) y_pred_dif y_pred_f0 y_pred_f1 (sd) y_pred_dif_sd=y_pred_dif y_pred_f0_sd=y_pred_f0 y_pred_f1_sd=y_pred_f1, by(sitecode)
merge m:1 sitecode using "${TEMP}/psweights_site.dta", nogen keep(3) keepusing(${wvar}) 

su y_pred_dif y_pred_dif_sd [aw=${wvar}]

su y_pred_dif [aw=${wvar}]
loc fgls_dif_mean = r(mean)
su y_pred_dif_sd [aw=${wvar}]
loc fgls_dif_sd = r(sd)
su y_pred_f0 [aw=${wvar}]
loc fgls_f0_mean = r(mean)
su y_pred_f0_sd [aw=${wvar}]
loc fgls_f0_sd = r(sd)
su y_pred_f1 [aw=${wvar}]
loc fgls_f1_mean = r(mean)
su y_pred_f1_sd [aw=${wvar}]
loc fgls_f1_sd = r(sd)


*******************************
** Results comparison table  **
*******************************

mat fgls_rf_cf=J(4,6,.)
mat fgls_rf_cf[1,1] = `fgls_f0_mean'
mat fgls_rf_cf[1,2] = `fgls_f0_sd'
mat fgls_rf_cf[2,1] = `fgls_f1_mean'
mat fgls_rf_cf[2,2] = `fgls_f1_sd'
mat fgls_rf_cf[3,1] = `fgls_dif_mean'
mat fgls_rf_cf[3,2] = `fgls_dif_sd'
mat fgls_rf_cf[4,1] = `fgls_rmse'

mat fgls_rf_cf[1,3] = `rf_f0_mean'
mat fgls_rf_cf[1,4] = `rf_f0_sd'
mat fgls_rf_cf[2,3] = `rf_f1_mean'
mat fgls_rf_cf[2,4] = `rf_f1_sd'
mat fgls_rf_cf[3,3] = `rf_dif_mean'
mat fgls_rf_cf[3,4] = `rf_dif_sd'
mat fgls_rf_cf[4,3] = `rmse'

mat fgls_rf_cf[3,5] = `cf_dif_mean'
mat fgls_rf_cf[3,6] = `cf_dif_sd'


frmttable using "${FIGS}/Comparison_Table", statmat(fgls_rf_cf) substat(1) ctitles("", "FGLS", "Random Forest", "Causal Forest") rtitles("Predicted yield (F=0)" \ "" \ "Predicted yield (F=1)" \ "" \ "Predicted fertilizer response" \ "" \ "RMSE") sdec(2) tex fragment replace 


** Variance decomposition:

use "${TEMP}/TrialSim_rf_processed.dta", clear
loc climvars temp_p1 temp_p2 temp_p3 precip_p1 precip_p2 precip_p3
loc soilvars soilcec soilph claypct siltpct bulkdens acidity soilom elevm
regress y_pred_dif `climvars'
regress y_pred_dif `soilvars'

correlate `climvars'
correlate `soilvars'
correlate `climvars' `soilvars'


** Variance decomposition:
use "${TEMP}/TrialSim_cf_processed.dta", clear
loc climvars temp_p1 temp_p2 temp_p3 precip_p1 precip_p2 precip_p3
loc soilvars soilcec soilph claypct siltpct bulkdens acidity soilom elevm
regress y_pred_dif `climvars'
regress y_pred_dif `soilvars'

correlate `climvars'
correlate `soilvars'
correlate `climvars' `soilvars'


cap log close