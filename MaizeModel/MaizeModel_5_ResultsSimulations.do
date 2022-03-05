***************************************************************************************
***************************************************************************************
******************************** RESULTS SIMULATIONS **********************************
***************************************************************************************
***************************************************************************************

** This file takes the Africa-wide predictions from the Random Forest model (from R) and analyzes them 

set more off
cap log close

clear mata
clear matrix
mata: mata clear

loc bootcumul_map = 1 /*bootstrap the cumulative distribution of vcr/irr*/
loc bootcumul_sd = 1 /*bootstrap the somers d stochastic dominance*/

loc keepsamesites = 1 /*don't draw a new shortlist of trial sites*/ 
loc randomsitelist "2506 922 15 2520 1020 501 929 1040 2510 1037"
if `keepsamesites' == 1 {
	forvalues j = 1/10 {
		loc site_`j': word `j' of `randomsitelist'
		}
	}

foreach type in rf cf {
	loc vcrThreshold_`type'robust=1.3
	loc irrThreshold_`type'robust=0.3
	}
loc vcrThreshold_naive=2
loc irrThreshold_naive = 1
loc probabilityThreshold = 0.7
loc probabilityThresholdList 10 20 30 40 50 60 70 80 90 


log using "${TEMP}/MaizeModel_ResultsSimulations.log", replace


***************
** Load Data **
***************

** Load RF data

use "${TEMP}/AfricaSites_all_sim_robust.dta", clear
expand 2, generate(fert)
tempfile data_robust
save `data_robust', replace

use "RandomForest/Results_SimData_rf_seg1.dta", clear
forvalues i = 2/10 {
	append using "RandomForest/Results_SimData_rf_seg`i'.dta"
	}

rename africaDataRobust_cellid cellid
rename africaDataRobust_iter iter
rename africaDataRobust_fert fert
rename predictions y_pred_f
rename variance_estimates y_pred_se_f
replace y_pred_se_f = sqrt(y_pred_se_f)

sort cellid iter fert
merge 1:1 cellid iter fert using `data_robust', nogen keep(1 3) keepusing(lat lon $xvars_cont $xvars_dum betamaizeprice betaureaprice betacvmaizeprice_var betaureaprice_var)

reshape wide y_pred_f y_pred_se_f, i(cellid iter) j(fert)
rename y_pred_se_f0 y_pred_f0_se
rename y_pred_se_f1 y_pred_f1_se

g y_pred_dif = y_pred_f1 - y_pred_f0
g y_pred_dif_se = (y_pred_f1_se^2 + y_pred_f0_se^2)^0.5

foreach type in f0 f1 dif {
	clonevar y_pred_`type'_orig = y_pred_`type'
	}

g temp_ave = (2*temp_p1 + temp_p2 + 2*temp_p3)/5
egen precip_tot = rowtotal(precip_p1 precip_p2 precip_p3)

* Screen out dropped sites
merge m:1 cellid using "${TEMP}/AfricaSites_screening.dta", nogen keepusing(toinclude_sim sitetrialprob)
keep if toinclude_sim==1
drop toinclude_sim

sort cellid
drop if iter==.

save "${TEMP}/Simdata_rfrobust.dta", replace


** Load CF data

use "${TEMP}/AfricaSites_all_sim_robust.dta", clear
tempfile data_robust
save `data_robust', replace

use "RandomForest/Results_SimData_cf_seg1.dta", clear
forvalues i = 2/10 {
	append using "RandomForest/Results_SimData_cf_seg`i'.dta"
	}

rename africaDataRobust_cellid cellid
rename africaDataRobust_iter iter
rename predictions y_pred_dif
rename variance_estimates y_pred_dif_se
replace y_pred_dif_se = sqrt(y_pred_dif_se)

sort cellid iter 
merge 1:1 cellid iter using `data_robust', nogen keep(1 3) keepusing(lat lon $xvars_cont $xvars_dum betamaizeprice betaureaprice  betacvmaizeprice_var betaureaprice_var)

foreach type in dif {
	clonevar y_pred_`type'_orig = y_pred_`type'
	}

g temp_ave = (2*temp_p1 + temp_p2 + 2*temp_p3)/5
egen precip_tot = rowtotal(precip_p1 precip_p2 precip_p3)

* Screen out dropped sites
merge m:1 cellid using "${TEMP}/AfricaSites_screening.dta", nogen keepusing(toinclude_sim sitetrialprob)
keep if toinclude_sim==1
drop toinclude_sim

sort cellid
drop if iter==.

foreach var in f0 f0_se f1 f1_se {
	g y_pred_`var' = .
}

save "${TEMP}/Simdata_cfrobust.dta", replace


** Load naive sim data

use "${TEMP}/AfricaSites_all_sim_naive.dta", clear
keep cellid iter lat lon yield_* $xvars_cont_fixed poordrain betamaizeprice betaureaprice betacvmaizeprice_var betaureaprice_var

rename yield_f0 y_pred_f0
g y_pred_f0_se = (yield_f0_sd)^0.5
rename yield_f1 y_pred_f1
g y_pred_f1_se = (yield_f1_sd)^0.5
drop yield_f0_s* yield_f1_s*

g y_pred_dif = y_pred_f1 - y_pred_f0
g y_pred_dif_se = (y_pred_f1_se^2 + y_pred_f0_se^2)^0.5

foreach type in f0 f1 dif {
	clonevar y_pred_`type'_orig = y_pred_`type'
	}

* Screen out dropped sites
merge m:1 cellid using "${TEMP}/AfricaSites_screening.dta", nogen keepusing(toinclude_sim sitetrialprob)
keep if toinclude_sim==1
drop toinclude_sim

merge 1:1 cellid iter using"${TEMP}/Simdata_cfrobust.dta", keep(3) keepusing() nogen 
drop if iter==.

save "${TEMP}/Simdata_naive.dta", replace


** Select a shortlist of trial sites for analysis  /*there is no ps_pro weighting in this*/

use "${TEMP}/TrialSim_cf_processed.dta", clear
merge m:1 sitecode using "${TEMP}/TrialSites_gridloc.dta", nogen keep(3)
g lat = real(latgrid)
g lon = real(longrid)

levelsof sitecode, local(sitecodelist) 
local nsites: word count `sitecodelist'

merge m:1 lat lon using "${TEMP}/AfricaSites_fixed", nogen keep(3) keepusing(betamaizeprice betaureaprice)

if "`keepsamesites'"=="0" {
	generate uidraw = floor(`nsites'*runiform())
	loc i = 1
	loc j = 1
	while `j' < 11 {
		loc k = uidraw in `i'
		loc site: word `k' of `sitecodelist'
		count if sitecode == `site' & betamaizeprice != . & betaureaprice != . 
		if r(N) > 90 {
			loc new = 1
			loc L = `j' - 1
			forvalues l = 1/`L' {
				if `site'==`site_`l'' {
					loc new = 0
					}
				}
			if `new'==1 {
				loc site_`j' = `site'
				loc j = `j' + 1	
				}
			}
		loc i = `i' + 1
		} 

	display "Randomly drawn sites: `site_1' `site_2' `site_3' `site_4' `site_5' `site_6' `site_7' `site_8' `site_9' `site_10'"
	}

display "Sites to use: `site_1' `site_2' `site_3' `site_4' `site_5' `site_6' `site_7' `site_8' `site_9' `site_10'"


* Generate a map of test sites
g sites_shortlist = 0
forvalues j = 1/10 {
	replace sites_shortlist = 1 if sitecode==`site_`j''
	}
keep if sites_shortlist==1
egen loc = mlabvpos(lat lon), matrix(12 12 12 1 1 \\ 12 12 12 12 3 \\ 12 12 12 12 4 \\ 12 12 3 3 3 \\ 12 12 12 3 3)
*replace loc = 11 if sitecode == 961
graph twoway (scatter lat lon if sites_shortlist==1, msize(small) mcolor(black) msymbol(square) mlabel(sitecode) mlabvpos(loc)), /*
	*/legend(off) xsc(r(-20 50)) xlabel(-20(10)50) saving("${TEMP}/TrialSiteShortlist", replace)
graph export  "${FIGS}/TrialSitesShortlist_map.pdf", replace

tempfile siteshortlist_robust
save `siteshortlist_robust', replace

	
** RF: reate Site-specific distributions of fert response and profitability

if `bootcumul_map'==1 {
	foreach sim in rfrobust naive {
		use "${TEMP}/Simdata_`sim'.dta", clear
		forvalues i = 1/$bootreps {
			foreach type in f0 f1 {
				cap drop err
				g err = rnormal(0,y_pred_`type'_se)
				g yp_`type'_`i' = y_pred_`type'_orig + (err*("`sim'"=="rfrobust"))
				}
			g yp_dif_`i' = yp_f1_`i' - yp_f0_`i' 
		
			cap drop err
			g err = rnormal(0,betacvmaizeprice_var^.5)
			g mp_`i' = $maizePrice * betamaizeprice + (err*("`sim'"=="rfrobust")) /*only sample stochastic prices if robust - only sample once per site-year*/				
			}
		
		g up = $ureaPrice * betaureaprice

		order cellid lat lon iter yp_f0_* yp_f1_* yp_dif_* mp_* up	
		keep cellid lat lon iter yp_f0_* yp_f1_* yp_dif_* mp_*  up 
		
		save "${TEMP}/Simdata_reps_`sim'.dta", replace 
		}


	** Reshape to long

	foreach sim in rfrobust naive {
		use "${TEMP}/Simdata_reps_`sim'.dta", clear

		reshape long yp_f0_ yp_f1_ yp_dif_ mp_ , i(cellid iter) j(draw)
 
		g vcr_ = yp_dif_ * mp_ /(up * ${fertilizerTreatmentVolume})
		g irr_ = (yp_dif_ * mp_ - (up * ${fertilizerTreatmentVolume} * (1+${int_rate}))) / (up * ${fertilizerTreatmentVolume} * (1+${int_rate}))
		
		foreach type in f0 f1 dif {
		    rename yp_`type'_ y_pred_sim_`type'
			}
		
		rename mp_ price_sim_maize 
		rename up price_sim_urea
		rename vcr_ vcr_sim
		rename irr_ irr_sim

		* g cr_sim = y_pred_sim_dif * 1/(`vcrThreshold_`sim''*${fertilizerTreatmentVolume}) /*f:m price ratio required for vcr>threshold*/
		g cr_sim = y_pred_sim_dif * 1/((1+`irrThreshold_`sim'')*${fertilizerTreatmentVolume}*(1+${int_rate})) /*f:m price ratio required for irr>threshold*/
		
		
		if "`sim'"=="naive" {
		    * g crnrt_sim = y_pred_sim_dif * 1/(`vcrThreshold_robust'*${fertilizerTreatmentVolume}) /*f:m price ratio required for vcr>threshold*/
			g crnrt_sim = y_pred_sim_dif * 1/((1+`irrThreshold_rfrobust')*${fertilizerTreatmentVolume}*(1+${int_rate})) /*f:m price ratio required for irr>threshold*/
			}
			else if "`sim'"=="rfrobust" {
			    g crnrt_sim = .
				}
		save "${TEMP}/Simdata_mean_`sim'.dta", replace
		}
	}
	
	
** CF -- Create Site-specific distributions of fert response and profitability

if `bootcumul_map'==1 {
	foreach sim in cfrobust {
		use "${TEMP}/Simdata_`sim'.dta", clear
		g err = rnormal(0,y_pred_dif_se)
		g yp_dif_1 = y_pred_dif_orig + (err*("`sim'"=="cfrobust"))
		
		g errpr = rnormal(0,betacvmaizeprice_var^.5)
		g mp_1 = $maizePrice * betamaizeprice + (errpr*("`sim'"=="cfrobust")) /*only sample stochastic prices if robust - only sample once per site-year*/				
		
		g up = $ureaPrice * betaureaprice

		order cellid lat lon iter yp_dif_* mp_* up err 
		keep cellid lat lon iter yp_dif_* mp_*  up err 
		
		save "${TEMP}/Simdata_reps_`sim'.dta", replace 
		}



	** Reshape to long

	foreach sim in cfrobust naive {
		use "${TEMP}/Simdata_reps_`sim'.dta", clear

		reshape long yp_dif_ mp_ , i(cellid iter) j(draw)
 
		g vcr_ = yp_dif_ * mp_ /(up * ${fertilizerTreatmentVolume})  
		g irr_ = (yp_dif_ * mp_ - (up * ${fertilizerTreatmentVolume} * (1+${int_rate}))) / (up * ${fertilizerTreatmentVolume} * (1+${int_rate}))
	
		foreach type in dif {
		    rename yp_`type'_ y_pred_sim_`type'
			}
		
		rename mp_ price_sim_maize 
		rename up price_sim_urea
		rename vcr_ vcr_sim
		rename irr_ irr_sim 

		g cr_sim = y_pred_sim_dif * 1/((1+`irrThreshold_`sim'')*${fertilizerTreatmentVolume}*(1+${int_rate})) /*f:m price ratio required for irr>threshold*/
		
		if "`sim'"=="naive" {
			g crnrt_sim = y_pred_sim_dif * 1/((1+`irrThreshold_cfrobust')*${fertilizerTreatmentVolume}*(1+${int_rate}))
			}
			else if "`sim'"=="cfrobust" {
			    g crnrt_sim = .
				}
		save "${TEMP}/Simdata_mean_`sim'.dta", replace
		}
	}

	
**************************************************************
** Describe Fert Response & Profitability, Robust and Naive **
**************************************************************

** Break even price ratio ...
* Solve for the Fertilizer:Maize price ratio that implies profitability > vcrT
* F:M > y_pred_dif / (vcrT * qfert) -->  profitability 

if `bootcumul_map'==1 {
	foreach sim in cfrobust naive {
		loc collapselist_`sim' ""
		use "${TEMP}/Simdata_mean_`sim'.dta", clear

		foreach var in y_pred_sim_dif vcr_sim irr_sim cr_sim crnrt_sim {
			cumul `var', g(cumul_`var') by(cellid) /*this is OK because all are increasing in fert response*/
			}

		foreach pt in `probabilityThresholdList' {
			loc pti = 100 - `pt'
			foreach var in y_pred_sim_dif vcr_sim irr_sim cr_sim crnrt_sim {
				egen `var'_`pti' =  pctile(`var'), p(`pti') by(cellid) 
				loc collapselist_`sim' `collapselist_`sim'' `var'_`pti'
				}
			}
		
		loc pt = 100-`probabilityThreshold'*100
		foreach var in y_pred_sim_dif vcr_sim irr_sim cr_sim crnrt_sim {
			g `var'_probT = `var'_`pt'
			loc collapselist_`sim' `collapselist_`sim'' `var'_probT 
			}
		
		foreach var in vcr irr {
			g `var'_sim_gtP = .
			replace `var'_sim_gtP = 1 if (`var'_sim >= ``var'Threshold_`sim'') & `var'_sim != .
			replace `var'_sim_gtP = 0 if (`var'_sim < ``var'Threshold_`sim'') & `var'_sim != .
			g `var'nrt_sim_gtP = .
			if "`sim'"=="naive" {
				replace `var'nrt_sim_gtP = 1 if (`var'_sim >= ``var'Threshold_cfrobust') & `var'_sim != .
				replace `var'nrt_sim_gtP = 0 if (`var'_sim < ``var'Threshold_cfrobust') & `var'_sim != . 
				}
			}
		
		g price_sim_ftm = price_sim_urea / price_sim_maize 
		
		collapse (mean) `collapselist_`sim'' vcr_sim_gtP vcrnrt_sim_gtP irr_sim_gtP irrnrt_sim_gtP price_sim_maize price_sim_urea price_sim_ftm (sd) vcr_sim_sd=vcr_sim irr_sim_sd=irr_sim y_pred_sim_dif_sd=y_pred_sim_dif, by(cellid lat lon)
		g sim="`sim'"

		save "${TEMP}/Simdata_cellid_`sim'.dta", replace
		}
	}


*** PROCEED WITH SELECTED RESULTS 

use "${TEMP}/Simdata_cellid_cfrobust.dta", clear  /*USE CF RESULTS NOT RF RESULTS*/
append using "${TEMP}/Simdata_cellid_naive.dta"

save "${TEMP}/Simdata_cellid.dta", replace

foreach var in y_pred_sim_dif_probT vcr_sim_probT irr_sim_probT cr_sim_probT crnrt_sim_probT vcrnrt_sim_gtP irrnrt_sim_gtP price_sim_ftm y_pred_sim_dif_sd {
	su `var' , detail 
	loc `var'_0 = r(min)
	loc `var'_100 = r(max)
	_pctile `var' if sim=="cfrobust", percentiles(10(10)90)
	forvalues i = 1/9 {
		loc pt = `i' * 10
		loc `var'_`pt' = r(r`i')
		}
	cap drop `var'_*
	forvalues pt = 0 (10) 90 {
		loc ptp = `pt' + 10
		g `var'_`pt' = `var' > ``var'_`pt'' & `var' <= ``var'_`ptp''
		loc intens_`pt' = (`ptp'/100)*3
		loc lab_`var'_`pt' = 0.5 * (``var'_`pt'' + ``var'_`ptp'')
		loc lab_`var'_`pt': di %6.2f `lab_`var'_`pt''
		display "var: `var'; `pt' cutoff: ``var'_`pt''; `ptp' cutoff: ``var'_`ptp'' midpoint label: `lab_`var'_`pt''"
		}
	}


** Heat map of profitability at threshold robust vs naive

loc var irr 
keep cellid lat lon sim `var'_sim_probT 
rename `var'_sim_probT `var'_sim_probT_
reshape wide `var'_sim_probT, i(cellid) j(sim) string

g profitable_robust = .
g profitable_naive = .
g profitablenrt_naive = .

replace profitable_robust = 1 if `var'_sim_probT_cfrobust >= ``var'Threshold_cfrobust' & `var'_sim_probT_cfrobust != .
replace profitable_robust = 0 if `var'_sim_probT_cfrobust < ``var'Threshold_cfrobust' & `var'_sim_probT_cfrobust != .
replace profitable_naive = 1 if `var'_sim_probT_naive >= ``var'Threshold_naive' & `var'_sim_probT_naive != .
replace profitable_naive = 0 if `var'_sim_probT_naive < ``var'Threshold_naive' & `var'_sim_probT_naive != .
replace profitablenrt_naive = 1 if `var'_sim_probT_naive >= ``var'Threshold_cfrobust' & `var'_sim_probT_naive != .
replace profitablenrt_naive = 0 if `var'_sim_probT_naive < ``var'Threshold_cfrobust' & `var'_sim_probT_naive != .

foreach type in nono yesyes type1 type2 {
	foreach class in class classnrt {
		g `class'_`type'_probT = .
		}
	}
	
replace class_nono_probT = 1 if (profitable_robust==0 & profitable_naive==0) & (`var'_sim_probT_cfrobust != . & `var'_sim_probT_naive != .)
replace class_nono_probT = 0 if !(profitable_robust==0 & profitable_naive==0) & (`var'_sim_probT_cfrobust != . & `var'_sim_probT_naive != .)
replace class_yesyes_probT = 1 if (profitable_robust==1 & profitable_naive==1) & (`var'_sim_probT_cfrobust != . & `var'_sim_probT_naive != .)
replace class_yesyes_probT = 0 if !(profitable_robust==1 & profitable_naive==1) & (`var'_sim_probT_cfrobust != . & `var'_sim_probT_naive != .)
replace class_type1_probT = 1 if (profitable_robust==1 & profitable_naive==0) & (`var'_sim_probT_cfrobust != . & `var'_sim_probT_naive != .)
replace class_type1_probT = 0 if !(profitable_robust==1 & profitable_naive==0) & (`var'_sim_probT_cfrobust != . & `var'_sim_probT_naive != .)
replace class_type2_probT = 1 if (profitable_robust==0 & profitable_naive==1) & (`var'_sim_probT_cfrobust != . & `var'_sim_probT_naive != .)
replace class_type2_probT = 0 if !(profitable_robust==0 & profitable_naive==1) & (`var'_sim_probT_cfrobust != . & `var'_sim_probT_naive != .)

replace classnrt_nono_probT = 1 if (profitable_robust==0 & profitablenrt_naive==0) & (`var'_sim_probT_cfrobust != . & `var'_sim_probT_naive != .)
replace classnrt_nono_probT = 0 if !(profitable_robust==0 & profitablenrt_naive==0) & (`var'_sim_probT_cfrobust != . & `var'_sim_probT_naive != .)
replace classnrt_yesyes_probT = 1 if (profitable_robust==1 & profitablenrt_naive==1) & (`var'_sim_probT_cfrobust != . & `var'_sim_probT_naive != .)
replace classnrt_yesyes_probT = 0 if !(profitable_robust==1 & profitablenrt_naive==1) & (`var'_sim_probT_cfrobust != . & `var'_sim_probT_naive != .)
replace classnrt_type1_probT = 1 if (profitable_robust==1 & profitablenrt_naive==0) & (`var'_sim_probT_cfrobust != . & `var'_sim_probT_naive != .)
replace classnrt_type1_probT = 0 if !(profitable_robust==1 & profitablenrt_naive==0) & (`var'_sim_probT_cfrobust != . & `var'_sim_probT_naive != .)
replace classnrt_type2_probT = 1 if (profitable_robust==0 & profitablenrt_naive==1) & (`var'_sim_probT_cfrobust != . & `var'_sim_probT_naive != .)
replace classnrt_type2_probT = 0 if !(profitable_robust==0 & profitablenrt_naive==1) & (`var'_sim_probT_cfrobust != . & `var'_sim_probT_naive != .)
	

** Pie chart robust vs naive

count if profitable_robust != . & profitable_naive != .
loc tot = r(N)
	
foreach type in nono yesyes type1 type2 {
	foreach thresh in class classnrt {
		count if `thresh'_`type'_probT == 1
		loc share_`thresh'_`type' = 100 * r(N) / `tot'
		loc share_`thresh'_`type' : di %3.1f `share_`thresh'_`type'' 
		display "thresh: `thresh' type: `type' share: `share_`thresh'_`type'' "
		}
	}
		
graph pie class_type1_probT class_type2_probT class_nono_probT class_yesyes_probT , /*
	*/pie(3, color(red)) pie(4, color(green)) pie(1, explode color(blue)) pie(2, explode color(orange)) /*
	*/plabel(3 "Never" "profitable" "(`share_class_nono'%)", gap(23)) plabel(4 "Always" "profitable" "(`share_class_yesyes'%)", gap(23)) plabel(1 /*"Type 1," */"Profitable" "Robust Only" "(`share_class_type1'%)", gap(23)) plabel(2 /*"Type 2," */"Profitable" "Naive Only" "(`share_class_type2'%)", gap(23)) legend(off)
graph export "${FIGS}/ProfitabilityAnalysis_probT.pdf", replace

graph pie classnrt_type1_probT classnrt_type2_probT classnrt_nono_probT classnrt_yesyes_probT , /*
	*/pie(3, color(red)) pie(4, color(green)) pie(1, explode color(blue)) pie(2, explode color(orange)) /*
	*/plabel(3 "Never" "profitable" "(`share_classnrt_nono'%)", gap(23)) plabel(4 "Always" "profitable" "(`share_classnrt_yesyes'%)", gap(23)) plabel(1 "Profitable" "Robust Only" "(`share_classnrt_type1'%)", gap(23)) plabel(2 "Profitable" "Naive Only" "(`share_classnrt_type2'%)", gap(23)) legend(off)
graph export "${FIGS}/ProfitabilityAnalysis_probT.pdf", replace

keep cellid lat lon profitable_naive profitable_robust class_* profitablenrt_naive classnrt_*

save "${TEMP}/Profitability_confusion.dta", replace


**************************
** Stochastic dominance **
**************************

use "${TEMP}/TrialSim_cf_processed.dta", clear 

merge m:1 sitecode using "${TEMP}/psweights_site.dta", nogen keep(1 3) 
merge m:1 sitecode using "${TEMP}/TrialSites_gridloc.dta", nogen keep(1 3)
g lat = real(latgrid)
g lon = real(longrid)
merge m:1 lat lon using "${TEMP}/AfricaSites_fixed.dta", nogen keep(1 3) keepusing(betamaizeprice betaureaprice  betacvmaizeprice_var betaureaprice_var)

if `bootcumul_sd'==1 {
   expand $bootreps
	g count = 1
	bysort sitecode iter: g draw = sum(count) 
	drop count
	
	g err = rnormal(0,betacvmaizeprice_var^.5)
	g mp_ = $maizePrice * betamaizeprice + err /*only sample stochastic prices for maize - each site gets same maize price for f0 and f1 */

	g up = $ureaPrice * betaureaprice  /*do not sample stochastic urea price -- assume fertilizer price known at time of planting*/

	cap drop err
	g err = rnormal(0,y_pred_dif_se)
	g yp_dif = y_pred_dif + err
		
	g pc_ = (yp_dif * mp_ - up * ${fertilizerTreatmentVolume} * (1+${int_rate}))/ (up * ${fertilizerTreatmentVolume} * (1+${int_rate}))
	
	levelsof sitecode, local(sitecodelist) 
	local nsites: word count `sitecodelist'
	
	foreach p in 1 5 10 25 50 75 90 95 99 {
		g prof_condp`p'_site = .
		}

	g sitename = ""

	forvalues j = 1/`nsites' {
		loc site: word `j' of `sitecodelist' 
		replace sitename = "`site'" in `j'
		count if betamaizeprice!=. & betaureaprice!=. & sitecode==`site'
		if r(N)> 90 {
			su pc_ if sitecode==`site', detail
			foreach p in 1 5 10 25 50 75 90 95 99 {
				replace prof_condp`p'_site = r(p`p') in `j'
				}
			}
		}
	
	loc J1 = `nsites'+1
	replace sitename = "All" in `J1'
	
	su pc_, detail
	foreach p in 1 5 10 25 50 75 90 95 99 {
		replace prof_condp`p'_site = r(p`p') in `J1'
		}
		
	save "${TEMP}/Stochasticity_reps_cf.dta", replace 
	}



use "${TEMP}/Stochasticity_reps_cf.dta", clear 
keep if sitename != ""
keep if prof_condp50_site != .

foreach p in 1 5 10 25 50 75 90 95 99 {
	replace prof_condp`p'_site = prof_condp`p'_site - `irrThreshold_cfrobust'
	}

keep sitename prof_condp*_site 
list sitename prof_condp*_site 

su prof_condp*_site if site != "All"

save "${TEMP}/StochasticDominance_cf.dta", replace


** Graph Shortlist of 10 CDFs
	
use "${TEMP}/Stochasticity_reps_cf.dta", clear
g sites_shortlist = 0

forvalues i = 1/10 {
	replace sites_shortlist=1 if sitecode==`site_`i''
	}

keep if sites_shortlist==1
drop pc_* prof_cond*
cap drop err

order sitecode iter y_pred_dif-betaureaprice_var yp_ mp_ up

g irr_sim = (yp_dif * mp_ - (up * ${fertilizerTreatmentVolume} * (1+${int_rate})))/(up * ${fertilizerTreatmentVolume} * (1+${int_rate}))
cumul irr_sim, generate(cumul_irr_sim) by(sitecode)

loc grcommand ""
loc legcommand ""
forvalues i = 1/10 {
	loc site = `site_`i''
	display "Iteration: `i'; Sitecode is: `site'; "
	su irr if sitecode==`site' & cumul_irr_sim <= 0.99
	loc xpos = r(max)
	count if sitecode == `site' & irr_sim != .
	loc obs = r(N)
	if `obs' > 0 {
		loc grcommand `grcommand' (line cumul_irr_sim irr_sim if sitecode==`site' & cumul_irr_sim >= 0.01 & cumul_irr_sim <= 0.99, text(.99 `xpos' "`site'"/*, place(`loc_`site'')*/))
		loc legcommand `legcommand' label(`i' "Site `site'")
		}
	}
	
loc xline = `irrThreshold_cfrobust'
loc xlinetitle = `xline' + .1
loc yline = 1-`probabilityThreshold'

sort sitecode irr_sim
graph twoway `grcommand',  /*
	*/ xtitle("Internal Rate of Return") ytitle("Cumulative Distribution") /*
	*/ legend(off/*`legcommand'*/) xline(`xline', lpattern(dash) lcolor(black)) yline(`yline', lpattern(dash) lcolor(black)) yscale(r(0 1.05)) ylabel(0 (0.2) 1) xscale(r(-2 10)) xlabel(-2 (2) 10) text(1.05 `xlinetitle' "Profitability threshold", placement(e) just(left)) text(`yline' 6 "Low-Profitability" "Tolerance (1-P)", placement(e)) saving("${TEMP}/TrialSiteShortlist_cumul_cfrobust", replace)
graph export "${FIGS}/TrialSites_Sim_Shortlist_cfrobust.pdf", replace

graph combine "${TEMP}/TrialSiteShortlist.gph" "${TEMP}/TrialSiteShortlist_cumul_cfrobust.gph", rows(1)
graph export "${FIGS}/TrialSiteShortlist_cumul.pdf", replace

keep sitecode iter irr_sim cumul_irr_sim
merge m:1 sitecode using "${TEMP}/TrialSites_gridloc.dta", nogen keep(3) keepusing(latgrid longrid)
g lat = real(latgrid)
g lon = real(longrid)

save "${RESULTS}/TrialSites_Sim_Shortlist", replace
outsheet using "${RESULTS}/TrialSites_Sim_Shortlist.csv",  delimiter(",") replace
	

****************************************
** Describe fert response by country  **
****************************************

insheet using "${DATA}/AfricaDatabase.csv", clear 
keep cellid country
sort cellid

replace country = "South Africa" if country=="Afrique du Sud"
replace country = "Cameroon" if country=="Cameroun"
replace country = "Cote d'Ivoire" if country=="C“te d'Ivoire"
replace country = "Ethiopia" if country=="Ethiopie"
replace country = "Equitorial Guinea" if country=="Guinee equatoriale"
replace country = "Namibia" if country=="Namibie"
replace country = "Uganda" if country=="Ouganda"
replace country = "Central African Republic" if country=="Republique centrafricaine"
replace country = "Democratic Republic of the Congo" if country=="République Democratique du Congo"
replace country = "Republic of the Congo" if country=="République du Congo"
replace country = "Somalia" if country=="Somalie"
replace country = "Sudan" if country=="Soudan"
replace country = "South Sudan" if country=="Sudan (Kenya administrated)"
replace country = "Tanzania" if country=="Tanzanie"
replace country = "Tanzania" if country=="Tanzania, United Rep of"
replace country = "Chad" if country=="Tchad"
replace country = "Zambia" if country=="Zambie"

tempfile countries
save `countries', replace

loc var irr 

use "${TEMP}/Simdata_cellid.dta", clear
keep if sim=="cfrobust"
merge 1:1 cellid using `countries', nogen keep(1 3)

g profitable = `var'_sim_gtP > `probabilityThreshold'
g count_ypred=1
g count_`var'=(`var'_sim_probT != .)
collapse (sum) n_ypred=count_ypred n_`var'=count_`var' (mean) y_pred_sim_dif_probT cr_sim_probT `var'_sim_probT `var'_sim_gtP profitable, by(country)

drop if n_ypred < 25
g sh_`var' = n_`var'/n_ypred
replace sh_`var' = . if sh_`irr' < 0.1
replace `var'_sim_probT = . if sh_`var'==.
replace profitable = . if sh_`var'==.
mkmat n_ypred y_pred_sim_dif_probT cr_sim_probT sh_`var' `var'_sim_probT profitable, matrix(data) rownames(country)

frmttable using  "${FIGS}/SimData_Country.tex", ctitles("", "Number", "Mean", "Mean", "Share cells", "Mean", "Share of price" \ "", "cells", "yield response", "F:M ratio", "with", "IRR", "data cells" \ "" ,"modeled", "(t/ha)", "required for", "price data", "", "robustly"\ "" ,"", "", "profitability", "", "", "profitable") sdec(0,2,2,2,2,2) vlines(1001000) statmat(data) tex replace


***************************
** Sensitivity analysis  **
***************************

foreach xvar in $xvars_cont {
    foreach hilo in hi lo {
		use "${TEMP}/Simdata_reps_cfrobust.dta", clear
		keep cellid iter err mp_1 up
		rename mp_1 mp_
		
		rename cellid africaDataRobust_cellid
		rename iter africaDataRobust_iter

		forvalues i = 1/10 {
			merge 1:1 africaDataRobust_cellid africaDataRobust_iter using "RandomForest/Results_SimData_SA_seg`i'_`xvar'_`hilo'.dta", update nogen 
			}

		rename africaDataRobust_cellid cellid
		rename africaDataRobust_iter iter
		rename predictions y_pred_dif
		rename variance_estimates y_pred_dif_se
		replace y_pred_dif_se = sqrt(y_pred_dif_se)


		* Screen out dropped sites
		sort cellid iter 
		merge m:1 cellid using "${TEMP}/AfricaSites_screening.dta", nogen keepusing(toinclude_sim sitetrialprob)
		keep if toinclude_sim==1
		drop toinclude_sim

		sort cellid
		drop if iter==.
		
		g yp_dif = y_pred_dif + err		
		g irr = (yp_dif * mp_ - (up * ${fertilizerTreatmentVolume} * (1+${int_rate})))/(up * ${fertilizerTreatmentVolume} * (1+${int_rate}))
		
		keep cellid iter yp_dif irr
		
		rename yp_dif yp_dif_`xvar'_`hilo'
		rename irr irr_`xvar'_`hilo'
		
		tempfile SA_cfrobust_`xvar'_`hilo'
		save `SA_cfrobust_`xvar'_`hilo'', replace
		}
	}
	
use "${TEMP}/Simdata_mean_cfrobust.dta", clear
keep cellid iter y_pred_sim_dif irr_sim 
rename y_pred_sim_dif yp_dif_baseline
rename irr_sim irr_baseline

foreach xvar in $xvars_cont {
	foreach hilo in hi lo {
	    merge 1:1 cellid iter using `SA_cfrobust_`xvar'_`hilo'', nogen 
		}
	}
	
* now add price sensitivity
merge 1:1 cellid iter using "${TEMP}/Simdata_reps_cfrobust.dta", nogen keep(1 3) keepusing(mp_1 up)
rename mp_1 mp_
rename up up_

g irr_mp_hi = (yp_dif_baseline * (mp_*1.05) - (up_ * ${fertilizerTreatmentVolume} * (1+${int_rate})))/(up_ * ${fertilizerTreatmentVolume} * (1+${int_rate}))
g irr_mp_lo =  (yp_dif_baseline * (mp_*.95) - (up_ * ${fertilizerTreatmentVolume} * (1+${int_rate})))/(up_ * ${fertilizerTreatmentVolume} * (1+${int_rate}))
g irr_up_hi =  (yp_dif_baseline * mp_ - ((up_*1.05) * ${fertilizerTreatmentVolume} * (1+${int_rate})))/((up_*1.05) * ${fertilizerTreatmentVolume} * (1+${int_rate}))
g irr_up_lo =  (yp_dif_baseline * mp_ - ((up_*.95) * ${fertilizerTreatmentVolume} * (1+${int_rate})))/((up_*.95) * ${fertilizerTreatmentVolume} * (1+${int_rate}))
g irr_int_hi = (yp_dif_baseline * mp_ - (up * ${fertilizerTreatmentVolume} * (1+(${int_rate}*1.05))))/(up * ${fertilizerTreatmentVolume} * (1+(${int_rate}*1.05)))
g irr_int_lo = (yp_dif_baseline * mp_ - (up * ${fertilizerTreatmentVolume} * (1+(${int_rate}*.95))))/(up * ${fertilizerTreatmentVolume} * (1+(${int_rate}*.95)))

cap drop mp_ up_

save "${TEMP}/Sensitivity_cfrobust.dta", replace

outsheet using "${RESULTS}/Sensitivity_cfrobust.csv",  delimiter(",") replace

	

******************
** Export Data  **
*******************

use "${TEMP}/Profitability_confusion.dta", clear

tostring *, replace force

sort cellid
tempfile profitability_confusion
save `profitability_confusion', replace

use "${TEMP}/Simdata_cellid.dta", clear
sort cellid sim

tostring *, replace force

foreach var in y_pred_sim_dif_probT `var'_sim_gtP `var'nrt_sim_gtP `var'_sim_probT cr_sim_probT crnrt_sim_probT { 
	replace `var' = "nan" if `var'=="."
	display "`var'"
	rename `var' `var'_
	}

loc var irr

keep cellid lat lon sim y_pred_sim_dif_probT_ `var'_sim_gtP_ `var'nrt_sim_gtP_ `var'_sim_probT_ cr_sim_probT_ crnrt_sim_probT_

reshape wide y_pred_sim_dif_probT_ `var'_sim_gtP_ `var'nrt_sim_gtP_ `var'_sim_probT_ cr_sim_probT_ crnrt_sim_probT_, i(cellid) j(sim) string

merge 1:1 cellid using `profitability_confusion', keepusing(class_* classnrt_*) keep(3) nogen
foreach type in nono yesyes type1 type2 {
	foreach class in class classnrt {
		replace `class'_`type'_probT = "nan" if `class'_`type'_probT=="."
		}
	}

outsheet using "${RESULTS}/Simdata_cellid.csv",  delimiter(",") replace


** Summary stats of confusion matrix

use "${TEMP}/Profitability_confusion.dta", clear
merge 1:1 cellid using "${TEMP}/AfricaSites_fixed.dta", nogen keep(1 3)
bysort country: su class_yesyes_probT class_nono_probT class_type1_probT class_type2_probT /*type 1 is robust only and type 2 is naive only*/ 
g profsame =  class_yesyes_probT==classnrt_yesyes_probT & class_nono_probT==classnrt_nono_probT & class_type1_probT== classnrt_type1_probT & class_type2_probT==classnrt_type2_probT
collapse (mean) class_yesyes_probT class_nono_probT class_type1_probT class_type2_probT profsame, by(country)
list 


use "${TEMP}/Simdata_cellid_cfrobust.dta", clear
merge 1:1 cellid using "${TEMP}/AfricaSites_fixed.dta", nogen keep(1 3)

g count = cr_sim_probT != .
g cr_lt2pt5 = cr_sim_probT < 2.5
g cr_bt2pt5and4 = cr_sim_probT >= 2.5 & cr_sim_probT < 4
g cr_bt4and6 = cr_sim_probT >= 4 & cr_sim_probT < 6
g cr_bt6and9 = cr_sim_probT >= 6 & cr_sim_probT < 9
g cr_gt9 = cr_sim_probT >= 9 

collapse (sum) count cr_lt* cr_bt* cr_gt*, by(country)

foreach var in cr_lt2pt5 cr_bt2pt5and4 cr_bt4and6 cr_bt6and9 cr_gt9 {
    g sh_`var' = `var'/count
	}
list country sh_* if count >= 10


use "${TEMP}/Simdata_cellid_cfrobust.dta", clear
merge 1:1 cellid using "${TEMP}/AfricaSites_fixed.dta", nogen keep(1 3)

g count = irr_sim_probT != .
collapse (sum) count (mean) irr_sim_probT, by(country)

list if count>10


cap log close
