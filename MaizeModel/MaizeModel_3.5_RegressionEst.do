
***************************************************************************************
***************************************************************************************
******************************** REGRESSION MODEL **********************************
***************************************************************************************
***************************************************************************************

** This file estimates a regression yield model after selecting variables with vselect (to compare with random forest)

set more off
cap log close

clear mata
clear matrix
mata: mata clear

set seed 123456

loc vselect = 0
loc K = 5

log using "${TEMP}/MaizeModel_RegressionEst.log", replace


foreach vlist in idvar fertvar depvarlog depvarlev xvars_dum xvars_fe xvars_cont_clim xvars_cont_fixed xvars_cont clustvar wvar {
	loc `vlist' "${`vlist'}"
	}

forvalues i = 1/5 {
	loc allvarslev`i' ""
	}

use "RandomForest/TrialData_act_lev", clear
g u = runiform()
sort $clustvar
xtile kg = u, nq(`K')
bysort $clustvar: egen kgroup = mode(kg)
collapse (median) kgroup ps_pro, by($clustvar)
replace kgroup = round(kgroup)
keep $clustvar kgroup ps_pro
loc level lev

tempfile psweights_site_kgroup
save `psweights_site_kgroup', replace

foreach set in act f0 f1 {
	use "RandomForest/TrialData_`set'_lev", clear	
	merge m:1 $clustvar using `psweights_site_kgroup', nogen keep(1 3) keepusing(kgroup ps_pro)
	sort $clustvar
	xtset sitecode
		
	loc k_xc: word count `xvars_cont' 
	loc k_xd: word count `fertvar' `xvars_dum'

	foreach var in `xvars_cont' {
		qui: su `var', detail
		loc `var'_med = r(p50)
		loc allvarslev1 `allvarslev1' `var'
		}

	foreach var in `xvars_dum' {
		loc `var'_med = 0
		loc allvarslev2 `allvarslev2' `var'
		}

	forvalues j = 1/`k_xc' {
		loc var1 `: word `j' of `xvars_cont' '
		forvalues k = `j'/`k_xc' {
			loc var2 `:word `k' of `xvars_cont' '
			g `var1'_x_`var2' = `var1' * `var2'
			loc allvarslev3 `allvarslev3' `var1'_x_`var2'
			}
			forvalues k = 1/`k_xd' {
				loc var2 `:word `k' of `fertvar' `xvars_dum''
				g `var1'_x_`var2' = `var1' * `var2'
				loc allvarslev4 `allvarslev4' `var1'_x_`var2'
				}
		}

		foreach var in `xvars_dum'   {
			g `fertvar'_x_`var' = `fertvar' * `var'
			loc allvarslev5 `allvarslev5' `fertvar'_x_`var' `xvars_fe_dum'
			}

	loc k_allvars: word count `allvarslev1' `allvarslev2' `allvarslev3' `allvarslev4' `allvarslev5'
	su `allvarslev1' `allvarslev2' `allvarslev3' `allvarslev4' `allvarslev5'


	tempfile est_dataset_`set'
	save `est_dataset_`set'', replace
	}


*** Now variable selection

if `vselect'==1 {
	forvalues k = 1/`K' {
	   use `est_dataset_act', clear
	   vselect `depvarlev' `fertvar' `allvarslev1' `allvarslev2' `allvarslev3' `allvarslev4' `allvarslev5' if kgroup != `k'  [aw=ps_pro], forward aic
	   indeplist
	   loc varlist_ols_lev_k`k' `r(X)'
	   loc varlist_fgls_lev_k`k' `varlist_ols_lev_k`k''
	   
		regress `depvarlev' `varlist_ols_lev_k`k'' if kgroup != `k' [iw=ps_pro], vce(cluster $clustvar)
		estimates save "${TEMP}\Vselect_ols_lev_k`k'", replace
		xtreg `depvarlev' `varlist_ols_lev_k`k'' if kgroup != `k' [iw=ps_pro], pa corr(ind) vce(robust)
		estimates save "${TEMP}\Vselect_fgls_lev_k`k'", replace
		}
	}

	
forvalues k = 1/`K' {
    use `est_dataset_act', clear
	estimates use "${TEMP}\Vselect_fgls_lev_k`k'"

	indeplist
	loc varlist_fgls_lev_k`k' `r(X)'

	loc p_fgls_lev_k`k' = e(rank) 
	
	xtreg `depvarlev' `varlist_fgls_lev_k`k'' if kgroup != `k' [iw=ps_pro], pa corr(ind) vce(robust)
	predict dep_k`k'
	
	su dep_k`k', detail
					
	keep if kgroup == `k' /*out of sample*/ 

	regress `depvarlev' dep_k`k' [aw=ps_pro], noconstant vce(cluster site_year_id) 
	loc r2_pseudo_noc_fgls_lev_k`k' = e(r2)
	regress `depvarlev' dep_k`k' [aw=ps_pro], vce(cluster site_year_id) 
	loc r2_pseudo_yesc_fgls_lev_k`k' = e(r2)
	
	g err2_k`k' = (`depvar`level'' - dep_k`k')^2   
	su err2_k`k' 
	loc count = r(N)
	loc sum = r(sum)
	loc rmse_k`k'_out = sqrt(`sum'/`count')
			
	cap drop dep_k`k'		
				
	display "depvar is `level' yield; model is fgls; k is `k';"
	display "k = `k'"
	display "p = `p_fgls_`level'_k`k''"
	display "rmse oos = `rmse_k`k'_out'"
	display "r2pseudo (no constant) = `r2_pseudo_noc_fgls_lev_k`k''"						
	display "r2pseudo (constant) = `r2_pseudo_yesc_fgls_lev_k`k''"						
	use `est_dataset_f1', clear
	estimates use "${TEMP}\Vselect_fgls_lev_k`k'"
	predict dep_k`k'_f1 
	keep id dep_k`k'_f1 
	tempfile pred_k`k'_f1
	save `pred_k`k'_f1', replace
	
	use `est_dataset_f0', clear
	estimates use "${TEMP}\Vselect_fgls_lev_k`k'"
	predict dep_k`k'_f0 
	keep id dep_k`k'_f0 
	
	merge 1:1 id using `pred_k`k'_f1', nogen keep(1 3) keepusing(dep_k`k'_f1)
	
	rename dep_k`k'_f0 y_pred_f0
	rename dep_k`k'_f1 y_pred_f1
	
	g y_pred_dif = y_pred_f1-y_pred_f0
	g kgroup = `k'
	save "${TEMP}\Vselect_fgls_lev_fertresp_k`k'.dta", replace
	}

loc rmse_ave_out = 0.2*(`rmse_k1_out' + `rmse_k2_out' + `rmse_k3_out' + `rmse_k4_out' + `rmse_k5_out')
loc r2pseudo_ave_noc = 0.2*(`r2_pseudo_noc_fgls_lev_k1' + `r2_pseudo_noc_fgls_lev_k2' + `r2_pseudo_noc_fgls_lev_k3' + `r2_pseudo_noc_fgls_lev_k4' + `r2_pseudo_noc_fgls_lev_k5')
loc r2pseudo_ave_yesc = 0.2*(`r2_pseudo_yesc_fgls_lev_k1' + `r2_pseudo_yesc_fgls_lev_k2' + `r2_pseudo_yesc_fgls_lev_k3' + `r2_pseudo_yesc_fgls_lev_k4' + `r2_pseudo_yesc_fgls_lev_k5')

display "Ave RMSE OOS: `rmse_ave_out'"
display "Ave R2 Pseudo OOS With Constant: `r2pseudo_ave_yesc'"
display "Ave R2 Pseudo OOS No Constant: `r2pseudo_ave_noc'"

forvalues k = 1/5 {
	estimates use "${TEMP}\Vselect_fgls_lev_k`k'"
	outreg using "${FIGS}\Vselect_fgls_lev", merge replace tex fragment ctitles("", "cv group `k'") addrows("OOS RMSE", "`rmse_k`k'_out'")
	}


clear
forvalues k = 1/5 {
    append using "${TEMP}\Vselect_fgls_lev_fertresp_k`k'.dta"
	}
	
sort id kgroup
collapse (mean) y_pred_f0 y_pred_f1 y_pred_dif, by(id)
g rmse_oos = `rmse_ave_out'
save "${TEMP}\Vselect_fgls_lev_fertresponse.dta", replace

cap log close
cap log close
