***************************************************************************************
***************************************************************************************
********************************** PREP DATA FOR R ************************************
***************************************************************************************
***************************************************************************************

** This file gets both the trial data and the Africa-wide data ready for estimation and simulation in R, respectively 

set more off
cap log close

clear mata
clear matrix
mata: mata clear

foreach vlist in idvar fertvar depvarlog depvarlev xvars_dum xvars_fe xvars_cont_clim xvars_cont_fixed xvars_cont clustvar wvar {
	loc `vlist' "${`vlist'}"
	}
log using "${TEMP}/MaizeModel_PrepDataForR.log", replace


************************
** Basic dataset prep **
************************

use "${TEMP}/TrialDataset", clear
drop lat lon
merge m:1 sitecode using "${TEMP}/TrialSites_gridloc.dta", nogen keep(1 3)
g lat = real(latgrid)
g lon = real(longrid)
merge m:1 lat lon using "${TEMP}/AfricaSites_fixed.dta", keep(1 3) /*we lose 225 observations because their locations don't match to grid locations, 146 obs (5 site-yr clusters) that would have been included*/
g toinclude_trial = $toinclude_trial 

drop _m 
keep if toinclude_trial

loc xvars_cont_orig ""
loc xvars_cont_lev ""

foreach xvar in `xvars_cont' {
	qui su `xvar' 
	loc `xvar'_orig_mean = r(mean)
	loc `xvar'_orig_sd = r(sd)
	clonevar `xvar'_orig = `xvar'
	replace `xvar' = (`xvar' - ``xvar'_orig_mean') / ``xvar'_orig_sd' /*Standardize the continuous x values*/
	loc xvars_cont_orig `xvars_cont_orig' `xvar'_orig
	loc xvars_cont_lev `xvars_cont_lev' `xvar'
	}
	*

rename yield levyield 

sort `idvar'  
g id = _n

summarize `depvarlog' `depvarlev' `fertvar' `xvars_cont' `xvars_cont_orig'  `xvars_dum' 
bysort `fertvar': su  `depvarlog' `depvarlev'  `xvars_cont' `xvars_cont_orig'  `xvars_dum' 

xtset sitecode

tempfile TrialDataset
save `TrialDataset', replace

keep id `xvars_cont_orig'
save "${TEMP}/TrialData_orig", replace


******************************
** Propensity Score Fitting **
******************************

** Estimate ps weights using historical mean climate in the site and soil chararacteristics

use `TrialDataset', clear

keep if toinclude_trial==1

sort sitecode

collapse (mean) sitehaslown `xvars_cont_fixed' (median) poordrain, by(sitecode lat lon)

tempfile sitedata_fixed
save `sitedata_fixed', replace

use "${TEMP}/AfricaSites_bt_tempprecip_all.dta", clear 
collapse (mean) `xvars_cont_clim', by(lat lon)  

tempfile sitedata_meanclim
save `sitedata_meanclim', replace


use `sitedata_fixed', clear
merge m:1 lat lon using `sitedata_meanclim', nogen keep(3)
 
probit sitehaslown `xvars_cont_fixed' `xvars_cont_clim' i.(poordrain)  

predict ps_pro 
su ps_pro if sitehaslown==1
loc min_ps_pro = r(min)

count if ps_pro < `min_ps_pro'
count if ps_pro < 0.1 /* see https://doi.org/10.1093/biomet/asn055 */
count if ps_pro > 0.9

sort sitecode
keep sitecode ps_*

tempfile psweights_site
save `psweights_site', replace

save "${TEMP}/psweights_site.dta", replace


use `TrialDataset', clear
keep if toinclude_trial==1

probit fert `xvars_cont_fixed' `xvars_cont_clim' i.(`xvars_dum')

predict ps_obs

su ps_obs if fert==0
loc min_ps_obs = r(min)

count if ps_obs < `min_ps_obs'
count if ps_obs < 0.1 /* see https://doi.org/10.1093/biomet/asn055 */
count if ps_obs > 0.9

tempfile psweights_obs
save `psweights_obs', replace

save "${TEMP}/psweights_obs.dta", replace


***********************************************
** Interactions and varlists - Trial Dataset **
***********************************************

use `TrialDataset', clear
rename fert fertact
expand 3
g count = 1
bysort id : g fertstatus = sum(count) 
g fert = fertact if fertstatus==1
g type = "act" if fertstatus==1
replace fert = 0 if fertstatus==2
replace type = "f0" if fertstatus==2
replace fert = 1 if fertstatus==3
replace type = "f1" if fertstatus==3

drop fertstatus count fertact

loc xvars_fe_dum ""
foreach fevar in `xvars_fe' {
	levelsof `fevar', local(fevals)
	foreach feval in `fevals' {
		g `fevar'_`feval' = `fevar'==`feval'
		loc xvars_fe_dum `xvars_fe_dum' `fevar'_`feval'
		}
	}

loc k_xc: word count `xvars_cont_lev' 
loc k_xd: word count `fertvar' `xvars_dum'

foreach var in `xvars_cont_lev' {
	qui: su `var', detail
	loc `var'_med = r(p50)
	loc allvarslev1 `allvarslev1' `var'
	}

foreach var in `xvars_dum' {
	loc `var'_med = 0
	loc allvarslev2 `allvarslev2' `var'
	}

loc k_allvars: word count `allvarslev1' `allvarslev2' 
display "`allvarslev1' `allvarslev2' 
su `allvarslev1' `allvarslev2' 

loc firstxvar: word 1 of `allvarslev1'
display "`firstxvar'"
loc ind: word count `allvarslev2'
loc lastvar: word `ind' of `allvarslev2'
display "`lastxvar'"

sort id 

merge m:1 sitecode using "${TEMP}/psweights_site.dta", nogen keep(1 3) keepusing(ps_pro)
merge m:1 id using "${TEMP}/psweights_obs.dta", nogen keep(1 3) keepusing(ps_obs)

keep id site_year_id sitecode type sitehaslown ps_pro ps_obs `depvarlev' `fertvar' `allvarslev1' `allvarslev2' 
order id site_year_id sitecode type sitehaslown ps_pro ps_obs `depvarlev' `fertvar' `allvarslev1' `allvarslev2' 

foreach var in site_year_id sitecode sitehaslown ps_pro ps_obs `depvarlev' `fertvar' `allvarslev1' `allvarslev2' {
	drop if `var'==.
	}
	
tempfile dataprep_all
save `dataprep_all', replace

foreach type in act f0 f1 {
	use `dataprep_all', clear
	keep if type=="`type'"
	sort id
	saveold "RandomForest/TrialData_`type'_lev", replace version(12)
	}

forvalues i = 1/5 {
    gl allvarslev`i' `allvarslev`i''
	}
	
	
********************************************************************
** Create a dataset that combines trial sites and historical clim **
********************************************************************

use "${TEMP}/AfricaSites_bt_tempprecip_all.dta", clear 
su simyr
loc yrmax = r(max)
loc yrmin = r(min)
loc nyears = 1 + `yrmax' - `yrmin'

use `sitedata_fixed', clear
expand `nyears'
g count = 1
bysort sitecode: g iter = sum(count)
g simyr = .
forvalues i = 1/`nyears' {
    replace simyr = `i' + `yrmin' - 1 if iter==`i'
	}
	
drop count 

merge m:1 lat lon simyr using "${TEMP}/AfricaSites_bt_tempprecip_all.dta", nogen keep(1 3) 

isid sitecode simyr
expand 2, generate(`fertvar')
g var_hy = 1
g wortmann = 0

* Generate x vars in simulated dataset
foreach xvar in `xvars_cont_clim' { /*the soil vars were already standardized*/
	clonevar `xvar'_orig = `xvar'
	replace `xvar' = (`xvar' - ``xvar'_orig_mean') / ``xvar'_orig_sd'
	su `xvar'_orig `xvar', detail
	display "Orig `xvar' mean is: ``xvar'_orig_mean' and sd is: ``xvar'_orig_sd'"
	}

loc k_xc: word count `xvars_cont_lev' 
loc k_xd: word count `fertvar' `xvars_dum' 

keep sitecode iter simyr `fertvar' `allvarslev1' `allvarslev2' 
order sitecode iter simyr `fertvar' `allvarslev1' `allvarslev2' 

saveold "RandomForest/SimData_trialsim_lev", replace version(12)

keep if `fertvar'==1
drop `fertvar'

saveold "RandomForest/SimData_trialsim_lev_cf", replace version(12)


******************************
** Prep Simulation datasets **
******************************

** Robust Simulation dataset prep AND sensitivity analysis dataset prep
/**/
foreach sim in robust {
    display "sim is: `sim'"
	forvalues seg = 1/10 {
		use "${TEMP}/AfricaSites_all_sim_`sim'.dta", clear
		sort cellid iter
		merge m:1 cellid using "${TEMP}/AfricaSites_screening.dta", keep(3) nogen 
		
		loc simlo = (`seg'-1)*100 + 1
		loc simhi = `simlo' + 99
		
		keep if toinclude_sim == 1 & iter >= `simlo' & iter <= `simhi'
		
		expand 2, generate(`fertvar')
		
		keep cellid iter simyr `fertvar' `xvars_dum' `xvars_fe' `xvars_cont'

		* Generate x vars in simulated dataset
		foreach xvar in `xvars_cont' {
			su `xvar', detail
			display "Orig `xvar' mean is: ``xvar'_orig_mean' and sd is: ``xvar'_orig_sd'"
			replace `xvar' = (`xvar' - ``xvar'_orig_mean') / ``xvar'_orig_sd'
			su `xvar', detail
			}

		loc k_xc: word count `xvars_cont_lev' 
		loc k_xd: word count `fertvar' `xvars_dum' 

			
		keep cellid iter simyr `fertvar' `allvarslev1' `allvarslev2' 
		order cellid iter simyr `fertvar' `allvarslev1' `allvarslev2' 
		
		saveold "RandomForest/SimData_robust_seg`seg'", replace version(12)
		
		keep if `fertvar'==1
		drop `fertvar'

		saveold "RandomForest/SimData_robust_cf_seg`seg'", replace version(12)
		
		foreach var in `xvars_cont'  {
		    cap drop v_orig
		    g v_orig = `var'*``var'_orig_sd'+``var'_orig_mean'
			replace `var' = (v_orig*.95 - ``var'_orig_mean')/``var'_orig_sd'
			saveold "RandomForest/SimData_SA_seg`seg'_`var'_lo", replace version(12)
			replace `var' = (v_orig*1.05 - ``var'_orig_mean')/``var'_orig_sd'
			saveold "RandomForest/SimData_SA_seg`seg'_`var'_hi", replace version(12)
			}
		}
	}
	
	
******************************
** Summarize the trial data **
******************************

use  `TrialDataset', clear
merge m:1 sitecode using "${TEMP}/psweights_site.dta", nogen keep(1 3) keepusing(ps_pro)
	
eststo nofert: quietly estpost summarize `depvarlev' `xvars_cont_orig'  if fert==0 [iw=ps_pro]
eststo yesfert: quietly estpost summarize `depvarlev' `xvars_cont_orig'  if fert==1 [iw=ps_pro]
eststo difffert: quietly estpost ttest `depvarlev' `xvars_cont_orig' , by(fert) unequal 

esttab nofert yesfert difffert using "${FIGS}/Summary_ByFert_tstat", ///
	onecell nodepvar tex frag ///
	cells("mean(pattern(1 1 0) fmt(2)) b(star pattern(0 0 1) fmt(2))" "sd(pattern(1 1 0) par fmt(2)) t(pattern(0 0 1) par fmt(2))") ///
	mtitles("No Fertilizer" "Optimal Fertilizer" "T-test") nonumbers ///
	substitute(degrees $^{\circ}\$ ) ///
	label replace 

mata: mata clear

foreach f in 0 1 {
	count if fert==`f'
	loc N_fert`f' = r(N)
	}
	
* normalized diff weighted only
foreach var in levyield `fertvar' `xvars_cont_orig' `xvars_dum'  {
	mat newrow = J(1,6,.)
	su `var' if fert==0 [iw=ps_pro]
	mat newrow[1,1] = r(mean)
	mat newrow[1,2] = r(sd)
	loc `var'_var0 = r(Var)
	loc N_fert0: display %6.0fc r(N)
	su `var' if fert==1 [iw=ps_pro]
	mat newrow[1,3] = r(mean)
	mat newrow[1,4] = r(sd)
	loc `var'_var1 = r(Var)
	loc N_fert1: display %6.0fc r(N)
	mat newrow[1,5] = (newrow[1,3] - newrow[1,1]) / sqrt(``var'_var0' + ``var'_var1')
	loc lbl: variable label `var'
	frmttable, append(Summary_ByFert_ndiff) statmat(newrow) substat(1) sdec(2) nocoltitl rtitles("`lbl'")
	}

outreg using "${FIGS}/Summary_ByFert_ndiff.tex", replay(Summary_ByFert_ndiff) /*
	*/ctitles("","No", "Optimal", "Normalized" \ "","Fertilizer", "Fertilizer", "Difference" ) /*
	*/addrows("N", "`N_fert0'", "`N_fert1'") /*
	*/tex fragment replace

mata: mata clear

* normalized diff unweighted and weighted
foreach var in levyield `fertvar' `xvars_cont_orig' `xvars_dum'  {
	mat newrow = J(1,12,.)
	
	* unweighted part
	su `var' if fert==0
	mat newrow[1,1] = r(mean)
	mat newrow[1,2] = r(sd)
	loc `var'_var0 = r(Var)
	su `var' if fert==1
	mat newrow[1,3] = r(mean)
	mat newrow[1,4] = r(sd)
	loc `var'_var1 = r(Var)
	mat newrow[1,5] = (newrow[1,3] - newrow[1,1]) / sqrt(``var'_var0' + ``var'_var1')
	
	* weighted part
	su `var' if fert==0 [iw=ps_pro]
	mat newrow[1,7] = r(mean)
	mat newrow[1,8] = r(sd)
	loc N_wt_fert0 = r(N)
	loc `var'_wt_var0 = r(Var)
	su `var' if fert==1 [iw=ps_pro]
	mat newrow[1,9] = r(mean)
	mat newrow[1,10] = r(sd)
	loc `var'_wt_var1 = r(Var)
	loc N_wt_fert1 = r(N)
	mat newrow[1,11] = (newrow[1,9] - newrow[1,7]) / sqrt(``var'_wt_var0' + ``var'_wt_var1')

	loc lbl: variable label `var'
	frmttable, append(Summary_ByFert_wtndiff) statmat(newrow) substat(1) sdec(2) nocoltitl rtitles("`lbl'")
	}

	
outreg using "${FIGS}/Summary_ByFert_ndiff_wt.tex", replay(Summary_ByFert_wtndiff) /*
	*/ctitles("","No", "Optimal", "Normalized", "No", "Optimal", "Normalized" \ "","Fertilizer", "Fertilizer", "Difference", "Fertilizer", "Fertilizer", "Difference" \ "", "", "", "", "(weighted)", "(weighted)", "(weighted)") /*
	*/addrows("N", "`N_fert0'", "`N_fert1'", "", "`N_wt_fert0'", "`N_wt_fert1'") /*
	*/tex fragment replace
	
clear
set obs 1
loc tex "${FIGS}/Summary_ByFert_ndiff.tex"
loc tex2  "${FIGS}/Summary_ByFert_ndiff_tab.tex"

generate strL s = fileread("`tex'") if fileexists("`tex'")
assert filereaderror(s)==0
replace s = subinstr(s,"\begin{tabular}{lccc}","\begin{tabular}{L{5.5cm}ccc}",1)
replace s = subinstr(s,"degrees","$^{\circ}$",.)
gen byte fw = filewrite("`tex2'",s,1)

clear
set obs 1
loc tex "${FIGS}/Summary_ByFert_ndiff_wt.tex"
loc tex2  "${FIGS}/Summary_ByFert_ndiff_wt_tab.tex"

generate strL s = fileread("`tex'") if fileexists("`tex'")
assert filereaderror(s)==0
replace s = subinstr(s,"degrees","$^{\circ}$",.)
gen byte fw = filewrite("`tex2'",s,1)

cap log close
	
	
