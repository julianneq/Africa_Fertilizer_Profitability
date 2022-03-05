***************************************************************************************
***************************************************************************************
********************************** CREATE SIM DATA ************************************
***************************************************************************************
***************************************************************************************

** This file compiles the simulation dataset (for all African sites)
** It replicates all of the estimation variables across the African sites  


set more off
cap log close

loc genAfricaBootTempPrecip=0 /*This is slow - set to 0 if not needed*/
loc genAfricaWGenTempPrecip=0 /*This is slow - set to 0 if not needed*/

loc mlength = 365/12

log using "${TEMP}/MaizeModel_CreateSimData.log", replace


********************************************************
** First, set up the fixed part of the Africa dataset **
********************************************************

insheet using "${DATA}/AfricaDatabase.csv", clear /*just for the country list*/
keep cellid country
sort cellid
tempfile countries
save `countries', replace

insheet using "${DATA}/AfricaDatabase.csv", clear
rename elev elevm
replace elevm = 0 if elevm < 0

label define aezs 0 "Subtropic - warm / arid" 1 "Subtropic - warm / semiarid" 2 "Subtropic - warm / subhumid" /*
*/3 "Subtropic - warm / humid" 4 "Subtropic - cool / arid" 5 "Subtropic - cool / semiarid" 6 "Subtropic - cool / subhumid" /*
*/7 "Tropic - warm / arid" 8 "Tropic - warm / semiarid" 9 "Tropic - warm / subhumid" 10 "Tropic - warm / humid" /*
*/11 "Tropic - cool /arid" 12 "Tropic - cool / semiarid" 13 "Tropic - cool / subhumid" 14 "Tropic - cool / humid", replace
rename aezid aez
la values aez aezs

*generate aez dummies
levelsof aez, local(aezvals)
foreach aezval in `aezvals' {
	g aez_`aezval' = aez==`aezval'
	}

rename bldsd2 bulkdens 
replace bulkdens = bulkdens/1000 /*See ISRC report*/

rename cecsd2 soilcec /*values seem well scaled according to ISRIC report*/

rename clyppt claypct
rename sltppt siltpct
rename sndppt sandpct
rename crfvol coarsefrag
foreach soilvar in claypct siltpct sandpct coarsefrag {
	replace `soilvar' = `soilvar' / 100
	}

g poordrain = median_drainfao == 1 | median_drainfao ==2 | median_drainfao==3 
/* Raw codes: 1:V, 2:P, 3:I, 4:M, 5:W, 6:S, 7:E from http://www.isric.org/sites/default/files/ISRIC_Report_%202015_02.pdf */

rename eackcl acidity /*there are only a few zero values but this seems feasible according to data cleaning protocol from the ISRIC report. Scaling consistent also */
rename nto soiln /*range is .17 to 6.2, this seems in line with ISRIC */
rename orcdrc soilom /*range is OK*/
rename phihox soilph /*range is 45 to 84, it should be 4.5 to 8.4 */
replace soilph = soilph / 10


** Process price data

rename betamaizeprize betamaizeprice
foreach var in betamaizeprice betacvmaizeprice betaureaprice betacvureaprice {
	recode `var' (-9999=.)
	}
g betacvmaizeprice_var = betacvmaizeprice * betamaizeprice
g betaureaprice_var = betacvureaprice * betaureaprice

drop betacvmaizeprice betacvureaprice


*Later we sample stochastically from prices

g price_maize = betamaizeprice * $maizePrice
g price_urea = betaureaprice * $ureaPrice


** Generate variety dummies (so this simulation is basically for quick maturing hybrids, or EIHY)

g var_EPOP = 0 
g var_ILPO = 0 
g var_ILHY = 0
g var_hy = 1
g wortmann = 0 /*make them cimmyt comparable*/


** Label variables

la var poordrain "Poor drainage (dummy)"
la var acidity "Soil exchangeable acidity (centimol charge per kg soil)"
la var soilph "Soil pH (pH determined in soil/water mixture)"
la var soilcec "Soil cation exchange capacity (centimol charge per kg soil)"
la var bulkdens "Soil bulk density (kg per cubic decimeter)"
la var coarsefrag "Soil coarse fragment contents (share by volume)"
la var soiln "Total soil nitrogen (g/kg soil)"
la var soilom "Soil organic matter (g/kg soil)"
la var claypct "Soil clay (share by volume)"
la var siltpct "Soil silt (share by volume)"
la var sandpct "Soil sand (share by volume)"
la var aez "Agroecological zone"
la var price_maize "Maize price (USD/mt)"
la var price_urea "Fertilizer price (USD/mt)"
la var elevm "Elevation (M)"

isid cellid

tempfile africasites_fixed
save `africasites_fixed', replace

save "${TEMP}/AfricaSites_fixed.dta", replace


**************************************************************
** Set up monthly temp and precip data for all of the years **
**************************************************************

use "${TEMP}/AfricaSites_fixed.dta", clear

keep cellid lat lon month day
sort cellid

count
loc sitesN = r(N)

format %9.3f lat lon

tempfile africasiteroster
save `africasiteroster', replace

loc grp1start = 1
loc grp1end = 4100
loc grp2start = 4101
loc grp2end = 8200
loc grp3start = 8201
loc grp3end = 12300
loc grp4start = 12301
loc grp4end = `sitesN'

if `genAfricaBootTempPrecip' == 1 { /*NOTE THIS IS VERY SLOW... */
	quietly {
		forvalues grp=1/1 {
			loc firstsite = 1
			forvalues i = `grp`grp'start' / `grp`grp'end'  {
				use `africasiteroster', clear
				sort cellid
				foreach latlon in lat lon {
					loc `latlon' = `latlon' in `i' 
					if ``latlon''<1 & ``latlon''>0 {
						loc `latlon'fmt 0``latlon''
						}
						else if ``latlon''>-1 & ``latlon''<0 {
							loc `latlon'fmt = substr("``latlon''",-3,.)
							loc `latlon'fmt -0.``latlon'fmt'
							}
							else if ``latlon'' > 1 | ``latlon''< -1 {
								loc `latlon'fmt ``latlon''
								}
					}

				loc month = month in `i'
				loc day = day in `i'

				cap import delimited using  "${CLIMDATA}/TempPrecip/lat`latfmt'_long`lonfmt'.txt", delimiters(" ", collapse) varnames(4) stringcols(1) numericcols(2 3 4 5 6 7) clear
				display "Lat is `lat'; lon is `lon'; formatted is `latfmt' `lonfmt'; return code is: `rc'"
				display _rc
				
				if _rc == 0 {
					g mdy = mdy(1,1,1980) + real(substr(date,3,.))
					g mo = month(mdy)
					g yr = year(mdy)
					g day = day(mdy)
					
					format mdy %d
					drop date
					rename mdy date
					rename tavg temp 
					
					recode temp tmin tmax rain wind (-99=.)
					
					su yr
					loc yrmin = r(min)
					loc yrmax = r(max)
					
					forvalues yr = `yrmin'/`yrmax' {
						loc date = mdy(`month',`day',`yr')
						display "Date: `date'"
						g gmonth = 0
						
						forvalues k = 1/5 {
							replace gmonth = `k' if (date >= int(`date'+(`k'-1)*`mlength')) & (date < int(`date'+`k'*`mlength'))
							}
							
						forvalues k = 1/5 {
							g precip_m`k'_`yr' = rain * (gmonth==`k')
							g temp_m`k'_`yr' = temp * (gmonth==`k')
							g gdd_m`k'_`yr' = (min(0.5*(tmin+tmax),${TM}) - ${TB})*(gmonth==`k')
							g hdd_m`k'_`yr' = (max(tmax,${TM}) - ${TM})*(gmonth==`k')
							}
						
						drop gmonth
						}
					
					recode temp_m* gdd_m* (0=.)
					
					collapse (sum) precip_m1_* precip_m2_* precip_m3_* precip_m4_* precip_m5_* gdd_m1_* gdd_m2_* gdd_m3_* gdd_m4_* gdd_m5_* hdd_m1_* hdd_m2_* hdd_m3_* hdd_m4_* hdd_m5_* (mean) temp_m1_* temp_m2_* temp_m3_* temp_m4_* temp_m5_* 
					
					g lat = `lat'
					g lon = `lon'
					
					reshape long precip_m1_ precip_m2_ precip_m3_ precip_m4_ precip_m5_ temp_m1_ temp_m2_ temp_m3_ temp_m4_ temp_m5_ gdd_m1_ gdd_m2_ gdd_m3_ gdd_m4_ gdd_m5_ hdd_m1_ hdd_m2_ hdd_m3_ hdd_m4_ hdd_m5_ , i(lat lon) j(simyr)

					if `firstsite' == 1 {
						tempfile africasites_bt_tempprecip_grp`grp'
						save `africasites_bt_tempprecip_grp`grp'', replace
						}

					else if `firstsite' == 0 {
						append using `africasites_bt_tempprecip_grp`grp''
						tempfile africasites_bt_tempprecip_grp`grp'
						save `africasites_bt_tempprecip_grp`grp'', replace
						}
					
					loc firstsite = 0
					}
				}

			save "${TEMP}/AfricaSites_bt_tempprecip_grp`grp'.dta", replace
			erase `africasites_bt_tempprecip_grp`grp''
			}
		}
	}	
			
use "${TEMP}/AfricaSites_bt_tempprecip_grp1.dta", clear
cap append using "${TEMP}/AfricaSites_bt_tempprecip_grp2.dta"
cap append using "${TEMP}/AfricaSites_bt_tempprecip_grp3.dta"
cap append using "${TEMP}/AfricaSites_bt_tempprecip_grp4.dta"

duplicates tag lat lon simyr, gen(dup)
tab dup
drop dup

duplicates drop lat lon simyr, force
merge m:1 lat lon using `africasiteroster', nogen keep(3)

egen temp_p1 = rowmean(temp_m1 temp_m2)
egen temp_p2 = rowmean(temp_m3)
egen temp_p3 = rowmean(temp_m4 temp_m5)

egen gdd_p1 = rowtotal(gdd_m1 gdd_m2)
egen gdd_p2 = rowtotal(gdd_m3)
egen gdd_p3 = rowtotal(gdd_m4 gdd_m5)

count if temp_p1 == 0
count if temp_p2 == 0
count if temp_p3 == 0

recode temp_p1 temp_p2 temp_p3 gdd_p1 gdd_p2 gdd_p3 (0=.) /*Mean temp can't feasibly be zero degrees celcius*/

egen precip_p1 = rowtotal(precip_m1 precip_m2)
egen precip_p2 = rowtotal(precip_m3)
egen precip_p3 = rowtotal(precip_m4 precip_m5) 
egen precip_tot = rowtotal(precip_p1 precip_p2 precip_p3)

count if precip_tot == 0
replace precip_p1 = . if precip_tot == 0
replace precip_p2 = . if precip_tot == 0
replace precip_p3 = . if precip_tot == 0

foreach var in temp precip gdd {
	forvalues m = 1/5 {
		rename `var'_m`m'_ `var'_m`m'
		}
	}

drop if simyr==2010 /*too many missing values*/

save "${TEMP}/AfricaSites_bt_tempprecip_all.dta", replace
*


**********************************************************************************************
** Set up synthetic growing season temp and precip data for the weather generated clim data **
**********************************************************************************************

if `genAfricaWGenTempPrecip' == 1 {
	quietly {
		forvalues grp=1/1 {
			loc firstsite = 1
			forvalues i = `grp`grp'start' / `grp`grp'end' {
				display "i is: `i'"
				use `africasiteroster', clear
				sort cellid
				
				foreach latlon in lat lon {
					loc `latlon' = `latlon' in `i' 
					if ``latlon''<1 & ``latlon''>0 {
						loc `latlon'fmt 0``latlon''
						}
						else if ``latlon''>-1 & ``latlon''<0 {
							loc `latlon'fmt = substr("``latlon''",-3,.)
							loc `latlon'fmt -0.``latlon'fmt'
							}
							else if ``latlon'' > 1 | ``latlon''< -1 {
								loc `latlon'fmt ``latlon''
								}
					}
				display "lat is: `lat' lon is: `lon'."
				
				cap import delimited using  "${CLIMDATA}/SyntheticTempPrecipCRU/lat`latfmt'_long`lonfmt'_temp.txt", delimiters(" ", collapse) clear
				if _rc != 0 {
					display "file not found"
					}
				else if _rc == 0 {
				    display "file found"
					g simyr = _n
					
					foreach m in 1 2 3 4 5 {
						rename v`m' temp_m`m'
						}

					keep in 1/$simreps 
					
					tempfile site_temp
					save `site_temp', replace
					
					cap noisily import delimited using  "${CLIMDATA}/SyntheticTempPrecipCRU/lat`latfmt'_long`lonfmt'_precip.txt", delimiters(" ", collapse) clear
					g simyr = _n

					foreach m in 1 2 3 4 5 {
						rename v`m' precip_m`m'
						}

					keep in 1/$simreps

					merge 1:1 simyr using `site_temp', nogen keep(1 3)
					
					recode temp_* precip_* (-99=.)
					egen precip_p1 = rowtotal(precip_m1 precip_m2)
					egen precip_p2 = rowtotal(precip_m3)
					egen precip_p3 = rowtotal(precip_m4 precip_m5)
					egen precip_tot = rowtotal(precip_p1 precip_p2 precip_p3)
					
					egen temp_p1 = rowmean(temp_m1 temp_m2)
					egen temp_p2 = rowmean(temp_m3)
					egen temp_p3 = rowmean(temp_m4 temp_m5)
	
					g lat = `latfmt'
					gen lon = `lonfmt'

					if `firstsite'==1 {
						tempfile africasites_wg_tempprecip_grp`grp'
						save `africasites_wg_tempprecip_grp`grp'', replace
						}

					else if `firstsite' == 0 {
						append using `africasites_wg_tempprecip_grp`grp''
						tempfile africasites_wg_tempprecip_grp`grp'
						save `africasites_wg_tempprecip_grp`grp'', replace
						}
						
					loc firstsite=0		
					}
				}
			
			save "${TEMP}/AfricaSites_wg_tempprecip_grp`grp'.dta", replace	
			erase `africasites_wg_tempprecip_grp`grp''
			}
		}
	}


use "${TEMP}/AfricaSites_wg_tempprecip_grp1.dta", clear
cap append using "${TEMP}/AfricaSites_wg_tempprecip_grp2.dta"
cap append using "${TEMP}/AfricaSites_wg_tempprecip_grp3.dta"
cap append using "${TEMP}/AfricaSites_wg_tempprecip_grp4.dta"
duplicates drop lat lon simyr, force

merge m:1 lat lon using `africasiteroster', nogen keep(3)

sort cellid simyr
isid cellid simyr

save "${TEMP}/AfricaSites_wg_tempprecip_all.dta", replace


********************************************************************
** Figure out which sites should be screened out by temp / precip **
********************************************************************

** Screening based on propensity score

use "${TEMP}/AfricaSites_bt_tempprecip_all.dta", clear 
egen temp_ave = rowmean(temp_m1 temp_m2 temp_m3 temp_m4 temp_m5)

collapse (mean) temp_ave precip_tot ${xvars_cont_clim} (sd) temp_sd=temp_ave precip_sd=precip_tot , by(cellid) 

merge 1:1 cellid using "${TEMP}/AfricaSites_fixed.dta", nogen keepusing(lat lon ${xvars_cont_fixed} poordrain aez) keep(1 3)
tempfile africa_all_mean
save `africa_all_mean', replace

use "${TEMP}/TrialSites_gridloc.dta", clear
rename latgrid lat
rename longrid lon
destring lat lon, replace
sort lat lon
merge m:1 lat lon using `africa_all_mean', nogen keep(3) keepusing(cellid)
isid sitecode
tempfile site_cell
save `site_cell', replace

use  "${TEMP}/TrialDataset", clear
merge m:1 sitecode using `site_cell', keepusing(cellid) keep(1 3)
g trialspersite=1
collapse (sum) trialspersite, by(cellid)
count
tab trialspersite
g sitehastrials = trialspersite >= 1 & trialspersite != .
tempfile trials_cell
save `trials_cell', replace

use `africa_all_mean', clear
merge 1:1 cellid using `trials_cell', nogen keep(1 3)
recode trialspersite (.=0)
recode sitehastrials (.=0)

* LPM
regress sitehastrials ${xvars_cont_clim} ${xvars_cont_fixed} temp_sd precip_sd poordrain 
predict sitetrialprob

scatter sitetrialprob sitehastrials
su sitetrialprob, detail
bysort sitehastrials: su sitetrialprob, detail
su sitetrialprob if sitehastrials==1
loc trialprob_min = r(min)
display "Min lp for fertilized site: `trialprob_min'"
su sitetrialprob if sitehastrials, detail /*we lose 25% of sites*/
g toinclude_sim = sitetrialprob >= `trialprob_min' & sitetrialprob != .

graph twoway (scatter lat lon if toinclude_sim==1, msize(vtiny) mcolor(black)) (scatter lat lon if !(toinclude_sim==1), msize(vtiny) mcolor(red))
heatplot toinclude_sim lat lon, cuts(0 (.1) 1) colors(hsv, sequential reverse) bwidth(0.5)
graph export "${TEMP}/AfricaSites_screening.pdf", replace

keep cellid lat lon sitehastrials sitetrialprob toinclude_sim
sort cellid

tempfile sitetrialprob
save `sitetrialprob', replace

sort cellid
save "${TEMP}/AfricaSites_screening.dta", replace


**************************************
** Create Robust Simulation Dataset **
**************************************

clear
set obs $simreps
g iter=_n
g simyr=int(runiform()*$simreps) + 1 /*note - we sample w replacement from synthetic climate dataset*/
tempfile simdraws_robust
save `simdraws_robust', replace

use `africasites_fixed', clear
expand $simreps
g count=1
bysort cellid: g iter=sum(count)
drop count

merge m:1 iter using `simdraws_robust', nogen keepusing(simyr)
merge m:1 cellid simyr using "${TEMP}/AfricaSites_wg_tempprecip_all.dta", nogen keep(1 3) /*note - we sample from the simulated data with replacement*/
tempfile africasites_all_sim_robust
save `africasites_all_sim_robust', replace

save "${TEMP}/AfricaSites_all_sim_robust.dta", replace


*************************************
** Create Naive Simulation Dataset **
*************************************

*First create yield data

use "${TEMP}/TrialDataset.dta", clear
merge m:1 sitecode using "${TEMP}/TrialSites_gridloc.dta", nogen keep(1 3) keepusing(latgrid longrid)
drop lat lon
g lat = real(latgrid)
g lon = real(longrid)

keep if $toinclude_trial
keep fert yield lat lon

foreach f in 0 1 {
	g yield`f' = yield if fert==`f'
	}
	
collapse (mean) yield0 yield1, by(lat lon)
g yd = yield1 - yield0
replace yield1 = . if yd==.
replace yield0 = . if yd==.
drop yd

merge 1:1 lat lon using "${TEMP}/AfricaSites_fixed.dta", /*nogen*/ keep(2 3) keepusing(aez cellid)
sort cellid
merge 1:1 cellid using `countries', nogen keep(1 3) keepusing(country)
sort cellid 

foreach f in 0 1 {
	bysort country aez: egen yield_1_f`f' = mean(yield`f')
	bysort country aez: egen yield_1_f`f'_sd = sd(yield`f')
	
	bysort country: egen yield_2_f`f' = mean(yield`f')
	bysort country: egen yield_2_f`f'_sd = sd(yield`f')

	bysort aez: egen yield_3_f`f' = mean(yield`f')
	bysort aez: egen yield_3_f`f'_sd = sd(yield`f')
	
	egen yield_4_f`f' = mean(yield`f')
	egen yield_4_f`f'_sd = sd(yield`f')

	g yield_f`f' = .
	g yield_f`f'_sd = .
	
	forvalues i = 1/4 {
		replace yield_f`f' = yield_`i'_f`f' if yield_f`f'==. & yield_`i'_f`f' != .
		replace yield_f`f'_sd = yield_`i'_f`f'_sd if yield_f`f'_sd==. & yield_`i'_f`f'_sd != .
		}
	}

g simyr=2050

sort cellid simyr
isid cellid simyr

save "${TEMP}/AfricaSites_meanyield.dta", replace


*Then create site-specific yield data
clear
set obs $simreps
g iter=_n
g simyr = 2050
tab simyr
tempfile simdraws_naive
save `simdraws_naive', replace

use `africasites_fixed', clear
expand $simreps
g count=1
bysort cellid: g iter=sum(count)
drop count

sort iter
merge m:1 iter using `simdraws_naive', nogen keepusing(simyr)
merge m:m cellid simyr using "${TEMP}/AfricaSites_meanyield.dta", nogen keep(1 3) keepusing(yield_f*)
	
cap drop var
	
tempfile africasites_all_sim_naive
save `africasites_all_sim_naive', replace
	
save "${TEMP}/AfricaSites_all_sim_naive.dta", replace



cap log close