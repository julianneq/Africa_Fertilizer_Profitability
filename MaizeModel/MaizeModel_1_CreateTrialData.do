***************************************************************************************
***************************************************************************************
********************************* CREATE TRIAL DATA ***********************************
***************************************************************************************
***************************************************************************************

** This file compiles the trial data and combines it with soils and climate data 


set more off
cap log close

loc mlength = 365/12

log using "${TEMP}/MaizeModel_CreateTrialData.log", replace

insheet using "${DATA}/maizedata.lobell.sep2011.csv", clear

foreach var in anthesisdate asi daystosilk anthprecip asi1 postsilkgdd30 presilktmin silktmin post45silktmin presilktmax silktmax post45silktmax presilktmax2 silktmax2 postsilktmax2 presilkgdd830 silkgdd830 postsilkgdd830 presilkgdd30 silkgdd30 {
	replace `var' = "" if `var'=="NA"
	destring `var', replace
	}

encode management, generate(mgmt)

g mgmt_dr = mgmt==1
g mgmt_lown = mgmt==2
g mgmt_lowph=mgmt==3
g mgmt_strk = mgmt==4
g mgmt_opt = mgmt==5

g var_EIHY = vargroup=="EIHY" /* early maturing hybrid*/
la var var_EIHY "Early maturing hybrid"
g var_EPOP = vargroup=="EPOP" /* early maturing opv*/
la var var_EPOP "Early maturing OPV"
g var_ILHY = vargroup=="ILHY" /* intermediate to late maturing hybrid */
la var var_ILHY "Late maturing hybrid"
g var_ILPO = vargroup=="ILPO" /* intermediate to late maturing OPVs */
la var var_ILPO "Late maturing OPV"

g var_hy = var_EIHY==1 | var_ILHY==1
replace var_hy = 0 if var_EPOP==1 | var_ILPO==1
la var var_hy "Hybrid variety"

g yield = exp(logyield)

format plantingdate %08.0f
tostring plantingdate, generate(plntdate)
g plntyr=substr(plntdate,1,4)
destring plntyr, replace
g plntmo=substr(plntdate,5,2)
destring plntmo, replace
g plntda=substr(plntdate,7,2)
destring plntda, replace

replace plntyr = 2000 if yrcode==2000 & plntyr==2005 /*28 obs appear to be miscoded*/
replace plntyr = 2004 if yrcode==2004 & plntyr==2005 /*90 obs appear to be miscoded*/

g date_plant=mdy(plntmo,plntda,plntyr)
g doy_plant=doy(date_plant)
g my_plant=mofd(date_plant)
replace my_plant = my_plant + 1 if day(date_plant) > `mlength'/2

tempfile trialdataset
save `trialdataset', replace


** Now prepare additional data

insheet using "${DATA}/additional.cimmyt.maizedata.jul2020.csv", clear

g mgmt_lown = management=="Low N"
g mgmt_opt = management=="Optimal"

g var_EIHY = vargroup=="EIHY" /* early maturing hybrid*/
la var var_EIHY "Early maturing hybrid"
g var_EPOP = vargroup=="EPOP" /* early maturing opv*/
la var var_EPOP "Early maturing OPV"
g var_ILHY = vargroup=="ILHY" /* intermediate to late maturing hybrid */
la var var_ILHY "Late maturing hybrid"
g var_ILPO = vargroup=="ILPO" /* intermediate to late maturing OPVs */
la var var_ILPO "Late maturing OPV"

g var_hy = var_EIHY==1 | var_ILHY==1
replace var_hy = 0 if var_EPOP==1 | var_ILPO==1
la var var_hy "Hybrid variety"

format plantingdate %08.0f
tostring plantingdate, generate(plntdate)
g plntyr=substr(plntdate,1,4)
destring plntyr, replace
g plntmo=substr(plntdate,5,2)
destring plntmo, replace
g plntda=substr(plntdate,7,2)
destring plntda, replace

g date_plant=mdy(plntmo,plntda,plntyr)
g doy_plant=doy(date_plant)
g my_plant=mofd(date_plant)
replace my_plant = my_plant + 1 if day(date_plant) > `mlength'/2

g yield = exp(logyield)

tempfile trialdataset_addl
save `trialdataset_addl', replace


** Now prepare Wortmann data

import excel using "${DATA}/additional.maizedata.wortmann.aug2021.xlsx", firstrow case(lower) clear

g var_hy = vargroup=="HY"
replace var_hy = 0 if vargrou=="OPV"

g mgmt_lown = management=="Low N"
g mgmt_opt = management=="Optimal"

format plantingdate %08.0f
tostring plantingdate, generate(plntdate)
g plntyr=substr(plntdate,1,4)
destring plntyr, replace
g plntmo=substr(plntdate,5,2)
destring plntmo, replace
g plntda=substr(plntdate,7,2)
destring plntda, replace

g date_plant=mdy(plntmo,plntda,plntyr)
g doy_plant=doy(date_plant)
g my_plant=mofd(date_plant)
replace my_plant = my_plant + 1 if day(date_plant) > `mlength'/2

g yield = exp(logyield)

g wortmann=1
la var wortmann "OFRA trial (Wortmann)"

tempfile trialdataset_wortmann
save `trialdataset_wortmann', replace


** Join the two trial datasets together

use `trialdataset', clear
append using `trialdataset_addl'
g wortmann=0

append using `trialdataset_wortmann'

tempfile trialdataset
save `trialdataset', replace


** Merge with temperature and precip data

* First match each trial site with a grid
insheet using "${DATA}/additional_EIL_site_lat_lon_files.csv", clear
keep locationid file1 
tempfile sites_addl
save `sites_addl', replace

insheet using "${DATA}/WortmannSitesSoilsData.csv", clear
keep locationid file1 
tempfile sites_wortmann
save `sites_wortmann', replace

insheet using "${DATA}/EIL_site_lat_lon_files.csv", clear
append using `sites_addl'
append using `sites_wortmann'

keep locationid file1
rename file1 gridstring
g gridstringlength = strlen(gridstring)
g gridpos1 = strpos(gridstring,"_")
g latlength = gridpos1-4
g gridpos2 = gridpos1+5
g lonlength = gridstringlength-gridpos2-3
g latgrid = substr(gridstring,4,latlength)
g longrid = substr(gridstring,gridpos2,lonlength)

keep locationid latgrid longrid
rename locationid sitecode

sort sitecode
save "${TEMP}/TrialSites_gridloc.dta", replace


* Get monthly precip aggregates for each site-planting date combo

insheet using "${DATA}/additional_EIL_site_lat_lon_files.csv", clear
drop file3 file4
tempfile filematch_addl
save `filematch_addl'

insheet using "${DATA}/WortmannSitesSoilsData.csv", clear
keep locationid country location region lat lon elevm file1 file2
tempfile filematch_wortmann
save `filematch_wortmann', replace

insheet using "${DATA}/EIL_site_lat_lon_files.csv", clear
append using `filematch_addl'
append using `filematch_wortmann'

rename locationid sitecode
rename latitude lat
rename longitude lon

format %9.3f lat lon
sort sitecode
tempfile filerefs
save `filerefs', replace

* drop sites not in study
insheet using "${DATA}/Site-years-withadditional.csv", clear
rename locationid sitecode
tempfile site_years_temp
save `site_years_temp', replace

use `filerefs', clear
merge 1:1 sitecode using `site_years_temp', nogen keepusing(instudy) keep(1 3)
drop if instudy==0
drop instudy

tempfile site_years
save `site_years', replace


use `trialdataset', clear
merge m:1 sitecode using `site_years', nogen keep(1 3)
g count = 1
collapse (sum) count, by(sitecode file1 file2 file3 file4 date_plant)
tempfile site_pldate
save `site_pldate', replace

loc firstsite = 1
levelsof sitecode, local(sites)

loc missingfiles=0
foreach site of local sites {
	display "`site'"
	use `site_pldate', clear
	keep if sitecode==`site'
	
	loc filecount = 0
	forvalues i = 1/4 {
		levelsof file`i', local(file`i') clean	
		if strlen("`file`i''") >= 1 {
			loc filecount = `filecount' + 1
			}
		}

	levelsof date_plant, local(dates)

	forvalues i = 1/`filecount' {
		*USE CRU FOR CIMMYT AND WORTMANN SITES
		display "File: ${CLIMDATA}/TempPrecipCRU/`file`i''"
		cap noisily import delimited using  "${CLIMDATA}/TempPrecipCRU/`file`i''", delimiters(" ", collapse) varnames(1) stringcols(1) 		numericcols(2 3 4 5) clear
		if _rc == 0 {
			g yr = real(substr(date,1,2))
			replace yr = yr + 100 if yr < 79
			replace yr = yr + 1900
			g mdy = mdy(1,1,yr) - 1 + real(substr(date,3,.))
			drop date
			format mdy %d
			rename mdy date
			rename tavg temp 
	
			recode temp tmin tmax rain (-99=.)
	
			g filename = `i'

			tempfile file_`i'
			save `file_`i'', replace
			}
			
			else if _rc != 0 {
				loc missingfiles = `missingfiles'+1
				}
	}
				
	use `file_1', clear
	if `filecount' > 1 {
		forvalues i = 2/`filecount' {
			cap append using `file_`i''
			}
		}
	
	collapse (mean) temp tmin tmax rain, by(date)
	
	foreach date of local dates {
			g gmonth = 0
			forvalues i = 1/5 {
				replace gmonth = `i' if (date >= int(`date'+(`i'-1)*`mlength')) & (date < int(`date'+`i'*`mlength'))
				}
			forvalues i = 1/5 {
				g precip_m`i'_`date' = rain * (gmonth==`i')
				g temp_m`i'_`date' = temp * (gmonth==`i')
				g gdd_m`i'_`date' = (min(0.5*(tmin+tmax),${TM}) - ${TB})*(gmonth==`i')
				g hdd_m`i'_`date' = (max(tmax,${TM}) - ${TM})*(gmonth==`i')
				g ztemp_m`i'_`date' = temp==0 * (gmonth==`i')				
				}
			drop gmonth
			}
	
	recode temp_m* gdd_m* /*hdd_m**/ (0=.)
	collapse (sum) precip_m1_* precip_m2_* precip_m3_* precip_m4_* precip_m5_* ztemp_m1_* ztemp_m2_* ztemp_m3_* ztemp_m4_* ztemp_m5_* gdd_m1_* gdd_m2_* gdd_m3_* gdd_m4_* gdd_m5_* hdd_m1_* hdd_m2_* hdd_m3_* hdd_m4_* hdd_m5_* (mean) temp_m1_* temp_m2_* temp_m3_* temp_m4_* temp_m5_* 
	g sitecode = `site'
	reshape long precip_m1_ precip_m2_ precip_m3_ precip_m4_ precip_m5_ temp_m1_ temp_m2_ temp_m3_ temp_m4_ temp_m5_ gdd_m1_ gdd_m2_ gdd_m3_ gdd_m4_ gdd_m5_ hdd_m1_ hdd_m2_ hdd_m3_ hdd_m4_ hdd_m5_ ztemp_m1_ ztemp_m2_ ztemp_m3_ ztemp_m4_ ztemp_m5_ , i(site) j(date)
	rename *_ *
	rename date date_plant

	if `firstsite'==1 {
		tempfile site_pldate_tempprecip
		save `site_pldate_tempprecip', replace
		}
	
		else if `firstsite' == 0 {
			append using `site_pldate_tempprecip'
			tempfile site_pldate_precip
			save `site_pldate_tempprecip', replace
			}

	loc firstsite=0	
	}	

display "Number of missing files: `missingfiles'"
su ztemp_m*, detail

forvalues i = 1/5 {
	tab sitecode if ztemp_m`i' > 0
	}
	
drop ztemp_m*

use `site_pldate_tempprecip', clear
save "${TEMP}/Site_pldate_tempprecip", replace


** Merge trial data with AEZ data 

* There is no aez_join for the original sites so the next 16 lines are a workaround
insheet using "${DATA}/AfricaDatabase.csv", clear
rename aezid aez 
keep lat lon aez
tempfile aezmatch
save `aezmatch', replace

use `sites_addl', clear
keep locationid
rename locationid sitecode
merge 1:1 sitecode using "${TEMP}/TrialSites_gridloc.dta", nogen keep(3) keepusing(latgrid longrid)
g lat = real(latgrid) 
g lon = real(longrid)

merge m:1 lat lon using `aezmatch', nogen keep(1 3)
rename sitecode locationid
keep locationid aez 
tempfile aez_addl
save `aez_addl', replace

insheet using "${DATA}/sites_aez_join.csv", clear
keep locationid  gridcode
rename gridcode aez

tempfile aez_join
save `aez_join', replace

insheet using "${DATA}/WortmannSitesAEZ.csv", clear
keep locationid aezcode
rename aezcode aez

append using `aez_join'
append using `aez_addl'

rename locationid sitecode

keep sitecode aez
label define aezs 0 "Subtropic - warm / arid" 1 "Subtropic - warm / semiarid" 2 "Subtropic - warm / subhumid" /*
*/3 "Subtropic - warm / humid" 4 "Subtropic - cool / arid" 5 "Subtropic - cool / semiarid" 6 "Subtropic - cool / subhumid" /*
*/7 "Tropic - warm / arid" 8 "Tropic - warm / semiarid" 9 "Tropic - warm / subhumid" 10 "Tropic - warm / humid" /*
*/11 "Tropic - cool /arid" 12 "Tropic - cool / semiarid" 13 "Tropic - cool / subhumid" 14 "Tropic - cool / humid" 
la values aez aezs

tempfile sites_aez
save `sites_aez', replace

use `trialdataset', clear
merge m:1 sitecode using `sites_aez', nogen keep(1 3) keepusing(aez)

tempfile trialdataset
save `trialdataset', replace


** Merge trial data with temperature and precip data

use `trialdataset', clear

sort sitecode date_plant
merge m:1 sitecode date_plant using "${TEMP}/Site_pldate_tempprecip.dta", nogen keep(1 3)

egen temp_p1 = rowmean(temp_m1 temp_m2)
egen temp_p2 = rowmean(temp_m3)
egen temp_p3 = rowmean(temp_m4 temp_m5)
recode temp_p1 temp_p2 temp_p3 (0=.) /*Mean temp can't feasibly be zero degrees celcius - this is no observations*/

egen gdd_p1 = rowtotal(gdd_m1 gdd_m2)
egen gdd_p2 = rowtotal(gdd_m3)
egen gdd_p3 = rowtotal(gdd_m4 gdd_m5)

egen hdd_p1 = rowtotal(hdd_m1 hdd_m2)
egen hdd_p2 = rowtotal(hdd_m3)
egen hdd_p3 = rowtotal(hdd_m4 hdd_m5)

egen precip_p1 = rowtotal(precip_m1 precip_m2)
egen precip_p2 = rowtotal(precip_m3)
egen precip_p3 = rowtotal(precip_m4 precip_m5)

egen precip_tot = rowtotal(precip_p1 precip_p2 precip_p3)
tab sitecode date_plant if precip_tot==0

foreach var in precip_p1 precip_p2 precip_p3 { /*5 months of zero growing season precip suggests missing data not extreme drought - this is 24 observations*/
	replace `var' = . if precip_tot==0
	}

** Additional data manipulations

g fert = (mgmt_lown == 0)

sort sitecode date_plant

egen site_year_id = group(sitecode plntyr)

tempfile trialdataset
save `trialdataset', replace


** Merge with soils data

insheet using "${DATA}/AdditionalSitesSoilsData.csv", clear
drop file1 file2 file3 file4
tempfile sitesoils_addl
save `sitesoils_addl', replace

insheet using "${DATA}/WortmannSitesSoilsData.csv", clear
drop file1 file2 file3 file4 
tempfile sitesoils_wortmann
save `sitesoils_wortmann', replace

insheet using "${DATA}/SitesSoilsData.csv", clear
drop file1 file2 file3 file4
append using `sitesoils_addl'
append using `sitesoils_wortmann'

rename locationid sitecode

loc soilvarlist bld cec clyppt crfvol drainfao eackcl nto orcdrc phihox sltppt sndppt elevm 

foreach soilvar in `soilvarlist' {
	recode `soilvar' (-9999=.)
	}

keep sitecode `soilvarlist' latitude longitude country elevm

tempfile sitessoilsdata
save `sitessoilsdata', replace


** Merge with trials data

use `trialdataset', clear
drop latitude longitude country elevm 
merge m:1 sitecode using `filerefs', nogen keep(1 3) keepusing(lat lon country elevm)
merge m:1 sitecode using `sitessoilsdata', nogen keep(1 3)

rename bld bulkdens 
recode bulkdens 0=. /*lots of 0 values and this makes no sense*/
replace bulkdens = bulkdens/1000 /*See ISRC report*/

rename cec soilcec /*there are some zero values but this seems feasible according to ISRIC report, values seem well scaled according to ISRIC report*/

rename clyppt claypct
rename sltppt siltpct
rename sndppt sandpct
rename crfvol coarsefrag
foreach soilvar in claypct siltpct sandpct coarsefrag {
	replace `soilvar' = `soilvar' / 100
	}
rename drainfao soildrain
g poordrain = soildrain==1 | soildrain==2 | soildrain==3
replace poordrain=. if soildrain==0
/* 1:V, 2:P, 3:I, 4:M, 5:W, 6:S, 7:E from http://www.isric.org/sites/default/files/ISRIC_Report_%202015_02.pdf */

rename eackcl acidity /*there are some zero values but this seems feasible according to data cleaning protocol from the ISRIC report. Scaling consistent also */

rename nto soiln /*range is 0 to 3.5, this seems in line with ISRIC */
rename orcdrc soilom 

rename phihox soilph /*range is 54 to 74, it should be 5.4 to 7.4 */
recode soilph (0=.)
replace soilph = soilph / 10


drop ztemp_*
 
loc templab "Temp"
loc tempunit "(mean, degrees C)"
loc gddlab "GDDs"
loc gddunit "(sum, degrees C)"
loc hddlab "HDDs"
loc hddunit "(sum, degrees C)"
loc preciplab "Precip"
loc precipunit "(tot, mm)"

loc p1per months 1-2
loc p2per month 3
loc p3per months 4-5
foreach var in temp hdd gdd precip {
   forvalues m = 1/5 {
		la var `var'_m`m' "``var'lab', month `m' ``var'unit'"
		}
	forvalues p = 1/3 {
		la var `var'_p`p' "``var'lab', `p`p'per' ``var'unit'"
		}
	}
la var yield "Yield (t/ha)"
la var logyield "Log yield (t/ha)"
la var fert "Fertilized plot (dummy)"
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
la var elevm "Elevation (M)"
la var wortmann "OFRA trial (Wortmann)"

drop location region sitecode_orig latitude longitude notes

cap drop count
g count = mgmt_lown==1
bysort sitecode: egen sitehaslown=sum(count)
cap drop count
replace sitehaslown = (sitehaslown >= 1 & sitehaslown != .)

sort sitecode
merge m:1 sitecode using `site_years', nogen keep(1 3) keepusing(country)

replace country = "KENYA" if country=="Kenya"
replace country = "MALAWI" if country=="Malawi"
replace country = "TANZANIA" if country=="Tanzania"
replace country = "DEMOCRATIC REPUBLIC OF THE CONGO" if country=="ZAIRE"
replace country = "RWANDA" if country == "Rwanda"
replace country = "SOUTH AFRICA REP." if country=="South Africa"
replace country = "ZAMBIA" if country=="Zambia"

tempfile trialdataset
save `trialdataset', replace

save "${TEMP}/TrialDataset", replace

cap log close

