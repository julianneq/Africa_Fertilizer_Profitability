/******************************************************************************/
/* FILENAME: Urea_Prices_for_Kriging.do*/
/******************************************************************************/

clear all
capture log close /* kill any open logs */
set more off

/* create time-date stamp for unique log file name */

local startdate "$S_DATE"
local starttime "$S_TIME"

input var1
1
end

gen starttime=subinstr("`starttime'",":","H",1)
replace starttime=subinstr(starttime,":","M",1)
gen startday=word("`startdate'",1)
gen startmonth=word("`startdate'",2)
gen startyear=word("`startdate'",3)

local starttime=starttime[1]
local startday=startday[1]
local startmonth=startmonth[1]
local startyear=startyear[1]

/* user definitions */

gl user=1
/* Andrew is 1 */

/* library definitions */

if $user==0 {
gl UREA = "Add Your Path Here"
}      

if $user==1 {
gl UREA = "/Users/andrewsimons/Google_Drive/Fordham/Research/VCR Paper/Fertilizer"
}  

gl RAW ${UREA}/raw_data
gl CLEAN ${UREA}/clean_data
gl SEMI ${UREA}/semi_data
gl DO ${UREA}/do_files
gl LOG ${UREA}/logs

log using "${LOG}/Fertilizer_Prices `startyear' `startmonth' `startday' `starttime'.log", replace text

* Bring in data
import excel "${RAW}/Urea_local_world_price_series.xls", sheet("Sheet1") firstrow case(lower) clear

* Change key variables to numeric to use as factor variables
encode market_location, generate(market_location2)
encode country, generate(country2)
drop market_location
drop country
rename market_location2 market_location
rename country2 country

save "${SEMI}/Urea_local_world_price_series.dta", replace

* Recover highest monthly beta within a country to see which month is highest price (peak fertilizer demand) (holding year, year^2 constant)
forvalues x = 1/17 {
	di `x'
    reg log_local_div_world i.market_location years_since_initial years_since_initial_squared i.month if country==`x'
	        }
	       
/* Country number, after country is peak month in fertilizer price in above regression 
          
           1 Benin 3 
           2 Burkina Faso month 5  
           3 Burundi month 8
           4 Cote d'Ivoire month 6
           5 Ghana month 3
           6 Kenya month 2
           7 Malawi month 5
           8 Mali month 5
           9 Mozambique 6
          10 Niger 12
          11 Nigeria month 4
          12 Rwanda month 8 
          13 Senegal month 5 
          14 Tanzania month 8
          15 Togo 5
          16 Uganda month 4
          17 Zambia month 6  */

* This holds month at 1 for all, but we want different peak month for each country
reg log_local_div_world i.market_location years_since_initial years_since_initial_squared i.month
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=1)

* Regression controlling for month/country fertilizer interaction
reg log_local_div_world i.market_location years_since_initial years_since_initial_squared country#month

* Margins command for each monthly/country combo for the highest fertilizer price for that country will give 
* the estimated price wedge for that market location at highest (peak demand) point
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=(3) country=(1))
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=(5) country=(2))
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=(8) country=(3))
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=(6) country=(4))
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=(3) country=(5))
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=(2) country=(6))
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=(5) country=(7))
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=(5) country=(8))
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=(6) country=(9))
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=(12) country=(10))
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=(4) country=(11))
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=(8) country=(12))
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=(5) country=(13))
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=(8) country=(14))
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=(5) country=(15))
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=(4) country=(16))
margins i.market_location, at(years_since_initial=9 years_since_initial_squared==81 month=(6) country=(17))

* These margins predict a wedge for each market location assigning the peak demand month (highest fertilzier price month) for that country
* I could not figure out a way to export these results, so I copied them from the log file into excel and then brought that excel file back 
* into stata with the following
import excel "${RAW}/predicted urea wedge by market location harvest month post estimation raw.xlsx", sheet("Sheet1") firstrow clear

* All responses just in one column, so split into difference columns on each space
split market_location

* Now drop all the (no estimate) responses
drop if market_location4=="(not"

* no duplicates
duplicates report

* keep and rename variables
keep market_location1 market_location3 market_location4
rename market_location1 market_location
rename market_location3 beta
rename market_location4 std_error
la var market_location "Name of market"
la var beta "Predicted percent increase over world price at high demand month in 2019"
la var std_error "Standard error of beta"
destring market_location, replace
destring beta, replace
destring std_error, replace

save "${SEMI}/Urea_predicted_beta_by_location_country.dta", replace

* now I need to merge these back into a file with lat and long, country name, etc.
use "${SEMI}/Urea_predicted_beta_by_location_country.dta", clear
merge m:m market_location using "${SEMI}/Urea_local_world_price_series.dta", keepusing(country commodity market_location latitude longitude)

drop _merge
keep in 1/102
label values market_location market_location2

* Switch back to string variables for excel export
decode country, gen(country_str)

decode market_location, gen(market_location_str)

* Order data so it looks nice
order country_str market_location_str commodity latitude longitude beta std_error, before(market_location)

* Export to excel
export excel country_str market_location_str commodity latitude longitude beta std_error using "${CLEAN}/urea_country_mkt_lat_long_beta_stderror.xls", firstrow(variables) replace


