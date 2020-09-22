/******************************************************************************/
/* FILENAME: Maize_Prices_for_Kriging.do*/
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
/* Andrew is 1, change user to 0 and insert your path below */

/* library definitions */

if $user==0 {
gl MAIZE = "Add Your Path Here"
}      

if $user==1 {
gl MAIZE = "/Users/andrewsimons/Google_Drive/Fordham/Research/VCR Paper/Maize"
}  

gl RAW ${MAIZE}/raw_data
gl CLEAN ${MAIZE}/clean_data
gl SEMI ${MAIZE}/semi_data
gl DO ${MAIZE}/do_files
gl LOG ${MAIZE}/logs

log using "${LOG}/Maize_Prices `startyear' `startmonth' `startday' `starttime'.log", replace text

* Bring in data from raw folder
import excel "${RAW}/country_market_lat_long_month_year_world_local_price.xls", sheet("Sheet1") firstrow clear

* regression can not differentiate between maize and maize meal because maize meal countries only have maize meal
* so it can't calculate the difference between types (all folds into country fixed effect)
* therefore we will drop maize meal from the analysis
drop if commodity=="Maize meal"

* need to make certain varaibles categorical (numeric) rather than string
encode market_location, generate(market_location2)
encode country, generate(country2)
encode commodity, generate(commodity2)
encode market_type, generate(market_type2)
drop market_location country commodity market_type
rename market_location2 market_location
rename country2 country
rename commodity2 commodity
rename market_type2 market_type

* Add labels
la var country "Country"
la var market_type "Retail or Wholesale Market"
la var market_location "Market location"
la var commodity "Commodity type"
la var unit "Price per unit"
la var commodity_code "FAO GIEWS original commodity code"
la var latitude "Latitude coordinate"
la var longitude "Longitude coordinate"
la var pricedate "Price on date"
la var month "Month"
la var year "Year"
la var world_price "World price white maize, USA (Kentucky) series, USD/kg"
la var price "Local price white maize, USD/kg"

* create marketid to differentiate markets that have both retail and wholesale price series
egen marketid = group(market_type market_location)
la var marketid "Unique id for market_type market_location"

* create new variables
gen local_div_world=price/world_price
la var local_div_world "Local price in USD divided by world price"

gen log_local_div_world=ln(local_div_world)
la var log_local_div_world "Natural log of local price divided by world price"

gen years_since_initial=year-2000
la var years_since_initial "Years since initial year (2000)"

gen years_since_initial_squared=years_since_initial^2
la var years_since_initial_squared "Years since initial squared"

* save file
save "${SEMI}/Maize_local_world_price_series.dta", replace

* Recover the lowest monthly beta within a country to see which month is lowest price 
* (post harvest flux of commodity onto market) (holding year, year^2 constant)
forvalues x = 1/20 {
	di `x'
    reg log_local_div_world i.market_location i.market_type years_since_initial years_since_initial_squared i.month if country==`x'
	        }
	   
/* Country number, after country is lowest month in maize price in above regression (representing the post-harvest
glut of commodity on market)

		   1 Benin 4
           2 Burundi 4
           3 Cameroon 11
           4 Central African Republic 1
           5 Chad 2
           6 Democratic Republic Congo 2
           7 Ethiopia 2
           8 Ghana 3
           9 Kenya 2
          10 Malawi 6
          11 Mozambique 5
          12 Namibia 2
          13 Nigeria 2
          14 Rwanda 5
          15 Somalia 1
          16 South Africa 6
          17 South Sudan 10
          18 Tanzania 7
          19 Togo 1
          20 Uganda 8  */

* Regression controlling for month/country maize interaction
reg log_local_div_world i.market_location i.market_type years_since_initial years_since_initial_squared country#month

* Margins command for each monthly/country combo for the lowest maize price for that country will give 
* the estimated price wedge for that market location at lowest (peak supply) point for wholesale
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(4) country=(1))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(4) country=(2))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(11) country=(3))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(1) country=(4))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(2) country=(5))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(2) country=(6))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(2) country=(7))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(3) country=(8))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(2) country=(9))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(6) country=(10))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(5) country=(11))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(2) country=(12))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(2) country=(13))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(5) country=(14))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(1) country=(15))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(6) country=(16))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(10) country=(17))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(7) country=(18))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(1) country=(19))
margins i.market_location, at(market_type=2 years_since_initial=19 years_since_initial_squared==361 month=(8) country=(20))

* These margins predict a wedge for each market location assigning the peak supply month (lowest maize price month) for that country
* I could not figure out a way to export these results, so I copied them from the log file into excel and then brought that excel file 
* back into stata with the following
import excel "${RAW}/predicted maize wedge by market location harvest month post estimation raw.xlsx", sheet("Sheet1") firstrow clear

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

* save predicted file
save "${SEMI}/Maize_predicted_beta_by_location_country.dta", replace

* now I need to merge these back into a file with lat and long, country name, etc.
use "${SEMI}/Maize_predicted_beta_by_location_country.dta", clear
merge m:m market_location using "${SEMI}/Maize_local_world_price_series.dta", keepusing(country commodity market_location latitude longitude)

drop _merge
keep in 1/99
label values market_location market_location2

decode country, gen(country_str)

decode market_location, gen(market_location_str)

order country_str market_location_str commodity latitude longitude beta std_error, before(market_location)

export excel country_str market_location_str commodity latitude longitude beta std_error using "${CLEAN}/maize_country_mkt_lat_long_beta_stderror.xls", firstrow(variables) replace


