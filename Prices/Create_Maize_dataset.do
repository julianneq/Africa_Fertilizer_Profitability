/******************************************************************************/
/* FILENAME: Create_Maize_dataset.do*/
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

gl ANTE ${MAIZE}/antecedent_files  
gl RAW ${MAIZE}/raw_data
gl CLEAN ${MAIZE}/clean_data
gl SEMI ${MAIZE}/semi_data
gl DO ${MAIZE}/do_files
gl LOG ${MAIZE}/logs

log using "${LOG}/Maize_Dataset `startyear' `startmonth' `startday' `starttime'.log", replace text

* Bring in data
import excel "${ANTE}/FAO_GIEWS_Food_Prices_All.xlsx", sheet("Sheet1") firstrow clear
split DateMonthly

* Various editing, rename variable, delete commas, change _ to space, and make variable numeric instead of string
rename DateMonthly1 country
replace country=subinstr(country, ",", "", 1)
replace country=subinstr(country, "_", " ", 2)
la var country "Country"

rename DateMonthly2 market_type
replace market_type=subinstr(market_type, ",", "", 1)
la var market_type "Retail or Wholesale Market"

rename DateMonthly3 market_location
replace market_location =subinstr(market_location, ",", "", 1)
replace market_location =subinstr(market_location, "_", " ", 2)
la var market_location "Market location"

* commodity requires some extra changes
/* Unfortuantely, there are some oddities in the naming of this variable and the underlying commodity codes 

                |                FAO GIEWS commodity code
 Commodity type |      1005    1005_10    1005_11    1005_12    1005_36 |     Total
----------------+-------------------------------------------------------+----------
          Maize |        54          0          0          4          0 |        58 
  Maize_(white) |         6          8          0         31          0 |        45 
 Maize_(yellow) |         0          0          1          0          0 |         1 
     Maize_meal |         0          0          0          0         15 |        15 
----------------+-------------------------------------------------------+----------
          Total |        60          8          1         35         15 |       119 



Based on documentation, it seems most of this is actually white maize.

90% of the world’s production is yellow maize, but in Africa, 90% of the total maize production is white maize. 
South Africa is currently the main maize producer of the African continent, and almost half of its production 
consists of white maize meant for human consumption. 
Taken from: http://www.vib.be/en/about-vib/Documents/VIB_MaizeInAfrica_EN_2017.pdf

Human consumption [of maize] is highest in Malawi, Tanzania, Kenya, Zambia, and Zimbabwe, while the crop has 
lower importance in the remaining countries. Consumers prefer white maize and use yellow maize primarily for 
[animal] feed manufacturing. Taken from: http://www.fao.org/3/CA2155EN/ca2155en.pdf

Ninety percent of white maize consumption is in Africa and Central America. It fetches premium prices in 
Southern Africa where it represents the main staple food. Yellow maize is preferred in most parts of 
South America and the Caribbean. Taken from: https://www.iita.org/cropsnew/maize/
*/

rename DateMonthly4 commodity
replace commodity =subinstr(commodity, ",", "", 1)
replace commodity =subinstr(commodity, "_", " ", 2)
la var commodity "Commodity type"

/* Discussion with Ellen (03Oct2019) based on above sources, we are going to combine Maize 
and Maize (white), and drop the one Maize (yellow) price series from S. Africa */

replace commodity="Maize" if commodity=="Maize (white)"
drop if commodity=="Maize (yellow)"

rename DateMonthly5 unit
la var unit "Price per unit"

rename DateMonthly8 commodity_code
la var commodity_code "FAO GIEWS original commodity code"

drop DateMonthly6 DateMonthly7

* Change key variables to numeric to use as factor variables
*encode market_location, generate(market_location2)
*encode country, generate(country2)
*encode commodity, generate(commodity2)
*drop market_location country commodity
*rename market_location2 market_location
*rename country2 country
*rename commodity2 commodity

* Merge in latitude and longitude
merge m:m market_location using "${ANTE}/latitude_longitude_All.dta"

* Need to make file with country, market_type, market_location, commodity, commodity_code, unit, latitude, longitude
keep country market_type market_location commodity commodity_code unit latitude longitude
order country market_type market_location commodity commodity_code unit latitude longitude
sort country market_type market_location

* use expand to make duplicate observations, and then export as excel (because 237 months from Jan 2000 to Sep 2019)
expand 237  
sort country market_type market_location


export excel using "${ANTE}/country_market_lat_long_237copies.xls", firstrow(variables) replace

* Make date month file repeat 118 times because there are 118 market_location market_type combinations
* this file will be exported into excel and pasted in excel next to "${SEMI}/country_market_lat_long_237copies.xls"

forvalues x = 1/119 {
	clear
    import excel "${ANTE}/month_date_world_price.xls", sheet("Sheet1") firstrow case(lower)
    save "${ANTE}/month_date_world_price`x'.dta", replace
        }

use "${ANTE}/month_date_world_price1.dta", clear

forvalues x = 2/119 {
    append using "${ANTE}/month_date_world_price`x'.dta"
        }
save "${ANTE}/month_date_world_price119copies.dta", replace

* save this in excel
export excel using "${ANTE}/month_date_world_price_119copies.xls", firstrow(variables) replace

* next take the month_date_world_price_119copies.xls (28204 lines long) and paste it in next to 
* country_market_lat_long_237copies (28204 lines long), then save as ${ANTE}/country_market_lat_long_month_year_world_price.xls
* then paste in each country's price series into that file one by one and save as
* ${ANTE}/country_market_lat_long_month_year_world_local_price.xls, this is the file to start the 
* Maize_Prices_for_Krigging.do and analysis 


