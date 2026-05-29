// Sub-group analyses
	// cancer-specific mortality
		// restricting entire cohort:
			// on tumour histology
			// BMI (<25 & >=25)
			// age (<70 and 70+ years old)
			// Stages 1, 2, 3 
			// high glucocorticoid use
			// prior osteoporosis diagnosis
		
********************************************************************************		
**# Endometrial cancer
********************************************************************************
			
clear
capture postutil clear
macro drop _all
set trace off

local model1 "" 

local model2 "ageatdx"

local model3 "ageatdx yeardx ncras_chemo6mafterdx ncras_surgery6mafterdx ncras_radio6mafterdx cc_mipriordx cc_chfpriordx cc_pvdpriordx cc_dempriordx cc_cpdpriordx cc_rheupriordx cc_pudpriordx cc_mldpriordx cc_sldpriordx cc_diabpriordx cc_diabcomppriordx cc_plegiapriordx cc_renalpriordx cc_cerebrvdpriordx deprivation_quintile"	

local model4 "ageatdx yeardx ncras_chemo6mafterdx ncras_surgery6mafterdx ncras_radio6mafterdx cc_mipriordx cc_chfpriordx cc_pvdpriordx cc_dempriordx cc_cpdpriordx cc_rheupriordx cc_pudpriordx cc_mldpriordx cc_sldpriordx cc_diabpriordx cc_diabcomppriordx cc_plegiapriordx cc_renalpriordx cc_cerebrvdpriordx deprivation_quintile i.stage"	


postfile results str10 Database str50 Cancer str20 Mortality str40 Subgroup str40 Analysis Outcomes PersonYears HR lCI UCI Pvalue HR2 lCI2 UCI2 Pvalue2 HR3 lCI3 UCI3 Pvalue3 HR4 lCI4 UCI4 Pvalue4 using "D:\CPRD data\Bisphos results\Both\cancerspecific_subgroups_ec.dta",replace	

********************************************************************************
** subgroups 
********************************************************************************

local sub1 "ncras_endometrioid==1"
local sub2 "ncras_nonendometrioid==1"		
local sub3 "latestbmipriordx<25" 
local sub4 "latestbmipriordx>=25 & latestbmipriordx~=."
local sub7 "ageatdx<70"
local sub8 "ageatdx>=70"
local sub9 "stage==1"
local sub10 "stage==2"
local sub11 "stage==3"
local sub11 "stage==3"
local sub12 "gp_m_corticosteroidprioryrddd>=63 & gp_m_corticosteroidprioryrddd~=."
local sub13 "osteopriordxyr==1"

********************************************************************************	

foreach c in 1  {

foreach s in 1 2 3 4 7 8 9 10 11 12 13 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c' & `sub`s''

* exclude duplicate GP practices
local characteristics="-Excluded: duplicate GP practices"
drop if inlist(pracid, 20024,20036, 20091,20171,20178,20202,20254,20389, 20430,20452, 20469,20487,20552,20554,20640,20717, 20734,20737,20740,20790,20803, 20822,20868,20912,20996, 21001,21015,21078,21112,21118,21172 ,21173, 21277,21281,21331,21334,21390,21430 ,21444,21451,21529,21553,21558,21585 ) & db==1

** cancer diagnosis between 1st April 2003 (start of outpatient) and 31st december 2018 (end of cancer reg data)
drop if ncras_diagnosisdatebest<mdy(4,1,2003) | ncras_diagnosisdatebest>mdy(12,31,2018)

* gold lesss than 1 year uts prior
gen startminusyear=start-365
drop if gp_uts>startminusyear & gp_uts~=.

// exclusion criteria

* drop <55 
drop if ageatdx<55

* drop stage 4
drop if stage==4

* drop missing stage
drop if stage==.

// set up dates 
gen lag=182.5
gen date_entry=start+lag
gen date_endfup=min(ons_deathdate,mdy(11,16,2020), gp_lastcollectiondate, gp_transferoutdate, gp_endrecords, gp_regenddate)

* less than 6 months of follow-up
drop if date_entry>date_endfup

* registered with gp less than 1 year 
drop if gp_regstartdate>startminusyear

gen date_outcome=.
replace date_outcome=ons_deathdate if ons_pc_ec==1 & gynaecancer==1
replace date_outcome=ons_deathdate if ons_pc_oc==1 & gynaecancer==2
gen outcome=1 if date_outcome~=.

* check if outcome occurs after end of fup
replace outcome=. if date_outcome>date_endfup

count 

// stset data 
stset date_endfup, fail(outcome=1) enter(time date_entry) id(patid) origin(time date_entry) scale(365.25) 

** stsplit to create time-varying exposure definition

* lag dates by six months first 
foreach i in bisphos  {
	replace gp_m_`i'afterdxdate=gp_m_`i'afterdxdate+182.5 if gp_m_`i'afterdxdate~=.
	replace gp_m_`i'afterdxdate=. if gp_m_`i'afterdxdate>=date_endfup
}

* then split 
foreach i in bisphos  {
	stsplit newtime`i', at(0) after(time=gp_m_`i'afterdxdate)
	
	gen ever`i'=1 if _t0>=(gp_m_`i'afterdxdate-_origin)/365.25 & gp_m_`i'afterdxdate~=.
	replace ever`i'=0 if ever`i'==.
}

* make variable for each time period
gen bisphos=1 if everbisphos==0
replace bisphos=2 if everbisphos==1
lab val bisphos bisphos 
lab define bisphos 1"Non-use" 2"Any use"

levelsof bisphos, local(bisphos)
	foreach i in `bisphos' {
	
		* database
		local database="CPRD"

		* cancer
		local cancer: label gynaecancer `c'

		* mortality
		local mortality="cancer-specific"

		* Subgroup type 
		local subgroup "`sub`s''"
		
		* Analysis 
		local analysis:label bisphos `i'

		* Outcomes 
		count if outcome==1 & bisphos==`i'
		local outcomes=r(N)
		
		* person time
		stptime if bisphos==`i'
		local ptime=r(ptime)

		* unadusted HR
		stcox i.bisphos
		matrix results=r(table)
		local hr=results[1,(2)]
		local lci=results[5,2]
		local uci=results[6,2]
		local pvalue=results[4,2]

		* age adjusted
		stcox i.bisphos `model2'
		matrix results=r(table)
		local 2hr=results[1,2]
		local 2lci=results[5,2]
		local 2uci=results[6,2]
		local 2pvalue=results[4,2]	
			
		stcox i.bisphos `model3'
		matrix results=r(table)
		local 3hr=results[1,2]
		local 3lci=results[5,2]
		local 3uci=results[6,2]
		local 3pvalue=results[4,2]	
			
		stcox i.bisphos `model4'
		matrix results=r(table)
		local 4hr=results[1,2]
		local 4lci=results[5,2]
		local 4uci=results[6,2]
		local 4pvalue=results[4,2]	

		post results ("`database'") ("`cancer'") ("`mortality'") ("`subgroup'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
}
}

postclose results	

clear
use "D:\CPRD data\Bisphos results\Both\cancerspecific_subgroups_ec.dta"

** formatting
gen unadjusted=string(HR, "%9.2f")+" ("+string(lCI,"%9.2f")+"-"+string(UCI,"%9.2f")+")"
gen adjusted2=string(HR2, "%9.2f")+" ("+string(lCI2,"%9.2f")+"-"+string(UCI2,"%9.2f")+")"
gen adjusted3=string(HR3, "%9.2f")+" ("+string(lCI3,"%9.2f")+"-"+string(UCI3,"%9.2f")+")"
gen adjusted4=string(HR4, "%9.2f")+" ("+string(lCI4,"%9.2f")+"-"+string(UCI4,"%9.2f")+")"

foreach i in unadjusted adjusted2 adjusted3 adjusted4 {
	replace `i'="ref" if Analysis=="Non-use"
}

keep Database Cancer Mortality Subgroup Analysis Outcomes PersonYears unadjusted adjusted2 adjusted3 adjusted4
order Database Cancer Mortality Subgroup Analysis Outcomes PersonYears unadjusted adjusted2 adjusted3 adjusted4

** keep track of covariates
gen model1="`model1'"
gen model2="`model2'"
gen model3="`model3'"
gen model4="`model4'"

save "D:\CPRD data\Bisphos results\Both\cancerspecific_subgroups_ec.dta", replace


********************************************************************************
**# Ovarian cancer
********************************************************************************
			
clear
capture postutil clear
macro drop _all
set trace off

local model1 "" 

local model2 "ageatdx"

local model3 "ageatdx yeardx ncras_chemo6mafterdx ncras_surgery6mafterdx ncras_radio6mafterdx cc_mipriordx cc_chfpriordx cc_pvdpriordx cc_dempriordx cc_cpdpriordx cc_rheupriordx cc_pudpriordx cc_mldpriordx cc_sldpriordx cc_diabpriordx cc_diabcomppriordx cc_plegiapriordx cc_renalpriordx cc_cerebrvdpriordx deprivation_quintile"	

local model4 "ageatdx yeardx ncras_chemo6mafterdx ncras_surgery6mafterdx ncras_radio6mafterdx cc_mipriordx cc_chfpriordx cc_pvdpriordx cc_dempriordx cc_cpdpriordx cc_rheupriordx cc_pudpriordx cc_mldpriordx cc_sldpriordx cc_diabpriordx cc_diabcomppriordx cc_plegiapriordx cc_renalpriordx cc_cerebrvdpriordx deprivation_quintile i.stage"	


postfile results str10 Database str50 Cancer str20 Mortality str40 Subgroup str40 Analysis Outcomes PersonYears HR lCI UCI Pvalue HR2 lCI2 UCI2 Pvalue2 HR3 lCI3 UCI3 Pvalue3 HR4 lCI4 UCI4 Pvalue4 using "D:\CPRD data\Bisphos results\Both\cancerspecific_subgroups_oc.dta",replace	

********************************************************************************
** subgroups 
********************************************************************************

local sub1 "ncras_epithelial==1"
local sub2 "ncras_borderline==1"	
local sub22 "ncras_nonepithelial==1"			
local sub3 "latestbmipriordx<25"
local sub4 "latestbmipriordx>=25 & latestbmipriordx~=."
local sub7 "ageatdx<70"
local sub8 "ageatdx>=70"
local sub9 "stage==1"
local sub10 "stage==2"
local sub11 "stage==3"
local sub12 "gp_m_corticosteroidprioryrddd>=63 & gp_m_corticosteroidprioryrddd~=."
local sub13 "osteopriordxyr==1"

********************************************************************************

foreach c in 2  {

foreach s in 1 2 3 4 7 8 9 10 11 12 13 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c' & `sub`s''

* exclude duplicate GP practices
local characteristics="-Excluded: duplicate GP practices"
drop if inlist(pracid, 20024,20036, 20091,20171,20178,20202,20254,20389, 20430,20452, 20469,20487,20552,20554,20640,20717, 20734,20737,20740,20790,20803, 20822,20868,20912,20996, 21001,21015,21078,21112,21118,21172 ,21173, 21277,21281,21331,21334,21390,21430 ,21444,21451,21529,21553,21558,21585 ) & db==1

** cancer diagnosis between 1st April 2003 (start of outpatient) and 31st december 2018 (end of cancer reg data)
drop if ncras_diagnosisdatebest<mdy(4,1,2003) | ncras_diagnosisdatebest>mdy(12,31,2018)

* gold lesss than 1 year uts prior
gen startminusyear=start-365
drop if gp_uts>startminusyear & gp_uts~=.

// exclusion criteria

* drop <55 
drop if ageatdx<55

* drop stage 4
drop if stage==4

* drop missing stage
drop if stage==.

// set up dates 
gen lag=182.5
gen date_entry=start+lag
gen date_endfup=min(ons_deathdate,mdy(11,16,2020), gp_lastcollectiondate, gp_transferoutdate, gp_endrecords, gp_regenddate)

* less than 6 months of follow-up
drop if date_entry>date_endfup

* registered with gp less than 1 year 
drop if gp_regstartdate>startminusyear

gen date_outcome=.
replace date_outcome=ons_deathdate if ons_pc_ec==1 & gynaecancer==1
replace date_outcome=ons_deathdate if ons_pc_oc==1 & gynaecancer==2
gen outcome=1 if date_outcome~=.

* check if outcome occurs after end of fup
replace outcome=. if date_outcome>date_endfup

count 

// stset data 
stset date_endfup, fail(outcome=1) enter(time date_entry) id(patid) origin(time date_entry) scale(365.25) 

** stsplit to create time-varying exposure definition

* lag dates by six months first 
foreach i in bisphos  {
	replace gp_m_`i'afterdxdate=gp_m_`i'afterdxdate+182.5 if gp_m_`i'afterdxdate~=.
	replace gp_m_`i'afterdxdate=. if gp_m_`i'afterdxdate>=date_endfup
}

* then split 
foreach i in bisphos  {
	stsplit newtime`i', at(0) after(time=gp_m_`i'afterdxdate)
	
	gen ever`i'=1 if _t0>=(gp_m_`i'afterdxdate-_origin)/365.25 & gp_m_`i'afterdxdate~=.
	replace ever`i'=0 if ever`i'==.
}

* make variable for each time period
gen bisphos=1 if everbisphos==0
replace bisphos=2 if everbisphos==1
lab val bisphos bisphos 
lab define bisphos 1"Non-use" 2"Any use"

levelsof bisphos, local(bisphos)
	foreach i in `bisphos' {
	
		* database
		local database="CPRD"

		* cancer
		local cancer: label gynaecancer `c'

		* mortality
		local mortality="cancer-specific"

		* Subgroup type 
		local subgroup "`sub`s''"
		
		* Analysis 
		local analysis:label bisphos `i'

		* Outcomes 
		count if outcome==1 & bisphos==`i'
		local outcomes=r(N)
		
		* person time
		stptime if bisphos==`i'
		local ptime=r(ptime)

		* unadusted HR
		stcox i.bisphos
		matrix results=r(table)
		local hr=results[1,(2)]
		local lci=results[5,2]
		local uci=results[6,2]
		local pvalue=results[4,2]

		* age adjusted
		stcox i.bisphos `model2'
		matrix results=r(table)
		local 2hr=results[1,2]
		local 2lci=results[5,2]
		local 2uci=results[6,2]
		local 2pvalue=results[4,2]	
			
		stcox i.bisphos `model3'
		matrix results=r(table)
		local 3hr=results[1,2]
		local 3lci=results[5,2]
		local 3uci=results[6,2]
		local 3pvalue=results[4,2]	
			
		stcox i.bisphos `model4'
		matrix results=r(table)
		local 4hr=results[1,2]
		local 4lci=results[5,2]
		local 4uci=results[6,2]
		local 4pvalue=results[4,2]	

		post results ("`database'") ("`cancer'") ("`mortality'") ("`subgroup'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
}
}

postclose results	

clear
use "D:\CPRD data\Bisphos results\Both\cancerspecific_subgroups_oc.dta"

** formatting
gen unadjusted=string(HR, "%9.2f")+" ("+string(lCI,"%9.2f")+"-"+string(UCI,"%9.2f")+")"
gen adjusted2=string(HR2, "%9.2f")+" ("+string(lCI2,"%9.2f")+"-"+string(UCI2,"%9.2f")+")"
gen adjusted3=string(HR3, "%9.2f")+" ("+string(lCI3,"%9.2f")+"-"+string(UCI3,"%9.2f")+")"
gen adjusted4=string(HR4, "%9.2f")+" ("+string(lCI4,"%9.2f")+"-"+string(UCI4,"%9.2f")+")"

foreach i in unadjusted adjusted2 adjusted3 adjusted4 {
	replace `i'="ref" if Analysis=="Non-use"
}

keep Database Cancer Mortality Subgroup Analysis Outcomes PersonYears unadjusted adjusted2 adjusted3 adjusted4
order Database Cancer Mortality Subgroup Analysis Outcomes PersonYears unadjusted adjusted2 adjusted3 adjusted4

** keep track of covariates
gen model1="`model1'"
gen model2="`model2'"
gen model3="`model3'"
gen model4="`model4'"

save "D:\CPRD data\Bisphos results\Both\cancerspecific_subgroups_oc.dta", replace


********************************************************************************
** combine into one
********************************************************************************

clear 
use "D:\CPRD data\Bisphos results\Both\cancerspecific_subgroups_ec.dta"
append using  "D:\CPRD data\Bisphos results\Both\cancerspecific_subgroups_oc.dta"

drop if Analysis=="Non-use"

compress 
save  "D:\CPRD data\Bisphos results\Both\cancerspecific_subgroups.dta",replace

export excel using "D:\CPRD data\Bisphos results\Both\CPRD_bisphos.xlsx", sheet("A&G.Sub-groups Cancer-specific",replace) firstrow(var)

 