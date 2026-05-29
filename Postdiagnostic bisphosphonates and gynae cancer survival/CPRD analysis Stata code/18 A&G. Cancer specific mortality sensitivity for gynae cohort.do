// Sensitivity analyses 
	// Cancer-specific mortality 
		// long term use (2, 3 & 4 years)
		// 1 year exposure lag 
		// 2 year exposure lag
		// adjust for prior bisphos use (3 years prior)
		// new users of bisphos (exclude use in 3 years prior)
		// adjust for aspirin, statin,  metformin, corticosteroids
		// adjust for bmi, smoking, alcohol (seperate then together)
		// adjust for grade
		// include women >40 years old
		// include all cancer stages
		// 2 prescriptions
		// prediagnostic bisphos use and survival
		// active comparator: nitrogen bisphos user vs non-nitrogen bisphos user
		// multiple imputation of missing values
		
		
clear
capture postutil clear
macro drop _all
set trace off

local model1 "" 

local model2 "ageatdx"

local model3 "ageatdx yeardx ncras_chemo6mafterdx ncras_surgery6mafterdx ncras_radio6mafterdx cc_mipriordx cc_chfpriordx cc_pvdpriordx cc_dempriordx cc_cpdpriordx cc_rheupriordx cc_pudpriordx cc_mldpriordx cc_sldpriordx cc_diabpriordx cc_diabcomppriordx cc_plegiapriordx cc_renalpriordx cc_cerebrvdpriordx deprivation_quintile"	

local model4 "ageatdx yeardx ncras_chemo6mafterdx ncras_surgery6mafterdx ncras_radio6mafterdx cc_mipriordx cc_chfpriordx cc_pvdpriordx cc_dempriordx cc_cpdpriordx cc_rheupriordx cc_pudpriordx cc_mldpriordx cc_sldpriordx cc_diabpriordx cc_diabcomppriordx cc_plegiapriordx cc_renalpriordx cc_cerebrvdpriordx deprivation_quintile i.stage"	

postfile results str10 Database str50 Cancer str20 Mortality str100 Sensitivity str40 Analysis Outcomes PersonYears HR lCI UCI Pvalue HR2 lCI2 UCI2 Pvalue2 HR3 lCI3 UCI3 Pvalue3 HR4 lCI4 UCI4 Pvalue4 using "D:\CPRD data\Bisphos results\Both\cancerspecific_sensitivity.dta",replace	

********************************************************************************
**# Long term use (2, 3 & 4 years)
********************************************************************************

foreach c in 1 2 {
    foreach d in 24 36 48 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c'

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
* lag date first
foreach j in gp_1_bisphosafterdxdate gp_`d'_bisphosafterdxdate  {
	replace `j'=`j'+182.5 if `j'~=.
	replace `j'=. if `j'>=date_endfup
	}

* then split 
foreach i in gp_1_bisphosafterdxdate gp_`d'_bisphosafterdxdate {
	stsplit newt`i', at(0) after(time=`i')
	
	gen ever`i'=1 if _t0>=(`i'-_origin)/365.25 & `i'~=.
	replace ever`i'=0 if ever`i'==.
	}

* make variable for duration  
gen bisphostime=0
replace bisphostime=1 if evergp_1_bisphosafterdxdate==1 & evergp_`d'_bisphosafterdxdate==0
replace bisphostime=2 if evergp_1_bisphosafterdxdate==1 & evergp_`d'_bisphosafterdxdate==1

lab val bisphostime bisphostime
lab define bisphostime 0"Non-use" 1"Use <1 months" 2"Use ≥`d' months" 

tab bisphostime,m
tab bisphostime outcome,m

levelsof bisphostime, local(bisphostime)
	foreach i in `bisphostime' {
	
		* database
		local database="CPRD"

		* cancer
		local cancer: label gynaecancer `c'

		* mortality
		local mortality="cancer-specific"
		
		* Sensitivity type 
		local sens "Long term user time varied"
		
		* Analysis 
		local analysis:label bisphostime `i'

		* Outcomes 
		count if outcome==1 & bisphostime==`i'
		local outcomes=r(N)
		
		* person time
		stptime if bisphostime==`i'
		local ptime=r(ptime)
		
		* unadusted HR
		stcox i.bisphostime
		matrix results=r(table)
		local hr=results[1,(1+`i')]
		local lci=results[5,(1+`i')]
		local uci=results[6,(1+`i')]
		local pvalue=results[4,(1+`i')]
		
		* age adjusted
		stcox i.bisphostime `model2'
		matrix results=r(table)
		local 2hr=results[1,(1+`i')]
		local 2lci=results[5,(1+`i')]
		local 2uci=results[6,(1+`i')]
		local 2pvalue=results[4,(1+`i')]	
		
		* 
		stcox i.bisphostime `model3'
		matrix results=r(table)
		local 3hr=results[1,(1+`i')]
		local 3lci=results[5,(1+`i')]
		local 3uci=results[6,(1+`i')]
		local 3pvalue=results[4,(1+`i')]	
		
		*  
		stcox i.bisphostime `model4'
		matrix results=r(table)
		local 4hr=results[1,(1+`i')]
		local 4lci=results[5,(1+`i')]
		local 4uci=results[6,(1+`i')]
		local 4pvalue=results[4,(1+`i')]

		post results ("`database'") ("`cancer'") ("`mortality'") ("`sens'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
	}
}
}

********************************************************************************
**# 1 year exposure lag 
********************************************************************************

foreach c in 1 2 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c'

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
gen lag=365.25
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
	replace gp_m_`i'afterdxdate=gp_m_`i'afterdxdate+lag if gp_m_`i'afterdxdate~=.
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

		* Sensitivity type 
		local sens "Exposure lag: 1 year"
		
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

		post results ("`database'") ("`cancer'") ("`mortality'") ("`sens'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
}


********************************************************************************
**# 2 year exposure lag 
********************************************************************************

foreach c in 1 2 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c'

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
gen lag=730.5
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
	replace gp_m_`i'afterdxdate=gp_m_`i'afterdxdate+lag if gp_m_`i'afterdxdate~=.
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

		* Sensitivity type 
		local sens "Exposure lag: 2 years"
		
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

		post results ("`database'") ("`cancer'") ("`mortality'") ("`sens'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
}

********************************************************************************
**# Adjust for prior bisphos use
********************************************************************************
// in three years prior to diagnosis
// among those with three years medication data 


foreach c in 1 2 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c'

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

* registered with gp less than 3 years
drop if gp_regstartdate>startminus3yrs

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

		* Sensitivity type 
		local sens "Adjust for bisphos use in 3 years prior to diagnosis"
		
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
			
		stcox i.bisphos `model4' gp_m_bisphospriordx3yr
		matrix results=r(table)
		local 4hr=results[1,2]
		local 4lci=results[5,2]
		local 4uci=results[6,2]
		local 4pvalue=results[4,2]	

		post results ("`database'") ("`cancer'") ("`mortality'") ("`sens'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
}



********************************************************************************
**# Exclude prior bisphos users
********************************************************************************
// in three years prior to diagnosis
// among those with three years med data 


foreach c in 1 2 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c'

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

* registered with gp less than 3 years
drop if gp_regstartdate>startminus3yrs

* prior user of bisphosphonates 
drop if gp_m_bisphospriordx3yr==1

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

		* Sensitivity type 
		local sens "New users of bisphos"
		
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

		post results ("`database'") ("`cancer'") ("`mortality'") ("`sens'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
}


********************************************************************************
**# Adjust for statin, aspirin and metformin & corticosteroid use 
********************************************************************************
// in year prior to diagnosis

foreach c in 1 2 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c'

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

		* Sensitivity type 
		local sens "Adjust for aspirin, statin, metformin and oral corticosteroid use prior yr"
		
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
			
		stcox i.bisphos `model4' gp_m_statinpriordxyr gp_m_aspirinpriordxyr gp_m_metforminpriordxyr gp_m_corticosteroidpriordxyr
		matrix results=r(table)
		local 4hr=results[1,2]
		local 4lci=results[5,2]
		local 4uci=results[6,2]
		local 4pvalue=results[4,2]	

		post results ("`database'") ("`cancer'") ("`mortality'") ("`sens'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
}


********************************************************************************
**# Adjust for smoking, alcohol and bmi
********************************************************************************

local sens1 "latestbmipriordx"
local sens2 "i.gp_lf_smokestatuspriordx"
local sens3 "i.gp_lf_drinkerstatuspriordx"
local sens4 "latestbmipriordx i.gp_lf_smokestatuspriordx i.gp_lf_drinkerstatuspriordx"

foreach c in 1 2 {
    
	foreach s in 1 2 3 4 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c'

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

		* Sensitivity type 
		local sens "additional covariate: `sens`s''"
		
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
			
		stcox i.bisphos `model4' `sens`s''
		matrix results=r(table)
		local 4hr=results[1,2]
		local 4lci=results[5,2]
		local 4uci=results[6,2]
		local 4pvalue=results[4,2]	

		post results ("`database'") ("`cancer'") ("`mortality'") ("`sens'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
	}
}



********************************************************************************
**# Adjust for grade
********************************************************************************
// complete case

foreach c in 1 2 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c'

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

		* Sensitivity type 
		local sens "additional covariate: i.grade"
		
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
			
		stcox i.bisphos `model4' i.grade
		matrix results=r(table)
		local 4hr=results[1,2]
		local 4lci=results[5,2]
		local 4uci=results[6,2]
		local 4pvalue=results[4,2]	

		post results ("`database'") ("`cancer'") ("`mortality'") ("`sens'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
	}


********************************************************************************
**# Include women >=40 years old
********************************************************************************

foreach c in 1 2 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c'

* exclude duplicate GP practices
local characteristics="-Excluded: duplicate GP practices"
drop if inlist(pracid, 20024,20036, 20091,20171,20178,20202,20254,20389, 20430,20452, 20469,20487,20552,20554,20640,20717, 20734,20737,20740,20790,20803, 20822,20868,20912,20996, 21001,21015,21078,21112,21118,21172 ,21173, 21277,21281,21331,21334,21390,21430 ,21444,21451,21529,21553,21558,21585 ) & db==1

** cancer diagnosis between 1st April 2003 (start of outpatient) and 31st december 2018 (end of cancer reg data)
drop if ncras_diagnosisdatebest<mdy(4,1,2003) | ncras_diagnosisdatebest>mdy(12,31,2018)

* gold lesss than 1 year uts prior
gen startminusyear=start-365
drop if gp_uts>startminusyear & gp_uts~=.

// exclusion criteria

* drop <40 
drop if ageatdx<40

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

		* Sensitivity type 
		local sens "Women aged 40+ years old"
		
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

		post results ("`database'") ("`cancer'") ("`mortality'") ("`sens'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
}


********************************************************************************
**# Include all cancer stages
********************************************************************************

foreach c in 1 2 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c'

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
*drop if stage==4

* drop missing stage
*drop if stage==.

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

		* Sensitivity type 
		local sens "Include all cancer stages"
		
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

		post results ("`database'") ("`cancer'") ("`mortality'") ("`sens'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
}

********************************************************************************
**# Cancer anywhere on death certificate (cancer-specific)
********************************************************************************

foreach c in 1 2 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c'

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
replace date_outcome=ons_deathdate if ons_ac_ec==1 & gynaecancer==1
replace date_outcome=ons_deathdate if ons_ac_oc==1 & gynaecancer==2
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

		* Sensitivity type 
		local sens "Cancer listed anywhere on death cert"
		
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

		post results ("`database'") ("`cancer'") ("`mortality'") ("`sens'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
}

********************************************************************************
**# 2 prescriptions to become user
********************************************************************************

foreach c in 1 2 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c'

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
foreach i in 2_bisphos  {
	replace gp_`i'afterdxdate=gp_`i'afterdxdate+182.5 if gp_`i'afterdxdate~=.
	replace gp_`i'afterdxdate=. if gp_`i'afterdxdate>=date_endfup
}

* then split 
foreach i in 2_bisphos  {
	stsplit newtime`i', at(0) after(time=gp_`i'afterdxdate)
	
	gen ever`i'=1 if _t0>=(gp_`i'afterdxdate-_origin)/365.25 & gp_`i'afterdxdate~=.
	replace ever`i'=0 if ever`i'==.
}

* make variable for each time period
gen bisphos=1 if ever2_bisphos==0
replace bisphos=2 if ever2_bisphos==1
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

		* Sensitivity type 
		local sens "Two prescriptions for user"
		
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

		post results ("`database'") ("`cancer'") ("`mortality'") ("`sens'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
}


********************************************************************************
**# Active comparator nitrogen bisphos vs non-nitrogen bisphos
********************************************************************************

foreach c in 1 2 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c'

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

tab bisphos,m
tab bisphos outcome,m

* nitrogen containing or not
// note: first need to lag to match up dates
foreach i in alendronic risedronate ibandronic etidronate tiludronic  {
	replace gp_m_`i'afterdxdate=gp_m_`i'afterdxdate+182.5 if gp_m_`i'afterdxdate~=.
	replace gp_m_`i'afterdxdate=. if gp_m_`i'afterdxdate>=date_endfup
	}

gen nbisphos=1 if bisphos==1
replace nbisphos=2 if bisphos==2 & (gp_m_bisphosafterdxdate==gp_m_alendronicafterdxdate | gp_m_bisphosafterdxdate==gp_m_risedronateafterdxdate | gp_m_bisphosafterdxdate==gp_m_ibandronicafterdxdate)
replace nbisphos=3 if bisphos==2 & (gp_m_bisphosafterdxdate==gp_m_etidronateafterdxdate | gp_m_bisphosafterdxdate==gp_m_tiludronicafterdxdate)
lab val nbisphos nbisphos 
lab define nbisphos 1"Non-use" 2"Nitrogen-containing" 3"Non-nitrogen containing"

tab nbisphos,m 
tab nbisphos outcome,m 

levelsof nbisphos, local(nbisphos)
	foreach i in `nbisphos' {
	
		* database
		local database="CPRD"

		* cancer
		local cancer: label gynaecancer `c'

		* mortality
		local mortality="cancer-specific"
		
		* Sensitivity type 
		local sens "Active comparator: Nitrogen vs non-nitrogen containing"
		
		* Analysis 
		local analysis:label nbisphos `i'

		* Outcomes 
		count if outcome==1 & nbisphos==`i'
		local outcomes=r(N)
		
		* person time
		stptime if nbisphos==`i'
		local ptime=r(ptime)

		* unadusted HR
		stcox ib3.nbisphos
		matrix results=r(table)
		local hr=results[1,(0+`i')]
		local lci=results[5,(0+`i')]
		local uci=results[6,(0+`i')]
		local pvalue=results[4,(0+`i')]

		* age adjusted
		stcox ib3.nbisphos `model2'
		matrix results=r(table)
		local 2hr=results[1,(0+`i')]
		local 2lci=results[5,(0+`i')]
		local 2uci=results[6,(0+`i')]
		local 2pvalue=results[4,(0+`i')]	
			
		stcox ib3.nbisphos `model3'
		matrix results=r(table)
		local 3hr=results[1,(0+`i')]
		local 3lci=results[5,(0+`i')]
		local 3uci=results[6,(0+`i')]
		local 3pvalue=results[4,(0+`i')]	
			
		stcox ib3.nbisphos `model4'
		matrix results=r(table)
		local 4hr=results[1,(0+`i')]
		local 4lci=results[5,(0+`i')]
		local 4uci=results[6,(0+`i')]
		local 4pvalue=results[4,(0+`i')]	

		post results ("`database'") ("`cancer'") ("`mortality'") ("`sens'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
}

postclose results	

clear
use "D:\CPRD data\Bisphos results\Both\cancerspecific_sensitivity.dta"

** formatting
gen unadjusted=string(HR, "%9.2f")+" ("+string(lCI,"%9.2f")+"-"+string(UCI,"%9.2f")+")"
gen adjusted2=string(HR2, "%9.2f")+" ("+string(lCI2,"%9.2f")+"-"+string(UCI2,"%9.2f")+")"
gen adjusted3=string(HR3, "%9.2f")+" ("+string(lCI3,"%9.2f")+"-"+string(UCI3,"%9.2f")+")"
gen adjusted4=string(HR4, "%9.2f")+" ("+string(lCI4,"%9.2f")+"-"+string(UCI4,"%9.2f")+")"

foreach i in unadjusted adjusted2 adjusted3 adjusted4 {
	replace `i'="ref" if Analysis=="Non-use"
}

keep Database Cancer Mortality Sensitivity Analysis Outcomes PersonYears unadjusted adjusted2 adjusted3 adjusted4
order Database Cancer Mortality Sensitivity Analysis Outcomes PersonYears unadjusted adjusted2 adjusted3 adjusted4

** keep track of covariates
gen model1="`model1'"
gen model2="`model2'"
gen model3="`model3'"
gen model4="`model4'"

replace model4="`model4'"+" latestbmipriordx" if Sensitivity=="additional covariate: latestbmipriordx"
replace model4="`model4'"+" i.gp_lf_smokestatuspriordx" if Sensitivity=="additional covariate: i.gp_lf_smokestatuspriordx"
replace model4="`model4'"+" i.gp_lf_drinkerstatuspriordx" if Sensitivity=="additional covariate: i.gp_lf_drinkerstatuspriordx"
replace model4="`model4'"+" latestbmipriordx i.gp_lf_smokestatuspriordx i.gp_lf_drinkerstatuspriordx" if Sensitivity=="additional covariate: latestbmipriordx i.gp_lf_smokestatuspriordx i.gp_lf_drinkerstatuspriordx"

replace model4="`model4'"+" gp_m_bisphospriordx3yr" if Sensitivity=="Adjust for bisphos use in 3 years prior to diagnosis"
replace model4="`model4'"+" gp_m_statinpriordxyr gp_m_aspirinpriordxyr gp_m_metforminpriordxyr " if Sensitivity=="Adjust for aspirin, statin and metformin use prior yr"

drop if strpos(Analysis,"<")>=1

save "D:\CPRD data\Bisphos results\Both\cancerspecific_sensitivity.dta", replace


********************************************************************************
**# Prediagnositic bisphos use and survival
********************************************************************************

clear
capture postutil clear
macro drop _all
set trace off

local model1 "" 

local model2 "ageatdx"

local model3 "ageatdx yeardx cc_mipriordx cc_chfpriordx cc_pvdpriordx cc_dempriordx cc_cpdpriordx cc_rheupriordx cc_pudpriordx cc_mldpriordx cc_sldpriordx cc_diabpriordx cc_diabcomppriordx cc_plegiapriordx cc_renalpriordx cc_cerebrvdpriordx deprivation_quintile"	


postfile results str10 Database str50 Cancer str20 Mortality str100 Sensitivity str40 Analysis Outcomes PersonYears HR lCI UCI Pvalue HR2 lCI2 UCI2 Pvalue2 HR3 lCI3 UCI3 Pvalue3 HR4 lCI4 UCI4 Pvalue4 using "D:\CPRD data\Bisphos results\Both\cancerspecific_sensitivity_before.dta",replace	


foreach c in 1 2 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c'

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
gen lag=0
gen date_entry=start+lag
gen date_endfup=min(ons_deathdate,mdy(11,16,2020), gp_lastcollectiondate, gp_transferoutdate, gp_endrecords, gp_regenddate)

* less than 6 months of follow-up
*drop if date_entry>date_endfup

* registered with gp less than 3 years
drop if gp_regstartdate>startminus3yrs

gen date_outcome=.
replace date_outcome=ons_deathdate if ons_pc_ec==1 & gynaecancer==1
replace date_outcome=ons_deathdate if ons_pc_oc==1 & gynaecancer==2
gen outcome=1 if date_outcome~=.

* check if outcome occurs after end of fup
replace outcome=. if date_outcome>date_endfup

count 

// stset data 
stset date_endfup, fail(outcome=1) enter(time date_entry) id(patid) origin(time date_entry) scale(365.25) 

* make variable
lab val gp_bisphospriordx3yr gp_bisphospriordx3yr 
lab define bisphos 0"Non-use" 1"Any use in 3 years prior"

levelsof gp_bisphospriordx3yr, local(bisphos)
	foreach i in `bisphos' {
	
		* database
		local database="CPRD"

		* cancer
		local cancer: label gynaecancer `c'

		* mortality
		local mortality="cancer-specific"

		* Sensitivity type 
		local sens "Prediagnostic bisphos use in 3 years prior to dx (no treatment adj)"
		
		* Analysis 
		local analysis:label bisphos `i'

		* Outcomes 
		count if outcome==1 & gp_bisphospriordx3yr==`i'
		local outcomes=r(N)
		
		* person time
		stptime if gp_bisphospriordx3yr==`i'
		local ptime=r(ptime)

		* unadusted HR
		stcox i.gp_bisphospriordx3yr
		matrix results=r(table)
		local hr=results[1,(2)]
		local lci=results[5,2]
		local uci=results[6,2]
		local pvalue=results[4,2]

		* age adjusted
		stcox i.gp_bisphospriordx3yr `model2'
		matrix results=r(table)
		local 2hr=results[1,2]
		local 2lci=results[5,2]
		local 2uci=results[6,2]
		local 2pvalue=results[4,2]	
			
		stcox i.gp_bisphospriordx3yr `model3'
		matrix results=r(table)
		local 3hr=results[1,2]
		local 3lci=results[5,2]
		local 3uci=results[6,2]
		local 3pvalue=results[4,2]	
			
		stcox i.gp_bisphospriordx3yr `model3'
		matrix results=r(table)
		local 4hr=results[1,2]
		local 4lci=results[5,2]
		local 4uci=results[6,2]
		local 4pvalue=results[4,2]	

		post results ("`database'") ("`cancer'") ("`mortality'") ("`sens'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
}


foreach c in 1 2 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c'

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
gen lag=0
gen date_entry=start+lag
gen date_endfup=min(ons_deathdate,mdy(11,16,2020), gp_lastcollectiondate, gp_transferoutdate, gp_endrecords, gp_regenddate)

* less than 6 months of follow-up
*drop if date_entry>date_endfup

* registered with gp less than 3 years
drop if gp_regstartdate>startminus3yrs

gen date_outcome=.
replace date_outcome=ons_deathdate if ons_pc_ec==1 & gynaecancer==1
replace date_outcome=ons_deathdate if ons_pc_oc==1 & gynaecancer==2
gen outcome=1 if date_outcome~=.

* check if outcome occurs after end of fup
replace outcome=. if date_outcome>date_endfup

count 

// stset data 
stset date_endfup, fail(outcome=1) enter(time date_entry) id(patid) origin(time date_entry) scale(365.25) 

* make variable
lab val gp_bisphospriordx3yr gp_bisphospriordx3yr 
lab define bisphos 0"Non-use" 1"Any use in 3 years prior"

*time-varying treatments
foreach i in chemo surgery radio  {
	stsplit newtime`i', at(0) after(time=ncras_`i'6mafterdxdate)
	
	gen tv`i'=1 if _t0>=(ncras_`i'6mafterdxdate-_origin)/365.25 & ncras_`i'6mafterdxdate~=.
	replace tv`i'=0 if tv`i'==.
}

levelsof gp_bisphospriordx3yr, local(bisphos)
	foreach i in `bisphos' {
	
		* database
		local database="CPRD"

		* cancer
		local cancer: label gynaecancer `c'

		* mortality
		local mortality="cancer-specific"

		* Sensitivity type 
		local sens "Prediagnostic bisphos use in 3 years prior to dx (stage & tv treatment adj)"
		
		* Analysis 
		local analysis:label bisphos `i'

		* Outcomes 
		count if outcome==1 & gp_bisphospriordx3yr==`i'
		local outcomes=r(N)
		
		* person time
		stptime if gp_bisphospriordx3yr==`i'
		local ptime=r(ptime)

		* unadusted HR
		stcox i.gp_bisphospriordx3yr
		matrix results=r(table)
		local hr=results[1,(2)]
		local lci=results[5,2]
		local uci=results[6,2]
		local pvalue=results[4,2]

		* age adjusted
		stcox i.gp_bisphospriordx3yr `model2'
		matrix results=r(table)
		local 2hr=results[1,2]
		local 2lci=results[5,2]
		local 2uci=results[6,2]
		local 2pvalue=results[4,2]	
			
		stcox i.gp_bisphospriordx3yr `model3' tvsurgery tvchemo tvradio
		matrix results=r(table)
		local 3hr=results[1,2]
		local 3lci=results[5,2]
		local 3uci=results[6,2]
		local 3pvalue=results[4,2]	
			
		stcox i.gp_bisphospriordx3yr `model3' tvsurgery tvchemo tvradio i.stage
		matrix results=r(table)
		local 4hr=results[1,2]
		local 4lci=results[5,2]
		local 4uci=results[6,2]
		local 4pvalue=results[4,2]	

		post results ("`database'") ("`cancer'") ("`mortality'") ("`sens'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
}


postclose results	

clear
use "D:\CPRD data\Bisphos results\Both\cancerspecific_sensitivity_before.dta"

** formatting
gen unadjusted=string(HR, "%9.2f")+" ("+string(lCI,"%9.2f")+"-"+string(UCI,"%9.2f")+")"
gen adjusted2=string(HR2, "%9.2f")+" ("+string(lCI2,"%9.2f")+"-"+string(UCI2,"%9.2f")+")"
gen adjusted3=string(HR3, "%9.2f")+" ("+string(lCI3,"%9.2f")+"-"+string(UCI3,"%9.2f")+")"
gen adjusted4=string(HR4, "%9.2f")+" ("+string(lCI4,"%9.2f")+"-"+string(UCI4,"%9.2f")+")"

foreach i in unadjusted adjusted2 adjusted3 adjusted4  {
	replace `i'="ref" if Analysis=="Non-use"
}

keep Database Cancer Mortality Sensitivity Analysis Outcomes PersonYears unadjusted adjusted2 adjusted3 adjusted4
order Database Cancer Mortality Sensitivity Analysis Outcomes PersonYears unadjusted adjusted2 adjusted3 adjusted4 

** keep track of covariates
gen model1="`model1'"
gen model2="`model2'"
gen model3="`model3'"
gen model4=""
replace model4="`model3' i.stage tvsurgery tvchemo tvradio" if strpos(Sensitivity,"tv")>=1  


save "D:\CPRD data\Bisphos results\Both\cancerspecific_sensitivity_before.dta", replace


********************************************************************************
**# MI for missing BMI, smoking, alcohol status, deprivation and grade
********************************************************************************
// Note: had to remove severe liver disease to prevent error from perfect prediction. Identified by huge confidence interval in output that failed.

// Endometrial cancer 

clear
capture postutil clear
macro drop _all
set trace off

local model1 "" 

local model2 "ageatdx"

local model3 "ageatdx yeardx ncras_chemo6mafterdx ncras_surgery6mafterdx ncras_radio6mafterdx cc_mipriordx cc_chfpriordx cc_pvdpriordx cc_cpdpriordx cc_rheupriordx cc_pudpriordx cc_mldpriordx cc_diabpriordx cc_dempriordx cc_renalpriordx cc_diabcomppriordx cc_plegiapriordx cc_cerebrvdpriordx deprivation_quintile"	

local model4 "ageatdx yeardx ncras_chemo6mafterdx ncras_surgery6mafterdx ncras_radio6mafterdx cc_mipriordx cc_chfpriordx cc_pvdpriordx cc_cpdpriordx cc_rheupriordx cc_pudpriordx cc_mldpriordx cc_diabpriordx  cc_dempriordx cc_renalpriordx cc_diabcomppriordx cc_plegiapriordx cc_cerebrvdpriordx deprivation_quintile i.stage i.grade"	

postfile results str10 Database str50 Cancer str20 Mortality str100 Sensitivity str40 Analysis Outcomes PersonYears HR lCI UCI Pvalue HR2 lCI2 UCI2 Pvalue2 HR3 lCI3 UCI3 Pvalue3 HR4 lCI4 UCI4 Pvalue4 using "D:\CPRD data\Bisphos results\Both\cancerspecific_sensitivity_mi_ec.dta",replace	


foreach c in 1 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c'

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

gen death=1 if outcome==1 
replace death=0 if outcome==.

count 

// stset data 
stset date_endfup, fail(outcome=1) enter(time date_entry) id(patid) origin(time date_entry) scale(365.25) 

// mi 
codebook latestbmipriordx gp_lf_smokestatuspriordx gp_lf_drinkerstatuspriordx
codebook grade deprivation_quintile

sts gen HT=na
mi set wide
mi register imputed latestbmipriordx gp_lf_smokestatuspriordx gp_lf_drinkerstatuspriordx grade deprivation_quintile
mi impute chained (regress) latestbmipriordx (mlogit) gp_lf_smokestatuspriordx gp_lf_drinkerstatuspriordx (ologit) deprivation_quintile grade  = death HT gp_m_bisphosafterdx ageatdx yeardx ncras_chemo6mafterdx ncras_surgery6mafterdx ncras_radio6mafterdx cc_mipriordx cc_chfpriordx cc_pvdpriordx cc_dempriordx cc_cpdpriordx cc_rheupriordx cc_pudpriordx cc_mldpriordx /*cc_sldpriordx*/ cc_diabpriordx cc_diabcomppriordx cc_plegiapriordx cc_renalpriordx cc_cerebrvdpriordx stage , noisily add(20) rseed(1000) 


** stsplit to create time-varying exposure definition

* lag dates by six months first 
foreach i in bisphos  {
	replace gp_m_`i'afterdxdate=gp_m_`i'afterdxdate+182.5 if gp_m_`i'afterdxdate~=.
	replace gp_m_`i'afterdxdate=. if gp_m_`i'afterdxdate>=date_endfup
}

* then split 
foreach i in bisphos  {
	mi stsplit newtime`i', at(0) after(time=gp_m_`i'afterdxdate)
	
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

		* Sensitivity type 
		local sens "Multiple imputation: BMI, smoking, alcohol, grade and deprivation"
		
		* Analysis 
		local analysis:label bisphos `i'

		* Outcomes 
		count if outcome==1 & bisphos==`i'
		local outcomes=r(N)
		
		* person time
		*stptime if bisphos==`i'
		local ptime=r(ptime)

		* unadusted HR
		mi estimate, hr: stcox i.bisphos
		matrix results=r(table)
		local hr=results[1,(2)]
		local lci=results[5,2]
		local uci=results[6,2]
		local pvalue=results[4,2]

		* age adjusted
		mi estimate, hr: stcox i.bisphos `model2'
		matrix results=r(table)
		local 2hr=results[1,2]
		local 2lci=results[5,2]
		local 2uci=results[6,2]
		local 2pvalue=results[4,2]	
			
		mi estimate, hr: stcox i.bisphos `model3'
		matrix results=r(table)
		local 3hr=results[1,2]
		local 3lci=results[5,2]
		local 3uci=results[6,2]
		local 3pvalue=results[4,2]	
			
		mi estimate, hr: stcox i.bisphos `model4' latestbmipriordx i.gp_lf_smokestatuspriordx i.gp_lf_drinkerstatuspriordx
		matrix results=r(table)
		local 4hr=results[1,2]
		local 4lci=results[5,2]
		local 4uci=results[6,2]
		local 4pvalue=results[4,2]	

		post results ("`database'") ("`cancer'") ("`mortality'") ("`sens'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
}

postclose results	

clear
use "D:\CPRD data\Bisphos results\Both\cancerspecific_sensitivity_mi_ec.dta"

** formatting
gen unadjusted=string(HR, "%9.2f")+" ("+string(lCI,"%9.2f")+"-"+string(UCI,"%9.2f")+")"
gen adjusted2=string(HR2, "%9.2f")+" ("+string(lCI2,"%9.2f")+"-"+string(UCI2,"%9.2f")+")"
gen adjusted3=string(HR3, "%9.2f")+" ("+string(lCI3,"%9.2f")+"-"+string(UCI3,"%9.2f")+")"
gen adjusted4=string(HR4, "%9.2f")+" ("+string(lCI4,"%9.2f")+"-"+string(UCI4,"%9.2f")+")"

foreach i in unadjusted adjusted2 adjusted3 adjusted4 {
	replace `i'="ref" if Analysis=="Non-use"
}

keep Database Cancer Mortality Sensitivity Analysis Outcomes PersonYears unadjusted adjusted2 adjusted3 adjusted4
order Database Cancer Mortality Sensitivity Analysis Outcomes PersonYears unadjusted adjusted2 adjusted3 adjusted4

** keep track of covariates
gen model1="`model1'"
gen model2="`model2'"
gen model3="`model3'"
gen model4="`model4'"

replace model4="`model4'"+" latestbmipriordx gp_lf_smokestatuspriordx gp_lf_drinkerstatuspriordx" 

save "D:\CPRD data\Bisphos results\Both\cancerspecific_sensitivity_mi_ec.dta", replace
 
 
** Ovarian MI 

// note: had to remove severe liver disease, paraplegia and dementia to prevent prefect prediction.

clear
capture postutil clear
macro drop _all
set trace off

local model1 "" 

local model2 "ageatdx"

local model3 "ageatdx yeardx ncras_chemo6mafterdx ncras_surgery6mafterdx ncras_radio6mafterdx cc_mipriordx cc_chfpriordx cc_pvdpriordx cc_cpdpriordx cc_rheupriordx cc_pudpriordx cc_mldpriordx cc_diabpriordx cc_renalpriordx cc_diabcomppriordx cc_sldpriordx cc_cerebrvdpriordx deprivation_quintile"	

local model4 "ageatdx yeardx ncras_chemo6mafterdx ncras_surgery6mafterdx ncras_radio6mafterdx cc_mipriordx cc_chfpriordx cc_pvdpriordx cc_cpdpriordx cc_rheupriordx cc_pudpriordx cc_mldpriordx cc_diabpriordx cc_renalpriordx cc_diabcomppriordx cc_cerebrvdpriordx deprivation_quintile i.stage i.grade"	

postfile results str10 Database str50 Cancer str20 Mortality str100 Sensitivity str40 Analysis Outcomes PersonYears HR lCI UCI Pvalue HR2 lCI2 UCI2 Pvalue2 HR3 lCI3 UCI3 Pvalue3 HR4 lCI4 UCI4 Pvalue4 using "D:\CPRD data\Bisphos results\Both\cancerspecific_sensitivity_mi_oc.dta",replace	

foreach c in 2 {

use "D:\CPRD data\Datasets\Both\gynaecohort_cleaned.dta", clear

keep if gynaecancer==`c'

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

gen death=1 if outcome==1 
replace death=0 if outcome==.

count 

// stset data 
stset date_endfup, fail(outcome=1) enter(time date_entry) id(patid) origin(time date_entry) scale(365.25) 

// mi 
codebook latestbmipriordx gp_lf_smokestatuspriordx gp_lf_drinkerstatuspriordx
codebook grade deprivation_quintile

sts gen HT=na
mi set wide
mi register imputed latestbmipriordx gp_lf_smokestatuspriordx gp_lf_drinkerstatuspriordx grade deprivation_quintile
mi impute chained (regress) latestbmipriordx (mlogit) gp_lf_smokestatuspriordx gp_lf_drinkerstatuspriordx (ologit) deprivation_quintile grade  = death HT gp_m_bisphosafterdx ageatdx yeardx ncras_chemo6mafterdx ncras_surgery6mafterdx ncras_radio6mafterdx cc_mipriordx cc_chfpriordx cc_pvdpriordx /*cc_dempriordx*/ cc_cpdpriordx cc_rheupriordx cc_pudpriordx cc_mldpriordx /*cc_sldpriordx*/ cc_diabpriordx cc_diabcomppriordx /*cc_plegiapriordx*/ cc_renalpriordx cc_cerebrvdpriordx stage , noisily add(20) rseed(1000) 


** stsplit to create time-varying exposure definition

* lag dates by six months first 
foreach i in bisphos  {
	replace gp_m_`i'afterdxdate=gp_m_`i'afterdxdate+182.5 if gp_m_`i'afterdxdate~=.
	replace gp_m_`i'afterdxdate=. if gp_m_`i'afterdxdate>=date_endfup
}

* then split 
foreach i in bisphos  {
	mi stsplit newtime`i', at(0) after(time=gp_m_`i'afterdxdate)
	
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

		* Sensitivity type 
		local sens "Multiple imputation: BMI, smoking, alcohol, grade and deprivation"
		
		* Analysis 
		local analysis:label bisphos `i'

		* Outcomes 
		count if outcome==1 & bisphos==`i'
		local outcomes=r(N)
		
		* person time
		*stptime if bisphos==`i'
		local ptime=r(ptime)

		* unadusted HR
		mi estimate, hr: stcox i.bisphos
		matrix results=r(table)
		local hr=results[1,(2)]
		local lci=results[5,2]
		local uci=results[6,2]
		local pvalue=results[4,2]

		* age adjusted
		mi estimate, hr: stcox i.bisphos `model2'
		matrix results=r(table)
		local 2hr=results[1,2]
		local 2lci=results[5,2]
		local 2uci=results[6,2]
		local 2pvalue=results[4,2]	
			
		mi estimate, hr: stcox i.bisphos `model3'
		matrix results=r(table)
		local 3hr=results[1,2]
		local 3lci=results[5,2]
		local 3uci=results[6,2]
		local 3pvalue=results[4,2]	
			
		mi estimate, hr: stcox i.bisphos `model4' latestbmipriordx i.gp_lf_smokestatuspriordx i.gp_lf_drinkerstatuspriordx
		matrix results=r(table)
		local 4hr=results[1,2]
		local 4lci=results[5,2]
		local 4uci=results[6,2]
		local 4pvalue=results[4,2]	

		post results ("`database'") ("`cancer'") ("`mortality'") ("`sens'") ("`analysis'") (`outcomes') (`ptime') (`hr') (`lci') (`uci') (`pvalue') (`2hr') (`2lci') (`2uci') (`2pvalue') (`3hr') (`3lci') (`3uci') (`3pvalue') (`4hr') (`4lci') (`4uci') (`4pvalue')
		}
}

postclose results	

clear
use "D:\CPRD data\Bisphos results\Both\cancerspecific_sensitivity_mi_oc.dta"

** formatting
gen unadjusted=string(HR, "%9.2f")+" ("+string(lCI,"%9.2f")+"-"+string(UCI,"%9.2f")+")"
gen adjusted2=string(HR2, "%9.2f")+" ("+string(lCI2,"%9.2f")+"-"+string(UCI2,"%9.2f")+")"
gen adjusted3=string(HR3, "%9.2f")+" ("+string(lCI3,"%9.2f")+"-"+string(UCI3,"%9.2f")+")"
gen adjusted4=string(HR4, "%9.2f")+" ("+string(lCI4,"%9.2f")+"-"+string(UCI4,"%9.2f")+")"

foreach i in unadjusted adjusted2 adjusted3 adjusted4 {
	replace `i'="ref" if Analysis=="Non-use"
}

keep Database Cancer Mortality Sensitivity Analysis Outcomes PersonYears unadjusted adjusted2 adjusted3 adjusted4
order Database Cancer Mortality Sensitivity Analysis Outcomes PersonYears unadjusted adjusted2 adjusted3 adjusted4

** keep track of covariates
gen model1="`model1'"
gen model2="`model2'"
gen model3="`model3'"
gen model4="`model4'"

replace model4="`model4'"+" latestbmipriordx gp_lf_smokestatuspriordx gp_lf_drinkerstatuspriordx" 

save "D:\CPRD data\Bisphos results\Both\cancerspecific_sensitivity_mi_oc.dta", replace
  

********************************************************************************
** combine 
********************************************************************************

clear
use "D:\CPRD data\Bisphos results\Both\cancerspecific_sensitivity.dta"

foreach i in before mi_oc mi_ec {
    append using "D:\CPRD data\Bisphos results\Both\cancerspecific_sensitivity_`i'.dta"
}

sort Cancer Sensitivity
drop if Analysis=="Non-use"

drop if strpos(Analysis,"Use <1 year")>=1 | strpos(Analysis,"Use ≥1 year")>=1
drop if strpos(Analysis,"Non-nitrogen containing")>=1

order Database Cancer Mortality Sensitivity Analysis 

* run analysis key do file (needed later for pooling)
do "Sensitivity analysis labels.do"

compress
save "D:\CPRD data\Bisphos results\Both\cancerspecific_sensitivity_all.dta",replace

export excel using "D:\CPRD data\Bisphos results\Both\CPRD_bisphos.xlsx", sheet("A&G.Sensitivity Cancer-specific",replace) firstrow(var)

