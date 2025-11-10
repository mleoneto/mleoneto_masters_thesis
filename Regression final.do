*log using stata_output.txt, text replace

clear all
 
cd /Users/leonorneto/Documents/Internship_lab/Testable/Results/Analysis/Stata


insheet using results_tab.csv, clear

*use "C:\Users\roolio\Dropbox\Research\Moral contagion\Results TESTABLE\DATA_Stata.dta", clear

describe
summarize

*tab lie 

replace pred_dec = "." if pred_dec == "NaN"
replace pred_true = "." if pred_true == "NaN"
destring pred_dec pred_true, replace



xtset subj_num trial_tot

gen block_type = 1 if block == 1
replace block_type = 2 if bad_group == 1 & block > 1
replace block_type = 3 if bad_group == 0 & block > 1
*gen bad_group = (block_type == 2)
*replace bad_group = . if block == 1
gen diff_pay = solofalsev - solotruev
gen abs_diff_pred = (pred_true - pred_dec)
gen pred_error = abs(pred_true - pred_dec)
replace pred_error = 0 if pred_error == . 
by subj_num block_type, s: egen prop_lie_block = mean(lie)
by subj_num, s: egen prop_lie_all = mean(lie)




* lag variable; but in fact this pertains to trial t-2 because of the structure of the data set. So I do not think this is really useful. 
gen prev_liars = pred_true[_n-1]



**** COMPARISONS BETWEEN GROUPS ****
preserve
collapse (mean) lie, by(subj_num block_type)
ranksum lie, by(b2.block_type b3.block_type)
restore


preserve
keep if block > 1
collapse (mean) lie, by(subj_num bad_group_first)
ranksum lie , by(bad_group_first)
restore


preserve 
collapse (mean) lie, by(subj_num lockd)
ranksum lie, by(lockd)
restore

xtlogit lie new_cases total_cases new_deaths total_deaths
margins, dydx(*)



**** REGRESSION ****


******** Without predictions *********
xtlogit lie female age seq q_score trial_tot i.bad_group_first b1.block_type c.solotruev##c.solofalsev b1.block_type##trial_tot, or vce(cluster subj_num)
margins, dydx(*)
estat ic
* BIC 5949.636

margins, over(solotruev solofalsev)
marginsplot

******* With predictions *************
preserve
keep if block_type > 1 
xtlogit lie female age seq q_score trial_tot i.bad_group_first b1.block_type c.solotruev##c.solofalsev pred_error, or vce(cluster subj_num)

restore

pwcorr lie pred_error, obs sig

******************** Checking psychology and student!  ***********************
preserve

keep if subj_num > 31
destring student psy, replace
*count if missing(psy)
xtlogit lie female age seq q_score trial_tot i.bad_group_first b1.block_type c.solotruev##c.solofalsev, or vce(cluster subj_num)
margins, dydx(*)
estat ic
* BIC 3649.441

keep if subj_num > 31
destring student psy, replace
xtlogit lie female age seq q_score trial_tot i.psy i.bad_group_first b1.block_type c.solotruev##c.solofalsev, or vce(cluster subj_num)
*margins, dydx(*)
estat ic
* BIC 3748.297 - not really useful to use psy

keep if subj_num > 31
destring student psy, replace
xtlogit lie female age seq q_score trial_tot i.student i.bad_group_first b1.block_type c.solotruev##c.solofalsev, or vce(cluster subj_num)
*margins, dydx(*)
estat ic
* BIC 3747.885 - not useful 

* inlcuding all demographics
keep if subj_num > 31
destring student psy, replace
describe
xtlogit lie female age seq q_score trial_tot i.student i.psy i.bad_group_first b1.block_type c.solotruev##c.solofalsev, or vce(cluster subj_num)
margins, dydx(*)
estat ic
* BIC 3756.756, none of them significant
restore



***** Regression on prediction *****
preserve
keep if block_type > 1
destring predtruev predfalsev, replace
xtreg pred_error female age seq q_score trial_tot i.bad_group_first b2.block_type c.predfalsev c.predtruev i.bad_group_first##b2.block_type part 
margins, dydx(bad_group_first) over(block_type)
marginsplot

describe
summarize
margins, over(solotruev solofalsev)
marginsplot

margins, dxdy(bad_group_first) over(blocktype)
marginsplot
restore


******* Regression on RT **********

