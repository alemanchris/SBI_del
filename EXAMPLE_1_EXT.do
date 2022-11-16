// install sbinorm.ado first
// First part runs SBI on life expectancy
clear all
set more off
clear programs

sysdescribe uslifeexp.dta
sysuse  uslifeexp.dta
// Food and Drug Administration (FDA) approved this vaccine in 1971
drop if year<1946
twoway (scatter le_w year)(scatter le_b year)
gen time = year-1945
lowess le_w time, bwidth(0.15) gen(le_ws) nog
lowess le_b time, bwidth(0.15) gen(le_bs) nog
twoway (line le_ws time)(line le_bs time)(scatter le_w time)(scatter le_b time), xline(26) name(f1) title(Before Normalization)
// 1971 translates into 26
// Normalize 
sbinorm le_ws le_bs time, tp(26)

// Manually get the figure
//drop stage
gen stage = 3.96+1.039*time
//drop cnorm
gen cnorms = 0.908*le_ws
gen cnorm = 0.908*le_w
mata tpn = 3.96+1.039*26
mata tpn
label variable cnorm "Life expenctancy (LE), white"
label variable cnorms "LE norm white (smooth)"
label variable le_bs "LE norm, blacks (smooth)"
label variable le_b "Life expenctancy (LE), blacks"
twoway (scatter cnorm stage)(scatter le_b time)(line cnorms stage)(line le_bs time), xline(26) xline(30.97) name(f2) title(After Normalization)






