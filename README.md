
### folders

wtreg folder contains '*_wtreg_*.RData' files that are results for each case study and window combination, includes original wq data and detided DO information created in 'source_me.r' of 'detiding_cases.RProj'

pdrnrm folder contains 'prdnrm_*.RData' that are results for simulations, created using 'source_me.r' in 'detiding_sim.RProj'

### data

'met_ls.RData' is list of metabolism estimates for each case/window combination before/after detiding, these were created using 'source_me.r' in 'detiding_cases.RProj'

'case_grds.RData' is data frame of window width combs used to eval case studies, note that only one combination is used in the manuscript, created using 'source_me.r' in 'detiding_cases.RProj'

'comb_grd.RData' is data frame of time series characteristics and window width combinations used to evaluate simulations, created in 'source_me.r' in 'detiding_sim.RProj'

'mod_perf.RData' is summary of simulations, created using 'source_me.r' in 'detiding_sim.RProj'

'cor_res.RData' is correlations of case study data with tide before/after detiding, created in 'cor_res' chunk 'swmp_noise.Rnw'

'met_comp.RData' is summary of mean/anom metabolism data before/after detiding of case studies, created in 'case_res' chunk in 'swmp_noise.Rnw'

'case_tab.RData' is data frame of summary informaiton for case studies, created in 'case_att' chunk in 'swmp_noise.Rnw'