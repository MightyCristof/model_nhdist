PRO load_model_nhdist


load_vars,'carroll20_data.sav','_data'
load_vars,'nhdist_obs.sav','_nhobs'
load_vars,'ctf_estimate.sav','_ctfest'
load_vars,'split_estimate.sav','_split'
load_vars,'rx_model.sav','_rxmod'
load_vars,'fx_estimate.sav','_fxest'


END

;edf, rand, xrand, yrand
;cdf = call_function('cdf', xrand, _extra = extra)
;delta = edf - cdf
