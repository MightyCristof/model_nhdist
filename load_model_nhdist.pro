PRO load_model_nhdist, subdir


;; assume current directory unless specified
if (n_elements(subdir) eq 0) then path = './' else $
                                  path = subdir+'/'

load_vars,'carroll21_data.sav','_data'
load_vars,'nhdist_obs.sav','_nhobs'
load_vars,'data_prep/rxz_scat01.sav','_rxnh'
pushd,path
load_vars,'ctf_estimate.sav','_ctfest'
load_vars,'split_estimate.sav','_split'
load_vars,'rx_model.sav','_rxmod'
load_vars,'fx_estimate.sav','_fxest'


END




