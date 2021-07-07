PRO load_model_nhdist, subdir


;; assume current directory unless specified
if (n_elements(subdir) eq 0) then path = './' else $
                                  path = subdir+'/'

;load_vars,'carroll21b_data.sav','_data'
load_vars,'carroll21_data.sav','_data'
load_vars,'observed_nhdist.sav','_nhdist'
pushd,path
load_vars,'select_nhobs.sav','_nhobs'
load_vars,'rx_conversion.sav','_rxnh'
load_vars,'select_group.sav','_group'
load_vars,'ctf_fixed.sav','_fixed'
load_vars,'ctf_free.sav','_free'
load_vars,'nh_split.sav','_split'
load_vars,'rx_model.sav','_model'


END




