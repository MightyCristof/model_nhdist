PRO model_nhdist, subdir, $
                  DATA = data, $
                  NHOBS = nhobs, $
                  CTFEST = ctfest, $
                  SPLIT = split, $
                  RXMOD = RXMOD, $
                  FXEST = fxest


;; check for keyword commands
nkeys = n_elements(data) + $
        n_elements(nhobs) + $
        n_elements(ctfest) + $
        n_elements(split) + $
        n_elements(rxmod) + $
        n_elements(fxest)
if (nkeys eq 0) then GOTO, NO_KEYS

;; assume current directory unless specified
if (n_elements(subdir) eq 0) then path = './' else $
                                  path = subdir+'/' & file_mkdir,path

;; pass data from XRAY_LACK_AGN to MODEL_NHDIST
if keyword_set(data) then begin
    pull_carroll21_data,mode='sec'
    nkeys--
endif
load_vars,'carroll21_data.sav','_data'

;; pull NH distributions from Lansbury+2015 and Ananna+2019
if keyword_set(nhobs) then begin
    prep_nhdist_obs
    nkeys--
endif
load_vars,'nhdist_obs.sav','_nhobs'
if (nkeys eq 0) then GOTO, NO_KEYS

;; load RX-NH conversions
load_vars,'data_prep/rxz_scat01.sav','_rxnh'

;; directory for output
pushd,path

;; estimate the distribution of CTF for modeling
if keyword_set(ctfest) then begin
    estimate_ctf
    nkeys--
endif
load_vars,'ctf_estimate.sav','_ctfest'
if (nkeys eq 0) then GOTO, NO_KEYS

;; estimate the distribution of NH=24-25 split for modeling
if keyword_set(split) then begin
    estimate_split
    nkeys--
endif
load_vars,'split_estimate.sav','_split'
if (nkeys eq 0) then GOTO, NO_KEYS

;; simulate the NH and RL distributions
if keyword_set(rxmod) then begin
    model_rxdist
    nkeys--
endif
load_vars,'rx_model.sav','_rxmod'
if (nkeys eq 0) then GOTO, NO_KEYS

;; estimate X-ray flux from model RL distribution
if keyword_set(fxest) then begin
    estimate_fx,/cha
    nkeys--
endif
load_vars,'fx_estimate.sav','_fxest'
if (nkeys eq 0) then GOTO, NO_KEYS



NO_KEYS:
END



;; lir vs z - like in chen17
;; lx vs lir - like in chen17
;; histogram LL
;; histogram NH + modeled
