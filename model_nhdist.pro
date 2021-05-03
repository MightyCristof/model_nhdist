PRO model_nhdist, sub_dir, $
                  DATA = data, $
                  NHDIST = nhdist, $
                  SETNH = setnh, $
                  SETRX = setrx, $
                  GROUP = group, $
                  FIXED = fixed, $
                  FREE = free, $
                  SPLIT = split, $
                  MODEL = model, $
                  TEST = test, $
                  POSTMOD = postmod

;; check for keyword commands
nkeys = n_elements(data) + $
        n_elements(nhdist) + $
        n_elements(setnh) + $
        n_elements(setrx) + $
        n_elements(group) + $
        n_elements(fixed) + $
        n_elements(free) + $
        n_elements(split) + $
        n_elements(model)
if (nkeys eq 0) then GOTO, NO_KEYS

;; assume current directory unless specified
if (n_elements(sub_dir) eq 0) then sub_dir = './' else $
                                   sub_dir = strlowcase(sub_dir)+'/' & file_mkdir,sub_dir

;; pass data from XRAY_LACK_AGN to MODEL_NHDIST
if keyword_set(data) then begin
    pull_carroll21_data
    nkeys--
endif
load_vars,'carroll21_data.sav','_data'

;; pull NH distributions from Lansbury+2015 and Ananna+2019
if keyword_set(nhdist) then begin
    prep_obs_nhdist
    nkeys--
endif
load_vars,'observed_nhdist.sav','_nhdist'
if (nkeys eq 0) then GOTO, NO_KEYS

;; directory for output
pushd,sub_dir

;; select data group WISE AGN, secondary sources, or all (WAC/SEC/ALL)
if keyword_set(setnh) then begin
    if (typename(setnh) ne 'STRING') then message, 'PLEASE SELECT OBSERVED NH: NH_LAN_NST/NH_LAN_CHA/NH_LAN_COR/NH_RIC_AVG'
    set_nhobs,nhobs=setnh
    nkeys--
endif
load_vars,'select_nhobs.sav','_nhobs'
if (nkeys eq 0) then GOTO, NO_KEYS

;; load RX-NH conversions
if keyword_set(setrx) then begin
    prep_dir = file_search('data_prep')
    up = ''
    while (prep_dir eq '') do begin
        up = up+'../'
        prep_dir = file_search(up+'data_prep')
    endwhile
    ;rx_path = prep_dir+'/rxz_scat01_nh2025.sav'
    ;rx_path = prep_dir+'/rxz_scat01_nh20.525.sav'
    ;rx_path = prep_dir+'/test_param_rx.sav'
    rx_path = prep_dir+'/rxz_scat01.sav'
    ;rx_path = prep_dir+'/rxz_scat01_power.sav'
    ;; -f force to overwrite symbolic link in case changing
    spawn, 'ln -sf '+rx_path+' rx_conversion.sav'
    nkeys--
endif
load_vars,'rx_conversion.sav','_rxnh'
if (nkeys eq 0) then GOTO, NO_KEYS

;; select data group WISE AGN, secondary sources, or all (WAC/SEC/ALL)
if keyword_set(group) then begin
    if (typename(group) ne 'STRING') then message, 'PLEASE SELECT INPUT GROUP: WAC/SEC/ALL/WAC_HIX/WAC_LOX/OFFSET'
    set_data_group,group=group,test=test
    nkeys--
endif
load_vars,'select_group.sav','_group'
if (nkeys eq 0) then GOTO, NO_KEYS

;; estimate the distribution of fixed CTF
if keyword_set(fixed) then begin
    estimate_fixed_ctf;_update
    nkeys--
endif
load_vars,'ctf_fixed.sav','_fixed'
if (nkeys eq 0) then GOTO, NO_KEYS

;; estimate the distribution of free CTF
if keyword_set(free) then begin
    estimate_free_ctf;_update
    nkeys--
endif
load_vars,'ctf_free.sav','_free'
if (nkeys eq 0) then GOTO, NO_KEYS

;; estimate the distribution of NH=24-25 split
if keyword_set(split) then begin
    estimate_split_nh
    nkeys--
endif
load_vars,'nh_split.sav','_split'
if (nkeys eq 0) then GOTO, NO_KEYS

;; simulate the NH and RL distributions
if keyword_set(model) then begin
    model_rxdist,postmod=postmod;_update
    nkeys--
endif
load_vars,'rx_model.sav','_model'
if (nkeys eq 0) then GOTO, NO_KEYS


NO_KEYS:
END



;; lir vs z - like in chen17
;; lx vs lir - like in chen17
;; histogram LL
;; histogram NH + modeled
