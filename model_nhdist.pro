PRO model_nhdist, sub_dir, $
                  DATA = data, $
                  NHDIST = nhdist, $
                  SETNH = setnh, $
                  SETRX = setrx, $
                  GROUP = group, $
                  UNIFORM = uniform, $
                  VARIABLE = variable, $
                  FINAL = final, $
                  POSTMOD = postmod
                  

chen = 0
;; check for keyword commands
nkeys = n_elements(data) + $
        n_elements(nhdist) + $
        n_elements(setnh) + $
        n_elements(setrx) + $
        n_elements(group) + $
        n_elements(uniform) + $
        n_elements(variable) + $
        n_elements(final)
if (nkeys eq 0) then GOTO, NO_KEYS

;; assume current directory unless specified
if (n_elements(sub_dir) eq 0) then sub_dir = './' else $
                                   sub_dir = strlowcase(sub_dir)+'/' & file_mkdir,sub_dir

;; pass data from XRAY_LACK_AGN to MODEL_NHDIST
if keyword_set(data) then begin
    pull_carroll21_data
    nkeys--
endif
if (chen eq 1) then load_vars,'carroll21b_data.sav','_data' else $
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
    rx_path = prep_dir+'/rxz_borus_gupta.sav'
    ;; -f force to overwrite symbolic link in case changing
    spawn, 'ln -sf '+rx_path+' rx_conversion.sav'
    nkeys--
endif
load_vars,'rx_conversion.sav','_rxnh'
if (nkeys eq 0) then GOTO, NO_KEYS

;; select data group WISE AGN, secondary sources, or all (WAC/SEC/ALL)
if keyword_set(group) then begin
    if (typename(group) ne 'STRING') then message, 'PLEASE SELECT INPUT GROUP: WAC/SEC/ALL/WAC_HIX/WAC_LOX/OFFSET'
    set_data_group,group=group
    nkeys--
endif
load_vars,'select_group.sav','_group'
if (nkeys eq 0) then GOTO, NO_KEYS

;; estimate the distribution of fixed CTF
if keyword_set(uniform) then begin
    ctf_uniform
    nkeys--
endif
load_vars,'fct_uniform.sav','_uniform'
if (nkeys eq 0) then GOTO, NO_KEYS

;; estimate the distribution of free CTF
if keyword_set(variable) then begin
    ctf_variable,/chisq
    nkeys--
endif
load_vars,'fct_variable.sav','_variable'
if (nkeys eq 0) then GOTO, NO_KEYS

;; simulate the NH and RL distributions
if keyword_set(final) then begin
    model_rxdist,postmod=postmod
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
