PRO pull_carroll21_data


;; Note 1: RXLIM: RX at flux limit, including detected sources. CMC 11-Dec-20 (after paper acceptace)
;dir = '/Users/ccarroll/Research/projects/xray_lack_agn/workspace/run_20201008_final_submission/'
dir = '/Users/ccarroll/Research/projects/xray_lack_agn/workspace/run_20201008_final_updated/'

;; variables for MODEL_NHDIST
vars = ['OBJID','RA','DEC','Z','ZERR', $                    ;; fits.sav
        'PERC_AGN', $                                       ;; resamp.sav
        'IIWAC', $                                          ;; detections_wac.sav
        'EBV','E_EBV','LIR','E_LIR','DL2', $                ;; src_luminosities.sav
        'LOGLIR','E_LOGLIR','LXIR','LOGLXIR', $             ;;
        'IIQUAL_DET','IIQUAL_NON', $                        ;; quality_src.sav
        'LX','E_LX','LOGLX','E_LOGLX', $                    ;; combined_lum.sav
        'LXLIM','E_LXLIM','LOGLXLIM','E_LOGLXLIM', $        ;;
        'LLDET','E_LLDET','LLNON','E_LLNON', $              ;;
        'XDET','XNON','LOGEBV', $                           ;;
        'RXLIM', $                                          ;; combined_lum.sav (UPDATED)
        'SIMRLW','WNH_POWER','WNH_BORUS', $                 ;; surv_anal.sav
        'FLUXV','FLUX_ERRV','LXV']                          ;; stack_output.sav

;; IDL save files containing above variables
files = ['fits.sav','resamp.sav','detections_wac.sav','src_luminosities.sav', $
         'quality_src.sav','combined_lum.sav','surv_anal.sav','stack_output/stack_output.sav']

for i = 0,n_elements(files)-1 do restore,dir+files[i]

;; correct RCH output to have BAGN then GAL
isort = [0,1,2,3,5,4]       ;; sort on Carroll+20 Table 6 [WDET,WNON,RDET,RNON,BAGN,GAL]
xsubg = ['WDET','WNON','RDET','RNON','BAGN','GAL']
fluxv = fluxv[isort,2]      ;; 2 == ENRANGEV of 2-7keV, fluxes scaled to 2-10keV for comparison
flux_errv = flux_errv[isort,2]
lxv = lxv[isort,2]

;; trim analysis subset
iq = where(iiqual)
for i = 0,n_elements(vars[0:-7])-1 do re = execute(vars[i]+' = '+vars[i]+'[iq]')
;; LL to RL
rl = ['RLDET','E_RLDET','RLNON','E_RLNON']
irl = where(strmatch(vars,'*LL*'))
for i = 0,n_elements(rl)-1 do re = execute(rl[i]+' = '+vars[irl[i]])
vars = [vars[0:irl[0]-1],rl,vars[irl[-1]+1:-1],'XSUBG']     ;; add X-ray stacking sub grouping

var_str = strjoin(vars,',')
re = execute('save,'+var_str+',file="carroll21_data.sav"')


END




