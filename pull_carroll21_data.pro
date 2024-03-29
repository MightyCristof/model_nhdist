PRO pull_carroll21_data


;; if not specified, use all sources
if (n_elements(mode) eq 0) then mode = 'ALL'

;; Note 1: RXLIM: RX at flux limit, including detected sources. CMC 11-Dec-20 (after paper acceptace)
;dir = '/Users/ccarroll/Research/projects/xray_lack_agn/workspace/run_20201008_final_submission/'
dir = '/Users/ccarroll/Research/projects/xray_lack_agn/workspace/run_20201008_final_updated/'

;; variables for MODEL_NHDIST
fits = ['OBJID','RA','DEC','Z','ZERR','ZSTR','ZTYPE', $                     ;; fits.sav
        'MAG','E_MAG','FLUX','E_FLUX','BIN', $                              ;;
        'PERC_AGN', $                                                       ;; resamp.sav
        'SDST_DET','SDST_NON', $                                            ;; infield_cha.sav,infield_xmm.sav,infield_nst.sav
        'IIWAC', $                                                          ;; detections_wac.sav
        'EBV','E_EBV','DL2', $                                              ;; src_luminosities.sav
        'LIR','E_LIR','LOGLIR','E_LOGLIR', $                                ;; 
        'FIR','E_FIR','LOGFIR','E_LOGFIR', $                                ;; 
        'LXIR','LXIR_SCAT','LOGLXIR','LOGLXIR_SCAT','FXIR','LOGFXIR', $     ;; 
        'IIQUAL_DET','IIQUAL_NON', $                                        ;; quality_src.sav
        'LX','E_LX','LOGLX','E_LOGLX', $                                    ;; combined_lum.sav
        'FX','E_FX','LOGFX','E_LOGFX', $                                    ;; 
        'LXLIM','E_LXLIM','LOGLXLIM','E_LOGLXLIM', $                        ;; 
        'FXLIM','E_FXLIM','LOGFXLIM','E_LOGFXLIM', $                        ;; 
        'LLDET','E_LLDET','LLNON','E_LLNON', $                              ;; 
        'XDET','XNON','LOGEBV', $                                           ;; 
        'RXLIM']                                                            ;; 
surv = ['SIMRLW','WNH_POWER','WNH_BORUS']                                   ;; surv_anal.sav
stak = ['ENRANGEV','FLUXV','FLUX_ERRV','LXV']                               ;; stack_output.sav

;; IDL save files containing above variables
files = ['fits.sav','resamp.sav','infield_cha.sav','infield_xmm.sav','infield_nst.sav', $
         'detections_wac.sav','src_luminosities.sav', $
         'quality_src.sav','combined_lum.sav','surv_anal.sav','stack_output/stack_output.sav']
;; restore all data
for i = 0,n_elements(files)-1 do restore,dir+files[i]

;; create combined array of off-axis angle for all three X-ray fields
fld = ['CHA','XMM','NST']
sdst_det = dblarr(size(xdet,/dim))-9999.
sdst_non = dblarr(size(xnon,/dim))-9999.
for i = 0,n_elements(fld)-1 do begin
    re = execute('sdst_det[where(xdet eq fld[i],/null)] = sdst_'+fld[i]+'[where(xdet eq fld[i],/null)]')
    re = execute('sdst_non[where(xnon eq fld[i],/null)] = sdst_'+fld[i]+'[where(xnon eq fld[i],/null)]')
endfor

;;; extract hard (2-7 keV) and soft (0.5-2 keV) band flux for Chandra sources
;chandra_hard_soft,[[fpl90s],[fpl90s_err]],[[fpl90m],[fpl90m_err]],[[fpl90h],[fpl90h_err]],softx,hardx
;fxs_cha = softx[*,0] & e_fxs_cha = softx[*,1]
;fxh_cha = hardx[*,0] & e_fxh_cha = hardx[*,1]
;fits = [fits,'FXS_CHA','E_FXS_CHA','FXH_CHA','E_FXH_CHA']

;; correct RCH output to have BAGN then GAL
isort = [0,1,2,3,5,4]       ;; sort on Carroll+20 Table 6 [WDET,WNON,RDET,RNON,BAGN,GAL]
xsubg = ['WDET','WNON','RDET','RNON','BAGN','GAL']
fluxv = fluxv[isort,*]
flux_errv = flux_errv[isort,*]
lxv = lxv[isort,*]
;; correct energy ranges [0.5-2,2-7,2-10] keV
;; remember RCH scaled 2-7 to 2-10 keV for comparison
;; 2-10 to 2-7 keV = 6.888E-01
ien = [1,0,2]
estak = enrangev[ien]
estak[1] = '27'
estak[2] = '210'
fxstak = fluxv[*,ien]
fxstak[*,1] = fxstak[*,2]*6.888E-01
e_fxstak = flux_errv[*,ien]
e_fxstak[*,1] = e_fxstak[*,2]*6.888E-01
logfxstak = alog10(fxstak)
e_logfxstak = e_fxstak/(alog(10.)*fxstak)
lxstak = lxv[*,ien]
lxstak[*,1] = lxstak[*,2]*6.888E-01
;stak  = ['XSUBG',stak]
stak = ['XSUBG','ESTAK','FXSTAK','E_FXSTAK','LOGFXSTAK','E_LOGFXSTAK']

;; source found outside Chandra field (Tonima)
iiflag = ra eq 122.32634735107422 and dec eq 28.25208282470703

;; trim analysis subset
iq = where(iiqual and ~iiflag)
for i = 0,n_elements(fits)-1 do begin
    re = execute('ndim = size('+fits[i]+',/n_dim)')
    if (ndim eq 1) then re = execute(fits[i]+' = '+fits[i]+'[iq]') else $
    if (ndim eq 2) then re = execute(fits[i]+' = '+fits[i]+'[*,iq]')
endfor
;; LL to RL
rx = ['RXDET','E_RXDET','RXNON','E_RXNON']
ill = where(strmatch(fits,'*LL*'),llct)
if (llct gt 0) then begin
    for i = 0,llct-1 do re = execute(rx[i]+' = '+fits[ill[i]])
    fits = [fits[0:ill[0]-1],rx,fits[ill[-1]+1:-1]]
endif

;rxdet[where(rxdet ne -9999)] = rxdet[where(rxdet ne -9999)]+loglxir[where(rxdet ne -9999)]-lxir_c17(loglir[where(rxdet ne -9999)]-0.2)

;; STDDEV observed in LX-LMIR relation of Chen+17
rx_scat = 0.23
;; separate WISE AGN, detections and non-detections, type 1/2
iwac = where(iiwac)
iiwd = iiwac and xdet ne ''
iiwn = iiwac and xnon ne ''
iisd = ~iiwac and xdet ne ''
iisn = ~iiwac and xnon ne ''
iiad = xdet ne ''
iian = xnon ne ''
iihix = loglxir gt 43.6
iilox = loglxir lt 43.6
iitype1 = rxdet gt 0. or rxnon gt 0.
iitype2 = (rxdet gt -9999. and rxdet lt 0.) or (rxnon gt -9999. and rxnon lt 0.); same as ~iitype1

;; add to variable list
fits = [fits,'RX_SCAT','IIWD','IIWN','IISD','IISN','IIAD','IIAN','IIHIX','IILOX','IITYPE1','IITYPE2','BAND']

;; save
var_str = strjoin([fits,surv,stak],',')
re = execute('save,'+var_str+',file="carroll21_data.sav"')


END




