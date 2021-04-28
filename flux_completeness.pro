PRO flux_completeness

pushd,'/Users/ccarroll/Research/projects/model_nhdist/workspace'
load_vars,'carroll21_data.sav','_data'
load_vars,'observed_nhdist.sav','_nhdist'
pushd,'variant_models/lan_ad'
load_vars,'select_nhobs.sav','_nhobs'
load_vars,'rx_conversion.sav','_rxnh'
load_vars,'select_group.sav','_group'
popd
popd

restore,'idlsave.dat'

runs = 'fx_est'+strtrim(indgen(100),2)
nruns = n_elements(runs)

;; estimate flux with contamination
for i = 0,nruns-1 do re = execute(runs[i]+' = test_fx_est_contamination(rx_mod,iimod,/cha,/iterate)')

;; concatenate fluxes into arrays
fullv = dblarr(nfrac,nruns)
hardv = dblarr(nfrac,nruns)
softv = dblarr(nfrac,nruns)
cfullv = dblarr(nfrac,nruns)
chardv = dblarr(nfrac,nruns)
csoftv = dblarr(nfrac,nruns)
for i = 0,nruns-1 do begin
    re = execute('fullv[*,i] = '+runs[i]+'.full')
    re = execute('hardv[*,i] = '+runs[i]+'.hard')
    re = execute('softv[*,i] = '+runs[i]+'.soft')
    re = execute('cfullv[*,i] = '+runs[i]+'.cfull')
    re = execute('chardv[*,i] = '+runs[i]+'.chard')
    re = execute('csoftv[*,i] = '+runs[i]+'.csoft')
endfor
;; estimate over FCT
full = mean(fullv,dim=2)
hard = mean(hardv,dim=2)
soft = mean(softv,dim=2)
cfull = mean(cfullv,dim=2)
chard = mean(chardv,dim=2)
csoft = mean(csoftv,dim=2)
;full = median(fullv,dim=2)
;hard = median(hardv,dim=2)
;soft = median(softv,dim=2)
;cfull = median(cfullv,dim=2)
;chard = median(chardv,dim=2)
;csoft = median(csoftv,dim=2)
;full = dblarr(nfrac)
;hard = dblarr(nfrac)
;soft = dblarr(nfrac)
;for i = 0,nfrac-1 do begin
;    full[i] = mode(reform(fullv[i,*]),bin=kde_bandwidth(fullv[i,*]))
;    hard[i] = mode(reform(hardv[i,*]),bin=kde_bandwidth(hardv[i,*]))
;    soft[i] = mode(reform(softv[i,*]),bin=kde_bandwidth(softv[i,*]))
;    cfull[i] = mode(reform(cfullv[i,*]),bin=kde_bandwidth(cfullv[i,*]))
;    chard[i] = mode(reform(chardv[i,*]),bin=kde_bandwidth(chardv[i,*]))
;    csoft[i] = mode(reform(csoftv[i,*]),bin=kde_bandwidth(csoftv[i,*]))
;endfor

;; compare to flux estimates with no removed sources
fr = cfull/full
hr = chard/hard
sr = csoft/soft

save,file='idlsave2.dat'



END