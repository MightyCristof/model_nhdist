PRO estimate_fixed_ctf, TEST = test


common _data
common _nhdist
common _nhobs
common _rxnh
common _group

;; AD test or JOINT (Fisher method)
if (n_elements(test) eq 0) then message, 'NO TEST STAT SPECIFIED.'

;; run this script NITER times and look at the distribution in fct
niter = 1000
fctv = dblarr(niter)
stat_fctv = dblarr(6,niter)

;; fixed CT fraction
step = 0.01d
fct = [step:1.-step:step]
nfrac = n_elements(fct)

;; counter for iteration alerts
ncount = ceil(niter/10.)*10.

imin1 = lonarr(100)
imin2 = lonarr(100)

for n = 0,niter-1 do begin
    ;; resample observed NH distribution to increase data density
    nsamp = nsrc*100.
    nh_samp = nh_mc(nh_obs,nsamp)
    ;; number of Compton-thin sources
    ithin = where(nh_samp lt 24.,nthin);,complement=ithick,ncomplement=nthick)
    fthin = nthin/nsamp
    
    ;; NH_RESAMP: structure of increased CT sources
    ;; NH_MOD: draw N sources from NH_RESAMP
    ;; RX_MOD: convert the model NH model to model RX with observed scatter in Chen+17 LX-LMIR
    ;; IIMOD: where RX model > RX limits from Carroll+21
    ;; AD: Anderson-Darling test statistic    
    nh_mod = dblarr(nsrc,nfrac)
    rx_mod = dblarr(nsrc,nfrac)
    iimod = bytarr(nsrc,nfrac)
    ad = dblarr(2,nfrac)
    ac = dblarr(7,nfrac)
    for i = 0,nfrac-1 do begin
        nct = round((nthin/(1.-fct[i]))*fct[i])
        nh_resamp = [nh_samp[ithin],24.+2.*randomu(seed,nct)]
        nresamp = n_elements(nh_resamp)
        nh_mod[*,i] = nh_resamp[randomi(nsrc,nresamp)]
        rx_mod[*,i] = rx2nh(nh_mod[*,i],/rx_out,scat=rx_scat)
        iimod[*,i] = rx_mod[*,i] gt rxl
        idet = where(iimod[*,i] eq 1,detct)
        if (detct ge 5) then begin
            ad[*,i] = ad_test(rxd,rx_mod[idet,i],/prob)
        endif else if (detct gt 0) then begin
            ad[*,i] = [-1.,-1.]
        endif else message, 'NO MODELED DETECTIONS.'
    endfor
    a2 = reform(ad[0,*])
    p_a2 = reform(ad[1,*])
    ;; estimate model fluxes
    fx_est = estimate_fx(rx_mod,iimod,/cha,/iterate)
    ;; X-ray stack data-model
    del_soft = fxstak[1,0]-fx_est.csoft
    del_hard = fxstak[1,1]-fx_est.chard
    ;; uncertainties
    ;sig_soft = e_fxstak[1,0]
    ;sig_hard = e_fxstak[1,1]
    sig_soft = abs(fxstak[1,0] * sqrt((e_fxstak[1,0]/fxstak[1,0])^2. + (fx_est.e_csoft/fx_est.csoft)^2.))
    sig_hard = abs(fxstak[1,1] * sqrt((e_fxstak[1,1]/fxstak[1,1])^2. + (fx_est.e_chard/fx_est.chard)^2.))
    ;; chi-square
    x2_soft = (del_soft/sig_soft)^2.
    x2_hard = (del_hard/sig_hard)^2.    
    x2 = x2_soft+x2_hard
    p_x2 = 1.-chisqr_pdf(x2,1.) ;; dof = 1 (2 X-ray data points - 1)
    ;; correct for p-value == 1
    pmin = min(p_x2[where(p_x2,null)])
    p_x2 = p_x2 > (pmin/10.)
    ;; combined test statistic
    x2_joint = -2.*(alog(p_a2)+alog(p_x2))
    p_joint = 1.-chisqr_pdf(x2_joint,4.)  ;; dof = 2k, where k == 2, the number of tests being combined
    ;; determine "best-fit"
    ;; QUESTION: method to determine best-fit?
    case strupcase(test) of 
        'AD': ibest = where(a2 eq min(a2[where(a2 gt 0,/null)]),nbest)
        'JOINT': ibest = where(x2_joint eq min(x2_joint) and a2 ge 0.,nbest)
        else: message, 'IMPROPER INPUT OF TEST METHOD.'
    endcase
    if (nbest gt 0) then ibest = ibest[where(fct[ibest] eq median(fct[ibest]))]
    fctv[n] = fct[ibest]
    stat_fctv[*,n] = [ad[*,ibest],x2[ibest],p_x2[ibest],x2_joint[ibest],p_joint[ibest]]
    ;; progress alert
    if (n eq 0) then begin
        print, '=============================================='
        print, 'BEGIN - FIXED FCT'
    endif else if ((n ne 0) and (n mod (ncount/10) eq 0)) then print, '      - '+string(100.*n/ncount,format='(i2)')+'% COMPLETE'
endfor
print, 'END   - FIXED FCT'
print, '=============================================='

sav_vars = ['FCTV','STAT_FCTV']
sav_inds = []

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="ctf_fixed.sav"')


END







