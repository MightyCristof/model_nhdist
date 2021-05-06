PRO estimate_free_ctf


common _data
common _nhdist
common _nhobs
common _rxnh
common _group
;common _fixed


;; run this script NITER times and look at the distribution in CTF
niter = 10000
fctv1 = dblarr(niter)
f24v1 = dblarr(niter)
f25v1 = dblarr(niter)
stat_fctv1 = dblarr(6,niter)
nh_modv1 = dblarr(nsrc,niter)

;; fixed CT fraction
step = 0.02d
fct = [step:1.-step:step]
nfrac = n_elements(fct)

;; free CT fraction split between NH=24-25 and 25-26
f25 = [step:1.-step:step]
f24 = 1.-f25
nfree = n_elements(f24)

;; counter for iteration alerts
ncount = ceil(niter/10.)*10.
nrej1 = 0l
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
    nh_mod = dblarr(nsrc,nfrac,nfree)
    rx_mod = dblarr(nsrc,nfrac,nfree)
    iimod = bytarr(nsrc,nfrac,nfree)
    ad = dblarr(2,nfrac,nfree)
    for i = 0,nfrac-1 do begin
        ;; vary CT fraction
        nct = round((nthin/(1.-fct[i]))*fct[i])
        for j = 0,nfree-1 do begin
            n25 = round(nct*f25[j])>1
            n24 = nct-n25
            nh_resamp = [nh_samp[ithin],24.+randomu(seed,n24),25.+randomu(seed,n25)]
            nresamp = n_elements(nh_resamp)
            nh_mod[*,i,j] = nh_resamp[randomi(nsrc,nresamp)]
            rx_mod[*,i,j] = rx2nh(nh_mod[*,i,j],/rx_out,scat=rx_scat)
            iimod[*,i,j] = rx_mod[*,i,j] gt rxl
            idet = where(iimod[*,i,j] eq 1,detct)
            if (detct ge 5) then begin
                ad[*,i,j] = ad_test(rxd,rx_mod[idet,i,j],prob=(test eq 'JOINT'))
            endif else if (detct gt 0) then begin
                ad[*,i,j] = [-1.,-1.,-1.]
            endif else message, 'NO MODELED DETECTIONS.'
        endfor
    endfor
    ;; account for extra dimensions    
    sz = size(ad)
    re = execute('a2 = reform(ad['+strjoin(['0',strarr(sz[0]-1)+'*'],',')+'])')
    re = execute('p_a2 = reform(ad['+strjoin(['1',strarr(sz[0]-1)+'*'],',')+'])')
    iia2 = p_a2 gt min(p_a2)
    ;; determine "best-fit"
    ;; QUESTION: method to determine best-fit?
    case strupcase(test) of 
        'AD': begin
            ibest = where(a2 eq min(a2[where(a2 gt 0,/null)]),nbest)
            if (nbest gt 1) then stop            
            stat = [a2[ibest],p_a2[ibest],0.,0.,0.,0.]
            end
        'JOINT': begin
            ;; constraints by X-ray stacked fluxes
            fx_est = estimate_fx(rx_mod,iimod,/cha,/iterate)
            ;; compare to X-ray stacked fluxes
            ;x2_soft = ((fxstak[1,0]-fx_est.soft)/e_fxstak[1,0])^2.
            ;x2_hard = ((fxstak[1,1]-fx_est.hard)/e_fxstak[1,1])^2.
            ;; uncertainties in quadrature
            unc_soft = abs(fxstak[1,0] * sqrt((e_fxstak[1,0]/fxstak[1,0])^2. + (fx_est.e_soft/fx_est.soft)^2.))
            unc_hard = abs(fxstak[1,1] * sqrt((e_fxstak[1,1]/fxstak[1,1])^2. + (fx_est.e_hard/fx_est.hard)^2.))
            x2_soft = ((fxstak[1,0]-fx_est.soft)/unc_soft)^2.
            x2_hard = ((fxstak[1,1]-fx_est.hard)/unc_hard)^2.
            x2 = x2_soft+x2_hard
            p_x2 = 1.-chisqr_pdf(x2,1.) ;; dof = 1 (2 X-ray data points - 1)
            ;; correct for p-value == 1
            p_x2 = p_x2 > min(p_x2[where(p_x2,/null)])/2.
            iix2 = p_x2 gt min(p_x2)
            ;; combined test statistic
            x2_joint = -2.*(alog(p_a2)+alog(p_x2))<99.
            p_joint = 1.-chisqr_pdf(x2_joint,4.) ;; dof = 2k, where k is the number of tests being combined            
            ibest = where(x2_joint eq min(x2_joint[where(iia2 and iix2,/null)]),nbest)
            if (nbest gt 1) then stop
            stat = [a2[ibest],p_a2[ibest],x2[ibest],p_x2[ibest],x2_joint[ibest],p_joint[ibest]]
            end
        else: message, 'IMPROPER INPUT OF TEST METHOD.'
    endcase    
    ;; skip if no good fits during iteration
    if (nbest eq 0) then begin
        nrej1++
        continue
    ;; record best fit statistics
    endif else begin
        ;; account for extra dimensions
        ind = array_indices(a2,ibest)
        fctv1[n] = fct[ind[0]]
        f24v1[n] = f24[ind[1]]
        f25v1[n] = f25[ind[1]]
        stat_fctv1[*,n] = stat
        nh_modv1[*,n] = nh_mod[*,ind[0],ind[1]]
    endelse
    ;; progress alert
    if (n eq 0) then begin
        print, '=============================================='
        print, 'BEGIN - FREE FCT, ROUND 1'
    endif else if (n mod (ncount/10) eq 0) then print, '      - '+string(100.*n/ncount,format='(i2)')+'% COMPLETE'
endfor
print, 'END   - FREE FCT, ROUND 1'
print, '=============================================='

;; uncertainties on NH bins
nhmv1 = hist2d_avg(nh_modv1,1.)
;sig_nhv1 = nh_mod.sig
;mad_nhv1 = nh_mod.mad

sav_vars = ['FCTV1','F24V1','F25V1','STAT_FCTV1','NREJ1','NHMV1'];'SIG_NHV1','MAD_NHV1']
sav_inds = []

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="ctf_free.sav"')


END










