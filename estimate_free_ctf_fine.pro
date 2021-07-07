PRO estimate_free_ctf


common _data
common _nhdist
common _nhobs
common _rxnh
common _group
;common _fixed


cvm = 0

;; run this script NITER times and look at the distribution in CTF
niter = 1000;0

;; ROUGH PASS
fctv1 = dblarr(niter)
f24v1 = dblarr(niter)
f25v1 = dblarr(niter)
statv1 = dblarr(6,niter)
nhmv1 = dblarr(nsrc,niter)
rxmv1 = dblarr(nsrc,niter)
iimv1 = bytarr(nsrc,niter)

;; fixed CT fraction
step = 0.10d
fct_rough1 = [step:1.-step:step]
nfrac = n_elements(fct_rough1)

;; free CT fraction split between NH=24-25 and 25-26
f25 = [0.05d:1.-0.05d:0.05d];[step:1.-step:step]
f24 = 1.-f25
nfree = n_elements(f24)

a2_rough1 = dblarr(nfrac,nfree,niter)

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
    a2 = dblarr(nfrac,nfree)
    p_a2 = dblarr(nfrac,nfree)
    rxd_scat = randomn(seed,ndet)*e_rxd
    ;rxm_scat = randomn(seed,nsrc)*rx_scat
    for i = 0,nfrac-1 do begin
        ;; vary CT fraction
        nct = round((nthin/(1.-fct_rough1[i]))*fct_rough1[i])
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
                ;; if comparing test statistics, need p_a2
                a2[i,j] = ad_test(rxd+rxd_scat,rx_mod[idet,i,j],permute=(test eq 'JOINT'),prob=p,cvm=cvm)
                p_a2[i,j] = p
                ;kstwo,rxd,rx_mod[idet,i,j],d,prob
                ;a2[i,j] = d
                ;p_a2[i,j] = prob
            endif else if (detct gt 0) then begin
                a2[i,j] = -1.
                p_a2[i,j] = -1.
            endif else message, 'NO MODELED DETECTIONS.'
        endfor
    endfor
    iia2 = finite(a2) and a2 gt 0.
    a2_rough1[*,*,n] = a2
    
    ;; determine "best-fit"
    ;; QUESTION: method to determine best-fit?
    case strupcase(test) of 
        'AD': begin
            ibest = where(a2 eq min(a2[where(a2 gt 0,/null)]),nbest)
            ;if (nbest gt 1) then stop
            if (nbest gt 1) then ibest = ibest[0]
            stat = [a2[ibest],p_a2[ibest],0.,0.,0.,0.]
            end
        'JOINT': begin
            ;; constraints by X-ray stacked fluxes
            fx_est = estimate_fx(rx_mod,iimod,/iterate)
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
        fctv1[n] = fct_rough1[ind[0]]
        f24v1[n] = f24[ind[1]]
        f25v1[n] = f25[ind[1]]
        statv1[*,n] = stat
        nhmv1[*,n] = nh_mod[*,ind[0],ind[1]]
        rxmv1[*,n] = rx_mod[*,ind[0],ind[1]]
        iimv1[*,n] = iimod[*,ind[0],ind[1]]
    endelse
    ;; progress alert
    if (n eq 0) then begin
        print, '=============================================='
        print, 'BEGIN - FREE FCT, ROUND 1, ROUGH PASS'
    endif else if (n mod (ncount/10) eq 0) then print, '      - '+string(100.*n/ncount,format='(i2)')+'% COMPLETE'
endfor
print, 'END   - FREE FCT, ROUND 1, ROUGH PASS'
print, '=============================================='

sav_vars = ['FCTV1','F24V1','F25V1','STATV1', $
            'NREJ1','NHMV1','RXMV1','IIMV1', $
            'FCT_ROUGH1','A2_ROUGH1']
sav_inds = []

;; FINE PASS
fctv1_fine = dblarr(niter)
f24v1_fine = dblarr(niter)
f25v1_fine = dblarr(niter)
statv1_fine = dblarr(6,niter)
nhmv1_fine = dblarr(nsrc,niter)
rxmv1_fine = dblarr(nsrc,niter)
iimv1_fine = bytarr(nsrc,niter)

;; fixed CT fraction
step = 0.05d
fct_fine1 = [(mode(fctv1)-2.*ceil(stddev(fctv1)/step)*step)>(0.+step):(mode(fctv1)+2.0*ceil(stddev(fctv1)/step)*step)<(1.-step):step]
nfrac = n_elements(fct_fine1)

;; free CT fraction split between NH=24-25 and 25-26
f25 = [0.05d:1.-0.05d:0.05d];[step:1.-step:step]
f24 = 1.-f25
nfree = n_elements(f24)

a2_fine1 = dblarr(nfrac,nfree,niter)

;; counter for iteration alerts
ncount = ceil(niter/10.)*10.
nrej1_fine = 0l
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
    a2 = dblarr(nfrac,nfree)
    p_a2 = dblarr(nfrac,nfree)
    rxd_scat = randomn(seed,ndet)*e_rxd
    ;rxm_scat = randomn(seed,nsrc)*rx_scat
    for i = 0,nfrac-1 do begin
        ;; vary CT fraction
        nct = round((nthin/(1.-fct_fine1[i]))*fct_fine1[i])
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
                ;; if comparing test statistics, need p_a2
                a2[i,j] = ad_test(rxd+rxd_scat,rx_mod[idet,i,j],permute=(test eq 'JOINT'),prob=p,cvm=cvm)
                p_a2[i,j] = p
                ;kstwo,rxd,rx_mod[idet,i,j],d,prob
                ;a2[i,j] = d
                ;p_a2[i,j] = prob
            endif else if (detct gt 0) then begin
                a2[i,j] = -1.
                p_a2[i,j] = -1.                
            endif else message, 'NO MODELED DETECTIONS.'
        endfor
    endfor
    iia2 = p_a2 gt min(p_a2)
    a2_fine1[*,*,n] = a2

    ;; determine "best-fit"
    ;; QUESTION: method to determine best-fit?
    case strupcase(test) of 
        'AD': begin
            ibest = where(a2 eq min(a2[where(a2 gt 0,/null)]),nbest)
            ;if (nbest gt 1) then stop
            if (nbest gt 1) then ibest = ibest[0]
            stat = [a2[ibest],p_a2[ibest],0.,0.,0.,0.]
            end
        'JOINT': begin
            ;; constraints by X-ray stacked fluxes
            fx_est = estimate_fx(rx_mod,iimod,/iterate)
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
        nrej1_fine++
        continue
    ;; record best fit statistics
    endif else begin
        ;; account for extra dimensions
        ind = array_indices(a2,ibest)
        fctv1_fine[n] = fct_fine1[ind[0]]
        f24v1_fine[n] = f24[ind[1]]
        f25v1_fine[n] = f25[ind[1]]
        statv1_fine[*,n] = stat
        nhmv1_fine[*,n] = nh_mod[*,ind[0],ind[1]]
        rxmv1_fine[*,n] = rx_mod[*,ind[0],ind[1]]
        iimv1_fine[*,n] = iimod[*,ind[0],ind[1]]
    endelse
    ;; progress alert
    if (n eq 0) then begin
        print, '=============================================='
        print, 'BEGIN - FREE FCT, ROUND 1, FINE PASS'
    endif else if (n mod (ncount/10) eq 0) then print, '      - '+string(100.*n/ncount,format='(i2)')+'% COMPLETE'
endfor
print, 'END   - FREE FCT, ROUND 1, FINE PASS'
print, '=============================================='

sav_vars = [sav_vars,'FCTV1_FINE','F24V1_FINE','F25V1_FINE','STATV1_FINE', $
                     'NREJ1_FINE','NHMV1_FINE','RXMV1_FINE','IIMV1_FINE', $
                     'FCT_FINE1','A2_FINE1']
sav_inds = [sav_inds]

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="ctf_free.sav"')


END











