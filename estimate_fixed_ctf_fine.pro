PRO estimate_fixed_ctf


common _data
common _nhdist
common _nhobs
common _rxnh
common _group


cvm = 1

;; run this script NITER times and look at the distribution in fct
niter = 1000;0

;; ROUGH PASS
fctv = dblarr(niter)
statv = dblarr(6,niter)
nhmv = dblarr(nsrc,niter)
rxmv = dblarr(nsrc,niter)
iimv = bytarr(nsrc,niter)

;; fixed CT fraction
step = 0.10d
fct_rough = [step:1.-step:step]
nfrac = n_elements(fct_rough)

a2_rough = dblarr(nfrac,niter)

;; counter for iteration alerts
ncount = ceil(niter/10.)*10.
nrej = 0l
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
    a2 = dblarr(nfrac)
    p_a2 = dblarr(nfrac)
    ;cdfmv = dblarr(n_elements(rxd),nfrac)
    rxd_scat = randomn(seed,ndet)*e_rxd
    ;rxm_scat = randomn(seed,nsrc)*rx_scat
    for i = 0,nfrac-1 do begin
        nct = round((nthin/(1.-fct_rough[i]))*fct_rough[i])
        nh_resamp = [nh_samp[ithin],24.+2.*randomu(seed,nct)]
        nresamp = n_elements(nh_resamp)
        nh_mod[*,i] = nh_resamp[randomi(nsrc,nresamp)]
        rx_mod[*,i] = rx2nh(nh_mod[*,i],/rx_out,scat=rx_scat)
        iimod[*,i] = rx_mod[*,i] gt rxl
        idet = where(iimod[*,i] eq 1,detct)
        if (detct ge 5) then begin
            ;; if comparing test statistics, need p_a2
            a2[i] = ad_test(rxd+rxd_scat,rx_mod[idet,i],permute=(test eq 'JOINT'),prob=p,cvm=cvm);,cdfm=cdfm)
            p_a2[i] = p
            ;cdfmv[*,i] = cdfm
            ;kstwo,rxd,rx_mod[idet,i],d,prob
            ;a2[i] = d
            ;p_a2[i] = prob
        endif else if (detct gt 0) then begin
            a2[i] = -1.
            p_a2[i] = -1.
        endif else message, 'NO MODELED DETECTIONS.'
    endfor
    iia2 = finite(a2) and a2 gt 0.
    a2_rough[*,n] = a2
    
    ;; determine "best-fit"
    ;; QUESTION: method to determine best-fit?
    case strupcase(test) of 
        'AD': begin
            ibest = where(a2 eq min(a2[where(iia2,/null)]),nbest)
            ;if (nbest gt 1) then stop
            if (nbest gt 1) then ibest = ibest[0]
            stat = [a2[ibest],p_a2[ibest],0.,0.,0.,0.]
            ;; estimate model fluxes
            ;fx_est = estimate_fx(rx_mod,iimod,/iterate)
            ;if (n eq 0) then fxmv = replicate(dup_struct((fx_est)[0]),niter)
            end
        'JOINT': begin
            ;; estimate model fluxes
            fx_est = estimate_fx(rx_mod,iimod,/iterate)
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
            ;; flag chisqr_pdf value == 1
            ;p_x2 = p_x2 > min(p_x2[where(p_x2,/null)])/2.
            ;iix2 = p_x2 gt min(p_x2)/2.
            iix2 = finite(p_x2) and p_x2 gt 0.
            ;; combined test statistic
            ;; Fisher
            x2_joint = -2.*(alog(p_a2)+alog(p_x2))
            p_joint = 1.-chisqr_pdf(x2_joint,4.)  ;; dof = 2k, where k == 2, the number of tests being combined
            ;; Pearson
            ;x2_joint = -2.*(alog(1.-p_a2)+alog(1.-p_x2))
            ;p_joint = 1.-chisqr_pdf(x2_joint,4.)
            ibest = where(x2_joint eq min(x2_joint[where(iia2 and iix2,/null)]),nbest)
            if (nbest gt 1) then stop
            stat = [a2[ibest],p_a2[ibest],x2[ibest],p_x2[ibest],x2_joint[ibest],p_joint[ibest]]
            end
        else: message, 'IMPROPER INPUT OF TEST METHOD.'
    endcase
    ;; skip if no good fits during iteration
    if (nbest eq 0) then begin
        nrej++
        continue
    ;; record best fit statistics
    endif else begin
        fctv[n] = fct_rough[ibest]
        statv[*,n] = stat
        nhmv[*,n] = nh_mod[*,ibest]
        rxmv[*,n] = rx_mod[*,ibest]
        iimv[*,n] = iimod[*,ibest]
        ;fxmv[n] = fx_est[ibest]
    endelse
    ;; progress alert
    if (n eq 0) then begin
        print, '=============================================='
        print, 'BEGIN - FIXED FCT, ROUGH PASS'
    endif else if ((n ne 0) and (n mod (ncount/10) eq 0)) then print, '      - '+string(100.*n/ncount,format='(i2)')+'% COMPLETE'
endfor
print, 'END   - FIXED FCT, ROUGH PASS'
print, '=============================================='

sav_vars = ['FCTV','STATV','NREJ','NHMV','RXMV','IIMV','FCT_ROUGH','A2_ROUGH']
sav_inds = []

;; FINE PASS
fctv_fine = dblarr(niter)
statv_fine = dblarr(6,niter)
nhmv_fine = dblarr(nsrc,niter)
rxmv_fine = dblarr(nsrc,niter)
iimv_fine = bytarr(nsrc,niter)

;; fixed CT fraction
step = 0.05d
fct_fine = [(mode(fctv)-2.*ceil(stddev(fctv)/step)*step)>(0.+step):(mode(fctv)+2.0*ceil(stddev(fctv)/step)*step)<(1.-step):step]
nfrac = n_elements(fct_fine)
a2_fine = dblarr(nfrac,niter)

;; counter for iteration alerts
ncount = ceil(niter/10.)*10.
nrej_fine = 0l
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
    a2 = dblarr(nfrac)
    p_a2 = dblarr(nfrac)
    ;dv = dblarr(nfrac)
    ;probv = dblarr(nfrac)
    rxd_scat = randomn(seed,ndet)*e_rxd
    ;rxm_scat = randomn(seed,nsrc)*rx_scat
    for i = 0,nfrac-1 do begin
        nct = round((nthin/(1.-fct_fine[i]))*fct_fine[i])
        nh_resamp = [nh_samp[ithin],24.+2.*randomu(seed,nct)]
        nresamp = n_elements(nh_resamp)
        nh_mod[*,i] = nh_resamp[randomi(nsrc,nresamp)]
        rx_mod[*,i] = rx2nh(nh_mod[*,i],/rx_out,scat=rx_scat)
        iimod[*,i] = rx_mod[*,i] gt rxl
        idet = where(iimod[*,i] eq 1,detct)
        if (detct ge 5) then begin
            ;; if comparing test statistics, need p_a2
            a2[i] = ad_test(rxd+rxd_scat,rx_mod[idet,i],permute=(test eq 'JOINT'),prob=p,cvm=cvm)
            p_a2[i] = p
            ;kstwo,rxd,rx_mod[idet,i],d,prob
            ;a2[i] = d
            ;p_a2[i] = prob            
        endif else if (detct gt 0) then begin
            a2[i] = -1.
            p_a2[i] = -1.
        endif else message, 'NO MODELED DETECTIONS.'
    endfor
    iia2 = finite(a2) and a2 gt 0.
    a2_fine[*,n] = a2

    ;; determine "best-fit"
    ;; QUESTION: method to determine best-fit?
    case strupcase(test) of 
        'AD': begin
            ibest = where(a2 eq min(a2[where(iia2,/null)]),nbest)
            ;if (nbest gt 1) then stop
            if (nbest gt 1) then ibest = ibest[0]
            stat = [a2[ibest],p_a2[ibest],0.,0.,0.,0.]
            ;; estimate model fluxes
            ;fx_est = estimate_fx(rx_mod,iimod,/iterate)
            ;if (n eq 0) then fxmv = replicate(dup_struct((fx_est)[0]),niter)
            end
        'JOINT': begin
            ;; estimate model fluxes
            fx_est = estimate_fx(rx_mod,iimod,/iterate)
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
            ;; flag chisqr_pdf value == 1
            ;p_x2 = p_x2 > min(p_x2[where(p_x2,/null)])/2.
            ;iix2 = p_x2 gt min(p_x2)/2.
            iix2 = finite(p_x2) and p_x2 gt 0.
            ;; combined test statistic
            ;; Fisher
            x2_joint = -2.*(alog(p_a2)+alog(p_x2))
            p_joint = 1.-chisqr_pdf(x2_joint,4.)  ;; dof = 2k, where k == 2, the number of tests being combined
            ;; Pearson
            ;x2_joint = -2.*(alog(1.-p_a2)+alog(1.-p_x2))
            ;p_joint = 1.-chisqr_pdf(x2_joint,4.)
            ibest = where(x2_joint eq min(x2_joint[where(iia2 and iix2,/null)]),nbest)
            if (nbest gt 1) then stop
            stat = [a2[ibest],p_a2[ibest],x2[ibest],p_x2[ibest],x2_joint[ibest],p_joint[ibest]]
            end
        else: message, 'IMPROPER INPUT OF TEST METHOD.'
    endcase
    ;; skip if no good fits during iteration
    if (nbest eq 0) then begin
        nrej_fine++
        continue
    ;; record best fit statistics
    endif else begin
        fctv_fine[n] = fct_fine[ibest]
        statv_fine[*,n] = stat
        nhmv_fine[*,n] = nh_mod[*,ibest]
        rxmv_fine[*,n] = rx_mod[*,ibest]
        iimv_fine[*,n] = iimod[*,ibest]
        ;fxmv[n] = fx_est[ibest]
    endelse
    ;; progress alert
    if (n eq 0) then begin
        print, '=============================================='
        print, 'BEGIN - FIXED FCT, FINE PASS'
    endif else if ((n ne 0) and (n mod (ncount/10) eq 0)) then print, '      - '+string(100.*n/ncount,format='(i2)')+'% COMPLETE'
endfor
print, 'END   - FIXED FCT, FINE PASS'
print, '=============================================='

sav_vars = [sav_vars,'FCTV_FINE','STATV_FINE','NREJ_FINE','NHMV_FINE','RXMV_FINE','IIMV_FINE','FCT_FINE','A2_FINE']
sav_inds = [sav_inds]

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="ctf_fixed.sav"')


END



;; difference in uncertainties on model flux estimates is negligible

