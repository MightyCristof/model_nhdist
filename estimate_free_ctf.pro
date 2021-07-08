PRO estimate_free_ctf


common _data
common _nhdist
common _nhobs
common _rxnh
common _group
;common _fixed


cvm = 0

;; run this script NITER times and look at the distribution in CTF
niter = 1000

;; modeling variables
fctv1 = dblarr(niter)
f24v1 = dblarr(niter)
f25v1 = dblarr(niter)
statv1 = dblarr(6,niter)
nhmv1 = dblarr(nsrc,niter)
rxmv1 = dblarr(nsrc,niter)
iimv1 = bytarr(nsrc,niter)
;; data variables
rxdv1 = dblarr(ndet,niter)

;; fixed CT fraction
step = 0.05d
fct_free1 = [step:1.-step:step]
nfrac = n_elements(fct_free1)

;; free CT fraction split between NH=24-25 and 25-26
f25 = [0.05d:1.-0.05d:0.05d];[step:1.-step:step]
f24 = 1.-f25
nfree = n_elements(f24)

a2_free1 = dblarr(nfrac,nfree,niter)

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
    mdetf = dblarr(nfrac,nfree)
    a2 = dblarr(nfrac,nfree)
    p_a2 = dblarr(nfrac,nfree)
    rxdv1[*,n] = rxd+randomn(seed,ndet)*0.23;rx_scat
    ;rxdv1[*,n] = rxd+randomn(seed,ndet)*e_rxd
    for i = 0,nfrac-1 do begin
        ;; vary CT fraction
        nct = round((nthin/(1.-fct_free1[i]))*fct_free1[i])
        tdetf = dblarr(nfree)       ;; temp detected fraction
        for j = 0,nfree-1 do begin
            n25 = round(nct*f25[j])>1
            n24 = nct-n25
            nh_resamp = [nh_samp[ithin],24.+randomu(seed,n24),25.+randomu(seed,n25)]
            nresamp = n_elements(nh_resamp)
            nh_mod[*,i,j] = nh_resamp[randomi(nsrc,nresamp)]
            rx_mod[*,i,j] = rx2nh(nh_mod[*,i,j],/rx_out,scat=rx_scat)
            iimod[*,i,j] = rx_mod[*,i,j] gt rxl
            idet = where(iimod[*,i,j] eq 1,moddet)
            mdetf[i,j] = 1.*moddet/nsrc
            if (moddet ge 5) then begin
                ;; if comparing test statistics, need p_a2
                a2[i,j] = ad_test(rxdv1[*,n],rx_mod[idet,i,j],permute=(test eq 'JOINT'),prob=p,cvm=cvm)
                p_a2[i,j] = p
                ;kstwo,rxd,rx_mod[idet,i,j],d,prob
                ;a2[i,j] = d
                ;p_a2[i,j] = prob
            endif else if (moddet gt 0) then begin
                a2[i,j] = -1.
                p_a2[i,j] = -1.
            endif else message, 'NO MODELED DETECTIONS.'
        endfor
    endfor
    ;; weight A2 test statistic by fractional detections
    mdetf = mean(mdetf,dim=2)
    dweight = ((mdetf-ddetf)/ddetf)^2.
    dw = rebin(dweight/total(dweight),nfrac,nfree)
    pw = a2/(a2+dw)
    a2 += dw
    p_a2 *= pw
    
    ;; penalize large disparity in NH split
    spenalty = (f24-f25)^2./nfree
    sp = rebin(reform(spenalty,1,nfree),nfrac,nfree)
    pp = a2/(a2+sp)
    a2 += sp
    p_a2 *= pp    
    
    ;; finite values only
    iia2 = finite(p_a2) and p_a2 gt 0.
    a2_free1[*,*,n] = a2
    
    ;; determine "best-fit"
    ;; QUESTION: method to determine best-fit?
    case strupcase(test) of 
        'AD': begin
            ibest = where(a2 eq min(a2[where(a2 gt 0,/null)]),nbest)
            if (nbest gt 1) then ibest = ibest[0]
            stat = [a2[ibest],p_a2[ibest],0.,0.,0.,0.]
            end
        'JOINT': begin
            ;; constraints by X-ray stacked fluxes
            fx_est = estimate_fx(rx_mod,iimod,/iterate)
            ;; X-ray stack data minus model
            del_soft = abs(fxstak[1,0]-fx_est.csoft)
            del_hard = abs(fxstak[1,1]-fx_est.chard)
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
            iix2 = finite(p_x2) and p_x2 gt 0.
            ;; combined test statistic
            x2_joint = -2.*(alog(p_a2)+alog(p_x2))
            p_joint = dblarr(nfrac,nfree)
            p_joint[where(iia2 and iix2,/null)] = 1.-chisqr_pdf(x2_joint[where(iia2 and iix2,/null)],4.)  ;; dof = 2k, where k == 2, the number of tests being combined
            ;; choose best model
            ibest = where(x2_joint eq min(x2_joint[where(iia2 and iix2,/null)]),nbest)
            if (nbest gt 1) then ibest = ibest[0]
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
        fctv1[n] = fct_free1[ind[0]]
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
        print, 'BEGIN - FREE FCT, ROUND 1'
    endif else if (n mod (ncount/10) eq 0) then print, '      - '+string(100.*n/ncount,format='(i2)')+'% COMPLETE'
endfor
print, 'END   - FREE FCT, ROUND 1'
print, '=============================================='

sav_vars = ['FCTV1','F24V1','F25V1','STATV1', $
            'NREJ1','NHMV1','RXMV1','RXDV1', $
            'FCT_FREE1','A2_FREE1']
sav_inds = ['IIMV1']

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="ctf_free.sav"')


END











