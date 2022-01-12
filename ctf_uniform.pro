PRO ctf_uniform


common _data
common _nhdist
common _nhobs
common _rxnh
common _group


;; AD vs CvM test
cvm = 0

;; run this script NITER times and look at the distribution in fct
niter = 1000

;; uniform CT fraction
step = 0.01d
fct_steps = [step:1.-step:step]
nfrac = n_elements(fct_steps)

;; data variables
rxdv = dblarr(ndet,niter)

;; modeling variables
fctv = dblarr(niter)
statv = dblarr(6,niter)
nhmv = dblarr(nsrc,niter)
rxmv = dblarr(nsrc,niter)
iimv = bytarr(nsrc,niter)
a2v = dblarr(nfrac,niter)

;; counter for iteration alerts
ncount = ceil(niter/10.)*10.
;; run the entire process NITER times
for n = 0,niter-1 do begin
    ;; resample observed NH distribution to increase density for draw
    nsamp = nsrc*100.
    nh_samp = nh_mc(nh_obs,nsamp)
    ;; number of Compton-thin sources
    ithin = where(nh_samp lt 24.,nthin);,complement=ithick,ncomplement=nthick)
    fthin = nthin/nsamp
    ;; iterative modeling arrays
    nh_mod = dblarr(nsrc,nfrac)
    rx_mod = dblarr(nsrc,nfrac)
    iimod = bytarr(nsrc,nfrac)
    a2 = dblarr(nfrac)
    p_a2 = dblarr(nfrac)
    rxdv[*,n] = rxd+randomn(seed,ndet)*0.23;rx_scat
    ;rxdv[*,n] = rxd+randomn(seed,ndet)*e_rxd
    ;; loop over CT fracions
    for i = 0,nfrac-1 do begin
        nct = round((nthin/(1.-fct_steps[i]))*fct_steps[i])     ;; determine the number of CT sources
        nh_resamp = [nh_samp[ithin],24.+2.*randomu(seed,nct)]   ;; create an NH distribution with specified CT fraction
        nresamp = n_elements(nh_resamp)                         ;; number of sources in the NH distribution
        nh_mod[*,i] = nh_resamp[randomi(nsrc,nresamp)]          ;; resample the NH distribution to match our number of sources
        rx_mod[*,i] = rx2nh(nh_mod[*,i],/rx_out,scat=rx_scat)   ;; convert model NH to model RX
        iimod[*,i] = rx_mod[*,i] gt rxl                         ;; simulate model detections/non-detections comparing to RX limits (from instrument flux limits)
        idet = where(iimod[*,i] eq 1,moddet,ncomplement=modnon) ;; flag model detections
        rxmn = dblarr(modnon)-9999.                             ;; array for model non-detections
        ;; determine test statistic
        if (moddet ge 5) then begin
            a2[i] = ad_test([rxdn,rxdv[*,n]],[rxmn,rx_mod[idet,i]],permute=0,prob=p,weight=1)
            p_a2[i] = p
        endif else if (moddet gt 0) then begin
            a2[i] = -1.
            p_a2[i] = -1.
        endif else message, 'NO MODELED DETECTIONS.'
    endfor
    ;; finite values only
    iia2 = finite(p_a2) and p_a2 gt 0.
    a2v[*,n] = a2
    ;; determine "best-fit" model per iteration
    ibest = where(a2 eq min(a2[where(iia2,/null)]),nbest)
    if (nbest gt 1) then ibest = ibest[0]
    stat = [a2[ibest],p_a2[ibest],0.,0.,0.,0.]
    ;; halt if no best fit found
    if (nbest eq 0) then stop
    ;; record best fit statistics
    fctv[n] = fct_steps[ibest]
    statv[*,n] = stat
    nhmv[*,n] = nh_mod[*,ibest]
    rxmv[*,n] = rx_mod[*,ibest]
    iimv[*,n] = iimod[*,ibest]
    ;; progress alert
    if (n eq 0) then begin
        print, '=============================================='
        print, 'BEGIN - UNIFORM CTF'
    endif else if ((n ne 0) and (n mod (ncount/10) eq 0)) then print, '      - '+string(100.*n/ncount,format='(i2)')+'% COMPLETE'
endfor
print, 'END   - UNIFORM CTF'
print, '=============================================='

sav_vars = ['FCT_STEPS', $
            'RXDV', $
            'FCTV','STATV','NHMV','RXMV','A2V']
sav_inds = ['IIMV']

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="fct_uniform.sav"')


END



;; difference in uncertainties on model flux estimates is negligible

