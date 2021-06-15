FUNCTION ad_test, data, $
                  model, $
                  PROB = p

                  
FORWARD_FUNCTION ad_test

sz = n_elements(data)
if sz lt 2 then message, 'data must contain at least two elements'

ndata = n_elements(data)
nmodel = n_elements(model)

edf = edf(data,xloc=xdata)
edf_model = edf(model,xloc=xmodel)
cdf = interpol(edf_model,xmodel,xdata)
;; correct interpolation outside of data bounds
ilt = where(cdf lt 0.,nlt)
if (nlt gt 0) then cdf[ilt] = interpol([0.,cdf[ilt[-1]+1]],[xdata[ilt[0]],xdata[ilt[-1]+1]],xdata[ilt])
igt = where(cdf gt 1.,ngt)
if (ngt gt 0) then cdf[igt] = interpol([cdf[igt[0]-1],1.],[xdata[igt[0]-1],xdata[igt[-1]]],xdata[igt])
delta = edf - cdf

;ks = max(abs(delta), /nan)
;ky = max(delta, /nan) + max(-delta, /nan)
ad = delta^2 / (cdf * (1 - cdf))
ad = total(ad, /nan) / total(finite(ad))
;mad = total(abs(delta), /nan) / total(finite(delta))
stop
if keyword_set(prob) then begin
    nreps = 1000
    ad_perm = dblarr(nreps)
    joint = [data,model]
    njoint = n_elements(joint)
    ;for i = 0,nreps-1 do ad_perm[i] = ad_test(data,joint[randomi(nmodel,njoint,/nodup)])

    ;; permutation test (well, my quick approximation of one anyway)
    for i = 0,nreps-1 do begin
        pjoint = joint[randomi(njoint,njoint,/nodup)]
        pdata = pjoint[0:ndata-1]
        pmodel = pjoint[ndata:-1]
        ad_perm[i] = ad_test(pdata,pmodel)
    endfor

    p = total(ad_perm gt ad)/nreps > 1./nreps
    ad = [ad,p]
endif

return, ad


END

                  

