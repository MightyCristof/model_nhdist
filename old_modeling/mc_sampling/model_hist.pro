function model_hist,xb,a,outsamp=fsampnh

hi_frac = a[0]
hi_sigma = a[1]
lo_frac = a[2]

nsamp=100000.

sampnh = nh_mc(nsamp);,/plt)
;sampnh = reform(samp[*,0])
;cdfsampnh = reform(samp[*,1])

;; number of observations to add
nhi = nsamp*hi_frac
nlo = nsamp*lo_frac

;; create array of observations
nhhi = hi_sigma*randomn(seed,nhi)+ 24.5
nhlo = fltarr(nlo)+19.;0.5*randomn(seed,40000)+ 20.
fsampnh = [nhlo,sampnh,nhhi]

lscat = 0.3
fsampll = ll2nh(fsampnh,/get_lum)+randomn(seed,n_elements(fsampnh)) * lscat

binsz = xb[1]-xb[0]

;; ?????????????????? isnt the output supposed to be the same size of input?

;yhist_out=histogram(fsampll,min=xb[0]-binsz/5.,binsize=binsz,locations=xhist_out)   ;; ????????? /5


yhist_out = histogram(fsampll,nbins=n_elements(xb),binsize=binsz,min=xb[0],locations=xhist_out)

;hist_out = [[nm(yhist_out)],[xhist_out]]
hist_out = nm(yhist_out)

return,hist_out

;plothist,ll[iagn_nox],bin=scott(ll[iagn_nox]),peak=1,col='red',/ov
;plothist,fsampll,bin=scott(fsampll),peak=1,col='purple',/ov



end