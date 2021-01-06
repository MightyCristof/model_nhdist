


plothist,ll[iagn_nux],bin=scott(ll[iagn_nux]),xhist,yhist,/noplot
e_yhist = sqrt(yhist)/max(yhist)
yhist = nm(yhist)

;plothist,ll[iagn_nux],bin=.3,xhist,yhist
;yhist_norm = (yhist*1d)/max(yhist)
;yhist_err = sqrt(yhist)/max(yhist)/5.

;parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
;                       limits:[0.D,0]}, 3)
;parinfo[*].value = [0.25,0.75,0.25]

;guessp = [0.,0.75,0.25]

;guessp= [0.7,0.3,0.25]
;guessp = [0.01,0.3,0.01]
guessp = [2e-5,0.3,2e-5]

;parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                       ;limits:[0 
;parinfo[0].limits(0) = 0.d
;parinfo[0]. 


p = mpfitfun('model_hist', xhist,yhist,e_yhist,guessp)
outy = model_hist(xhist,p,outsamp=outsamp)

;plothist,ll[iagn_nox],bin=scott(ll[iagn_nox]),peak=1,col='red;',/ov
;plothist,fsampll,bin=scott(fsampll),peak=1,col='purple',/ov

;; ???????????????????
;outy= model_hist(xhist,p,outsamp=outsamp)

;ploterror,xhist,yhist_norm,yhist_err
;oplot,xhist,model_hist(xhist,p),col=fscr()

sampnh = nh_mc(100000)
sampll = ll2nh(sampnh,/get_lum) + randomn(seed,n_elements(sampll))*lscat
plothist,sampll,xll,hll,bin=scott(sampll),peak=1,/noplot

plt  = histplot(xll,hll,/fill_b)
plt1 = histplot(xhist,yhist,e_yhist,col='dodger blue',/fill_b,/ov,name='AGN Nu-X')
plt2 = histplot(xhist,outy,col='teal',/fill_b,/ov,name='Corr. L15 Model ')
lg = legend(target=[plt,plt1,plt2],sample_width=0,/auto_text_color,position=[0.95,0.85])


END

;; lets figure this thing out

;; test 1
;plothist,ll[iagn_nux],bin=scott(ll[iagn_nux]),xhist,yhist,/noplot
;e_yhist = sqrt(yhist)/max(yhist)
;yhist = nm(yhist)
;guessp= [0.7,0.3,0.25]
;p = mpfitfun('model_hist',xhist,yhist,e_yhist,guessp);;;

;; is ll2nh working correctly?
;plothist,sampnh,xsampnh,hsampnh,peak=1
;sampll = ll2nh(sampnh,/get_lum)
;plothist,sampll,xsampll,hsampll,peak=1





;start = [0.5d,0.1d,0.5d]
;pi = replicate({fixed:0, limited:[1,1], limits:[0.1d,0.9d]},3)
;mp = mpfitfun('model_hist',xllnu,hllnu_nm,e_hllnu,start,parinfo=pi)
;mod_out = model_hist(xllnu,mp,outsamp=outsampnh)
;outsampy = mod_out[*,0]
;outsampx = mod_out[*,1]
;outsampll = ll2nh(outsampnh,/reverse_calc)+randomn(seed,n_elements(outsampnh))*0.2
;plt = errorplot(xllnu,hllnu_nm,e_hllnu,'o',color='purple',linestyle='-',xr=[-2.,1.],yr=[0.,1.2],/stairstep)
;plt = plot(outsampx,outsampy,/ov)
;stop
;step = 0.1
;lim = 1./step
;bnorm = 99999.
;for i = 1,lim-1 do begin
;	for j = 1,lim-1 do begin
;		for k = 1,lim-1 do begin
;			guessp = step*[i,j,k]
;			par = mpfitfun('model_hist',xllnu,hllnu_nm,e_hllnu,guessp,bestnorm=fnorm)
;			if (fnorm lt bnorm) then begin
;				bnorm = fnorm
;				bestp = par
;			endif
;		endfor
;	endfor
;endfor
;IDL> print, bnorm
;       62.047083
;IDL> print, bestp
;     0.050000000      0.80000003     0.099999993

;outmod = model_hist(xllnu,bestp,outsamp=outmodnh)
;houtmodll = outmod[*,0]
;xoutmodll = outmod[*,1]



