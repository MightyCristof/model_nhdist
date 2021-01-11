FUNCTION mc_rldist, arr, ndraw


;; take XY from the input NH structure
yh = histogram(arr,bin=scott(arr),locations=xh)
;; create finer grid for x-axis (simulate continuous function)
dx = [xh[0]:xh[-1]:width(minmax(xh))/9999.]
;; interpolate and normalize to create probability distribution function
pdf = spline(xh,yh,dx)
pdf >= 0.
pdf /= total(pdf)
;; sum to create cumulative distribution function
cdf = total(pdf,/cumulative)
;; draw random values from 0-1 from CDF
draw_cdf = randomu(seed,ndraw)
;; project CDF back to PDF
pdf_samp = interpol(dx,cdf,draw_cdf)

return, pdf_samp


END




;;; create finer grid
;dx = [xnon[0]:xnon[-1]:diff(minmax(xnon))/999.]
;;; interpolate grid to create PDF
;pdf = spline(xnon,ynon,dx) > 0.
;pdf /= total(pdf)
;;; sum PDF to get CDF
;cdf = total(pdf,/cumulative)
;;; unique elements of CDF
;iu = uniq(cdf)
;udx = dx[iu]
;cdf = cdf[iu]
;;; set seed
;if keyword_set(seed) then seed = seed
;;; draw random values from 0-1 from CDF
;draw_cdf = randomu(seed,nobj)
;;; project CDF back to PDF
;;; remember newY=interpol(Y,X,newX)
;samp_pdf = interpol(udx,cdf,draw_cdf)
;;; histogram sampled PDF
;ypdf = histogram(samp_pdf,locations=xpdf,bin=binsz,min=xnon[0],max=xnon[-1])
;
;if keyword_set(plt) then begin
;    p = plot(wbinc,wdc,/stairstep,xra=[-3.5,2.])
;    p = plot(xrld,yrld,/stairstep,/ov,col='dodger blue')
;    p = plot(xnon,ynon,/stairstep,/ov,col='orange')
;    p = plot(xpdf,nm(ypdf)*max(wdc),/stairstep,/ov,col='red')
;    ;; plot fine grid
;    ;samp_plt = interpol(udx,cdf,randomu(seed,100000))
;    ;yplt = histogram(samp_plt,locations=xplt,bin=width(dx),min=dx[0],max=dx[-1])
;    ;p = plot(xplt,nm(yplt)*max(wdc),/stairstep,/ov,col='red',transparency=75)
;
;    ;p = plot(xpdf,ypdf,/stairstep)
;    ;plothist,samp_pdf,bin=0.5,xhist,yhist
;    ;p = plot(xhist,yhist,/stairstep,/ov,col='red')
;endif
;
;return, samp_pdf










