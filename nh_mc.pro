;-----------------------------------------------------------------------------------------
; NAME:                                                                       IDL Function
;   nh_mc.pro
;                                                                    
; PURPOSE:
;   
; CALLING SEQUENCE:
;   
; INPUTS:
;  
; OPTIONAL INPUTS:
; 
; KEYWORDS:
;   
; 
; OUTPUTS:
;   
; PROCEDURE:
;
; RESTRICTIONS: 
; 
;
; cmc
;-----------------------------------------------------------------------------------------
FUNCTION nh_mc, nh_obs, $
                niter, $
                SEED = seed, $
                PLT = plt


;; take XY from the input NH structure
xh = nh_obs.xh
yh = nh_obs.yh
;; create finer grid for x-axis
dx = [xh[0]:xh[-1]:diff(minmax(xh))/1000.]

;; interpolate and normalize to create probability distribution function
;; remember newY=spline(X,Y,newX)
pdf = spline(xh,yh,dx)
pdf >= 0.
pdf /= total(pdf)
;; sum to create cumulative distribution function
cdf = total(pdf,/cumulative)
;; unique elements of CDF
iu = uniq(cdf)
xdx = dx[iu]
cdf = cdf[iu]
;; set seed
if keyword_set(seed) then seed = seed
;; draw random values from 0-1 (same as resamp_nhling from CDF)
randcdf = randomu(seed,niter)
;; project to resamp_nhle X input variable
;; remember newY=interpol(Y,X,newX)
resamp_nh = interpol(xdx,cdf,randcdf)
if keyword_set(plt) then begin
    hproj = histogram(resamp_nh,bin=scott(resamp_nh),locations=xproj)
    e = {xra:minmax(xh), margin:0.2}
    p = histplot(xh,nm(yh),xtitle='$Original NH Dist.$',layout=[1,3,1],_extra=e,dimension=[400,600])
    p = histplot(dx,pdf,col='blue',xtitle = 'Probability Dist.',/current,layout=[1,3,2],_extra=e)
    p = histplot(xproj,nm(hproj),col='green',xtitle = 'Resampled NH Dist.',/current,layout=[1,3,3],_extra=e)
endif 

return, resamp_nh


END



