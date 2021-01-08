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
dx = [xh[0]:xh[-1]:width(minmax(xh))/1000.]

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

return, resamp_nh


END



