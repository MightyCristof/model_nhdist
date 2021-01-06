;-----------------------------------------------------------------------------------------
; NAME:                                                                       IDL Function
;	mcsamp
;	
; PURPOSE:
;	Randomly sample input vector with MC.
;   
; CALLING SEQUENCE:
;	samphist = mcsamp( data, nsamp, [, SEED= ])   
;	
; INPUTS:
;	data			- Input vector to be sampled.
;	nsamp			- Number of samples to be pulled from data.
;	   
; OPTIONAL INPUTS:
;   SEED			- Set the value of seed for call to RANDOMU().
;	PLT				- Keyword to show histogram of proj.
; OUTPUTS:
;   proj			- Projection of data, resampled.
;	
; OPTIONAL OUTPUTS:
;  
; COMMENTS:
;   
; EXAMPLES:
;	Resample a random normal distribution at higher resolution.
;	 
;	in_arr = randomu(seed,100)
;	out_arr = samp(in_arr,10000,/plt)
;	
; REVISION HISTORY:
;   2018-Jun-19  Written by Christopher M. Carroll (Dartmouth)
;-----------------------------------------------------------------------------------------
FUNCTION samp, data, $
               nsamp, $
               SEED = seed, $
               PLT = plt
                 
                 
;; histogram the input data
binsz = scott(data)
yhist = histogram(data,bin=binsz,locations=xhist)
;; add endpoints to the histogram
yhist = [0.,yhist,0.]
xhist = [xhist[0]-binsz,xhist,xhist[-1]+binsz]

;; increase resolution of the data by nsamp
xq = [xhist[0]:xhist[-1]:binsz/nsamp]
;; probability distribution function
pdf = spline(xhist,yhist,xq)
pdf >= 0.
pdf /= total(pdf)
;; cumulative distribution function
cdf = total(pdf,/cumulative)
;; uniq the CDF
iu = uniq(cdf)
x_cdf = xq[iu]
y_cdf = cdf[iu]

;; backwards re-sample the CDF
cdf_ysamp = randomu(seed,nsamp)
samp_pdf = interpol(x_cdf,y_cdf,cdf_ysamp)
;; if set, display the input and sampled histograms
if keyword_set(plt) then begin
	y_samp = histogram(samp_pdf,bin=scott(cdf_xsamp),locations=x_samp)
	p = plot(xhist,nm(yhist),/stairstep,yr=[0.,1.05])           ;; original data
	p = plot(x_samp,nm(y_samp),'red',/ov)                      ;; resampled data
endif

return, pdf_samp


END


