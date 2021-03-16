FUNCTION hist2d_avg, arr, $
                     bin, $
                     IIDET = iidet


;; bin size must be double precision to maintain histogram bins
if (typename(bin) ne 'DOUBLE') then message, 'BIN SIZE MUST BE DOUBLE PRECISION'

;; determine decimal place of bin size
digit = 0
escape = 0
while (escape eq 0) do begin
    escape = rnd(bin,digit)/bin
    digit++
endwhile
digit--

;; bounds of array, pushed by factor x2 bin
mm = rnd(minmax(arr)+1.*[-bin,bin],digit)
sz = size(arr,/dim)
yh = dblarr(n_elements([mm[0]:mm[1]:bin]),sz[1])
yh_det = dblarr(n_elements([mm[0]:mm[1]:bin]),sz[1])
yh_non = dblarr(n_elements([mm[0]:mm[1]:bin]),sz[1])
for i = 0,sz[1]-1 do begin
    yh[*,i] = histogram(arr[*,i],locations=xh,bin=bin,min=mm[0],max=mm[1])
    if keyword_set(iidet) then begin
        yh_det[*,i] = histogram(arr[where(iidet[*,i]),i],bin=bin,min=mm[0],max=mm[1])
        yh_non[*,i] = histogram(arr[where(~iidet[*,i]),i],bin=bin,min=mm[0],max=mm[1])
    endif
endfor
;; XH is bin start, but PLOT() assumes centered abscissa values.
;; need XH to normalize histogram (i.e., Ananna+2019)
tags = 'xh:xh,xoff:xh*0.+width(xh)/2.,yh:mean(yh,dim=2),sig:stddev(yh,dim=2)'
if keyword_set(iidet) then tags = tags+',yh_det:mean(yh_det,dim=2),sig_det:stddev(yh_det,dim=2),yh_non:mean(yh_non,dim=2),sig_non:stddev(yh_non,dim=2)'
re = execute('struct = soa2aos({'+tags+'})')

return, struct


END






