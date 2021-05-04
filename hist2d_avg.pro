FUNCTION hist2d_avg, arr, $
                     binsz, $
                     IIDET = iidet


sz = size(arr,/dim)
xh = [floor(min(arr)):ceil(max(arr)):binsz]
nbins = n_elements(xh)
xhoff = xh+binsz/2.
yh = dblarr(nbins,sz[1])
for i = 0,sz[1]-1 do yh[*,i] = histogram(arr[*,i],bin=binsz,min=xh[0],max=xh[-1])
nm = total(yh[where(xh lt 24.),*],1)
nm = rebin(transpose(nm),nbins,sz[1])
yh /= nm

tags = 'xh:xh,xhoff:xh+width(xh,/med)/2.,yh:mean(yh,dim=2),sig:stddev(yh,dim=2),mad:medabsdev(yh,dim=2)'

if keyword_set(iidet) then begin
    yh_det = dblarr(n_elements(xh),sz[1])
    yh_non = dblarr(n_elements(xh),sz[1])
    for i = 0,sz[1]-1 do begin
        yh_det[*,i] = histogram(arr[where(iidet[*,i] eq 1),i],bin=binsz,min=xh[0],max=xh[-1])
        yh_non[*,i] = histogram(arr[where(iidet[*,i] eq 0),i],bin=binsz,min=xh[0],max=xh[-1])
    endfor
    yh_det /= nm
    yh_non /= nm
    tags = tags+',yh_det:mean(yh_det,dim=2),sig_det:stddev(yh_det,dim=2),mad_det:medabsdev(yh_det,dim=2),' $
                +'yh_non:mean(yh_non,dim=2),sig_non:stddev(yh_non,dim=2),mad_non:medabsdev(yh_non,dim=2)'
endif

re = execute('struct = soa2aos({'+tags+'})')
return, struct


END






;;; binsz size must be double precision to maintain histogram binszs
;if (typename(binsz) ne 'DOUBLE') then message, 'binsz SIZE MUST BE DOUBLE PRECISION'
;
;;; determine decimal place of binsz size
;digit = 0
;escape = 0
;while (escape eq 0) do begin
;    escape = rnd(binsz,digit)/binsz
;    digit++
;endwhile
;digit--
;
;;; bounds of array, pushed by factor x2 binsz
;mm = rnd(minmax(arr)+1.*[-binsz,binsz],digit)
;mm = [floor(min(arr)),ceil(max(arr))]
;sz = size(arr,/dim)
;yh = dblarr(n_elements([mm[0]:mm[1]:binsz]),sz[1])
;yh_det = dblarr(n_elements([mm[0]:mm[1]:binsz]),sz[1])
;yh_non = dblarr(n_elements([mm[0]:mm[1]:binsz]),sz[1])
;for i = 0,sz[1]-1 do begin
;    yh[*,i] = histogram(arr[*,i],locations=xh,binsz=binsz,min=mm[0],max=mm[1])
;    if keyword_set(iidet) then begin
;        yh_det[*,i] = histogram(arr[where(iidet[*,i]),i],binsz=binsz,min=mm[0],max=mm[1])
;        yh_non[*,i] = histogram(arr[where(~iidet[*,i]),i],binsz=binsz,min=mm[0],max=mm[1])
;    endif
;endfor
;;; XH is binsz start, but PLOT() assumes centered abscissa values.
;;; need XH to normalize histogram (i.e., Ananna+2019)
;tags = 'xh:xh,xhoff:xh+width(xh,/med)/2.,yh:mean(yh,dim=2),sig:stddev(yh,dim=2)'
;if keyword_set(iidet) then tags = tags+',yh_det:mean(yh_det,dim=2),sig_det:stddev(yh_det,dim=2),yh_non:mean(yh_non,dim=2),sig_non:stddev(yh_non,dim=2)'
;re = execute('struct = soa2aos({'+tags+'})')
;
;return, struct

