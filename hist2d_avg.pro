FUNCTION hist2d_avg, arr, $
                     binsz, $
                     NORMALIZE = normalize, $
                     CONFIDENCE = confidence, $
                     IIDET = iidet


sz = size(arr,/dim)
xh = [floor(min(arr)-binsz):ceil(max(arr)+binsz):binsz]
nbins = n_elements(xh)
xhoff = xh+binsz/2.
yh = dblarr(nbins,sz[1])
for i = 0,sz[1]-1 do yh[*,i] = histogram(arr[*,i],bin=binsz,min=xh[0],max=xh[-1])
if keyword_set(normalize) then begin
    if (typename(normalize) eq 'STRING') then re = execute('nm = total(yh['+normalize+',*],1)') else $
                                              nm = total(yh,1)
    nm = rebin(transpose(nm),nbins,sz[1])
    yh /= nm
endif

tags = 'xh:xh,xhoff:xh+width(xh,/med)/2.,yh:mean(yh,dim=2),mad:medabsdev(yh,dim=2),sig:stddev(yh,dim=2)'
if keyword_set(confidence) then begin
    sig_up = dblarr(nbins,3)
    sig_lo = dblarr(nbins,3)
    for i = 0,nbins-1 do begin
        edf = edf(yh[i,*],xloc=xedf)
        sig_up[i,*] = xedf[value_locate(edf,[0.68,0.95,0.99])]
        sig_lo[i,*] = xedf[value_locate(edf,[0.32,0.05,0.01])]
    endfor
    tags = tags+',sig1u:sig_up[*,0],sig2u:sig_up[*,1],sig3u:sig_up[*,2],sig1l:sig_lo[*,0],sig2l:sig_lo[*,1],sig3l:sig_lo[*,2]'
endif

if keyword_set(iidet) then begin
    yh_det = dblarr(n_elements(xh),sz[1])
    yh_non = dblarr(n_elements(xh),sz[1])
    for i = 0,sz[1]-1 do begin
        yh_det[*,i] = histogram(arr[where(iidet[*,i] eq 1),i],bin=binsz,min=xh[0],max=xh[-1])
        yh_non[*,i] = histogram(arr[where(iidet[*,i] eq 0),i],bin=binsz,min=xh[0],max=xh[-1])
    endfor
    if keyword_set(normalize) then begin
        yh_det /= nm
        yh_non /= nm
    endif
    tags = tags+',yh_det:mean(yh_det,dim=2),sig_det:stddev(yh_det,dim=2),mad_det:medabsdev(yh_det,dim=2),' $
                +'yh_non:mean(yh_non,dim=2),sig_non:stddev(yh_non,dim=2),mad_non:medabsdev(yh_non,dim=2)'

    if keyword_set(confidence) then begin    
        sig_up_det = dblarr(nbins,3)
        sig_lo_det = dblarr(nbins,3)
        sig_up_non = dblarr(nbins,3)
        sig_lo_non = dblarr(nbins,3)
        for i = 0,nbins-1 do begin
            edf_det = edf(yh_det[i,*],xloc=xedf_det)
            sig_up_det[i,*] = xedf_det[value_locate(edf_det,[0.68,0.95,0.99])]
            sig_lo_det[i,*] = xedf_det[value_locate(edf_det,[0.32,0.05,0.01])]
            edf_non = edf(yh_non[i,*],xloc=xedf_non)
            sig_up_non[i,*] = xedf_non[value_locate(edf_non,[0.68,0.95,0.99])]
            sig_lo_non[i,*] = xedf_non[value_locate(edf_non,[0.32,0.05,0.01])]
        endfor
        tags = tags+',sig1u_det:sig_up_det[*,0],sig2u_det:sig_up_det[*,1],sig3u_det:sig_up_det[*,2],' $
                    +'sig1l_det:sig_lo_det[*,0],sig2l_det:sig_lo_det[*,1],sig3l_det:sig_lo_det[*,2],' $
                    +'sig1u_non:sig_up_non[*,0],sig2u_non:sig_up_non[*,1],sig3u_non:sig_up_non[*,2],' $
                    +'sig1l_non:sig_lo_non[*,0],sig2l_non:sig_lo_non[*,1],sig3l_non:sig_lo_non[*,2]'
    endif
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

