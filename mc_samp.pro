FUNCTION mc_samp, arr, $
                  ndraw


;; finite values only
ifin = where(finite(arr),sz)
fin = arr[ifin]
;; sort and create EDF
is = sort(fin)
sarr = fin[is]
edf = (findgen(sz))/sz
;; sample
draw = randomu(seed,ndraw)
samp = interpol(sarr,edf,draw)

return, samp


END






