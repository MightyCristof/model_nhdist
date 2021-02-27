FUNCTION hist2d_avg, arr, $
                     ndraw, $
                     ALT = alt


irand = randomi(ndraw,n_elements(arr))
avg_arr = arr[irand]

if keyword_set(alt) then begin
    sz = size(arr,/dim)
    sort_arr = dblarr(sz)
    num = n_elements(arr[0,*])
    for i = 0,num-1 do sort_arr[*,i] = arr[sort(arr[*,i]),i]
    avg_arr = mean(sort_arr,dim=2)
endif

return, avg_arr


END

