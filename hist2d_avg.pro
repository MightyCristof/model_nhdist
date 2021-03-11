FUNCTION hist2d_avg, arr, $
                     SRT = srt, $
                     HIST = hist, $
                     DRAW = draw


if keyword_set(srt) then begin
    sz = size(arr,/dim)
    sort_arr = dblarr(sz)
    for i = 0,sz[1]-1 do sort_arr[*,i] = arr[sort(arr[*,i]),i]
    avg_arr = mean(sort_arr,dim=2)
endif

if keyword_set(hist) then begin
    bin = hist
    sz = size(arr,/dim)
    ;; determine decimal place of bin size
    digit = 0
    escape = 0
    while (escape eq 0) do begin
        escape = rnd(bin,digit)/bin
        digit++
    endwhile
    digit--
    ;; bounds of array, pushed by factor x1 bin
    mm = rnd(minmax(arr)+2.*[-bin,bin],digit)
    yh = dblarr(n_elements([mm[0]:mm[1]:bin]),sz[1])
    for i = 0,sz[1]-1 do yh[*,i] = histogram(arr[*,i],locations=xh,bin=bin,min=mm[0],max=mm[1])
    avg_arr = mean(yh,dim=2)
    avg_arr = [[xh],[avg_arr]]
endif

if keyword_set(draw) then begin
    irand = randomi(draw,n_elements(arr))
    avg_arr = arr[irand]
endif

return, avg_arr


END

