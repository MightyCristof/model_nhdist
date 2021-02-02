FUNCTION mc_samp, arr, $
                  ndraw
                  
                  
edf,arr,arr_sort,arr_edf
draw_edf = randomu(seed,ndraw)
arr_samp = interpol(arr_sort,arr_edf,draw_edf)

return, arr_samp


END






