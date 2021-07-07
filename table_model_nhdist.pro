PRO table_model_nhdist, FLUXES = fluxes    


common _data
common _nhdist
common _nhobs
common _rxnh
common _group
common _fixed
common _free
common _split
common _model


file_mkdir,'tables'

if keyword_set(fluxes) then begin
    ;; stacking estimates
    stk_non = fxstak[1,0:1]
    e_stk_non = e_fxstak[1,0:1]
    stk_det = fxstak[0,0:1]
    e_stk_det = e_fxstak[0,0:1]
    ;; flux estimates - fixed
    ;mod_non = [mode(fx_non.soft,bin=scott(fx_non.soft)),mode(fx_non.hard,kde=kde_bandwidth(fx_non.hard))]
    ;e_mod_non = sqrt([mode(fx_non.e_soft,kde=kde_bandwidth(fx_non.e_soft)),mode(fx_non.e_hard,kde=kde_bandwidth(fx_non.e_hard))]^2.+([stddev(fx_non.soft),stddev(fx_non.hard)])^2.)
    ;mod_det = [mode(fx_det.soft,bin=scott(fx_det.soft)),mode(fx_det.hard,kde=kde_bandwidth(fx_det.hard))]
    ;e_mod_det = sqrt([mode(fx_det.e_soft,kde=kde_bandwidth(fx_det.e_soft)),mode(fx_det.e_hard,kde=kde_bandwidth(fx_det.e_hard))]^2.+([stddev(fx_det.soft),stddev(fx_det.hard)])^2.)
;
    ;print, '================================================================='
    ;print, 'FIXED MODELING'
    ;print, ''
    ;print, '                              SOFT                       HARD'
    ;print, 'STACK DET:   '+strjoin(string(stk_det,format='(e8.2)')+'$\pm$'+string(e_stk_det,format='(e8.2)'),'      ')
    ;print, 'MODEL DET:   '+strjoin(string(mod_det,format='(e8.2)')+'$\pm$'+string(e_mod_det,format='(e8.2)'),'      ')
    ;print, 'STACK NON:   '+strjoin(string(stk_non,format='(e8.2)')+'$\pm$'+string(e_stk_non,format='(e8.2)'),'      ')
    ;print, 'MODEL NON:   '+strjoin(string(mod_non,format='(e8.2)')+'$\pm$'+string(e_mod_non,format='(e8.2)'),'      ')
    ;print, '================================================================='
    
    ;; FREE modeling
    mod_non = [mode(fx_non_.soft,bin=scott(fx_non_.soft)),mode(fx_non_.hard,kde=kde_bandwidth(fx_non_.hard))]
    e_mod_non = [mode(fx_non_.e_soft,kde=kde_bandwidth(fx_non_.e_soft)),mode(fx_non_.e_hard,kde=kde_bandwidth(fx_non_.e_hard))]
    ;e_mod_non = [stddev(fx_non_.soft),stddev(fx_non_.hard)]
    mod_det = [mode(fx_det_.soft,bin=scott(fx_det_.soft)),mode(fx_det_.hard,kde=kde_bandwidth(fx_det_.hard))]
    e_mod_det = [mode(fx_det_.e_soft,kde=kde_bandwidth(fx_det_.e_soft)),mode(fx_det_.e_hard,kde=kde_bandwidth(fx_det_.e_hard))]
    ;e_mod_det = [stddev(fx_det_.soft),stddev(fx_det_.hard)]
    
    print, '================================================================='
    print, 'FREE MODELING'
    print, ''
    print, '                              SOFT                       HARD'
    print, 'STACK DET:   '+strjoin(string(stk_det,format='(e8.2)')+'$\pm$'+string(e_stk_det,format='(e8.2)'),'      ')
    print, 'MODEL DET:   '+strjoin(string(mod_det,format='(e8.2)')+'$\pm$'+string(e_mod_det,format='(e8.2)'),'      ')
    print, 'STACK NON:   '+strjoin(string(stk_non,format='(e8.2)')+'$\pm$'+string(e_stk_non,format='(e8.2)'),'      ')
    print, 'MODEL NON:   '+strjoin(string(mod_non,format='(e8.2)')+'$\pm$'+string(e_mod_non,format='(e8.2)'),'      ')
    print, '================================================================='
    
    frac_det = mod_det/stk_det
    frac_non = mod_non/stk_non
    
    print, '================================================================='
    print, 'DIFFERENCE'
    print, ''
    print, '                              SOFT                     HARD                     MEAN'
    print, 'RATIO DET:   '+strjoin(strtrim(frac_det,2)+'('+strtrim(alog10(frac_det),2)+')','    ')+'    '+strtrim(mean(frac_det),2)+'('+strtrim(alog10(mean(frac_det)),2)+')'
    print, 'RATIO NON:   '+strjoin(strtrim(frac_non,2)+'('+strtrim(alog10(frac_non),2)+')','    ')+'    '+strtrim(mean(frac_non),2)+'('+strtrim(alog10(mean(frac_non)),2)+')'
    print, '================================================================='



endif


END





