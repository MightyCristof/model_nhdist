PRO xspec_model_rxz, XCM = xcm, $
                     CONV = conv, $
                     SCAT = scat, $
                     GUPTA = gupta, $
                     POWER = power, $
                     SUFFIX = suffix


;; energy ranges and string prefix
engv = ['0.5 2','2 7','2 10']
prefv = ['052','27','210']
;; redshift bins for model
zv = [0.:0.8:0.01]
if keyword_set(power) then model = 'power' else $
                           model = 'borus'
;; cut off NH at 25, vs go to 25.5
nhlim = 1

;; prep NH array
nh = [21.:25.5:0.25]
nhlo = strtrim(10.^[21.:21.75:0.25]/1e22,2)
nhhi = string([22.:25.5:0.25],format='(f5.2)')
nnh = n_elements(nh)
nlo = n_elements(nhlo)
nhi = n_elements(nhhi)
;; values for fill
if (model eq 'power') then nhval = strtrim(10.^nh/1e22,2) else $
                           nhval = [nhlo,nhhi]
if (n_elements(suffix) eq 0) then suffix = ''

;; determine scattering fraction
if keyword_set(scat) then begin
    fs = replicate(scat,nnh)
endif else if keyword_set(gupta) then begin
    ;readcol,'Gupta+2021_Fig01_model_LR.csv',nhr,fsr,format='d,d'
    readcol,'Gupta+2021_Fig01_model_MC.csv',nhr,fsr,format='d,d'
    ;; interpolate, convert to decimal, string
    fs = string(10.^interpol(alog10(fsr),nhr,nh)/100.,format='(d6.4)')
    suffix = '_gupta'
endif

;; create XSPEC model
if keyword_set(xcm) then begin
    ;; XSPEC script setup
    openw,1,'model_'+model+'_scat'+fs[0]+suffix+'.xcm'
    printf,1,'abund angr'
    printf,1,'cosmo 70 0 0.73'
    printf,1,'method leven 10 0.01'
    printf,1,'xsect vern'
    printf,1,'xset delta 0.01'
    ;; run for chosen model
    case model of
        'borus': begin
            for i = 0,n_elements(engv)-1 do begin
                eng = engv[i]
                pref = prefv[i]
                if (pref eq '210') then zra = zv[0] else zra = zv
                printf,1,'set id [open spectra_'+pref+'_'+model+'_scat'+fs[0]+suffix+'.dat a]'
                for z = 0,n_elements(zra)-1 do begin
                    printf,1,'model  constant*phabs*zphabs(atable{/Users/ccarroll/heasoft-6.25/localmodels/borus/borus02_v170323a.fits} + zphabs*cabs*cutoffpl + constant*cutoffpl)'
                    printf,1,'              1         -1          0          0      1e+10      1e+10'
                    printf,1,'           0.01         -1          0          0     100000      1e+06'
                    printf,1,'           0.01         -1          0          0     100000      1e+06'
                    printf,1,'            0.0      -0.01     -0.999     -0.999         10         10'
                    printf,1,'            1.8       0.05        1.4       1.45       2.55        2.6'
                    printf,1,'            300         -1         20         50        500       2000'
                    printf,1,'             22        0.1         22      22.25         25       25.5'
                    printf,1,'/'
                    printf,1,'           0.05         -1'
                    printf,1,'              1         -1'
                    printf,1,'= p4'
                    printf,1,'              1       0.01          0          0      1e+20      1e+24'
                    printf,1,'              1      0.001          0          0     100000      1e+06'
                    printf,1,'= p4'
                    printf,1,'= p13'
                    printf,1,'= p5'
                    printf,1,'= p6'
                    printf,1,'= p12'
                    printf,1,'           '+fs[0]+'       '+fs[0]+'          0          0      1e+10      1e+10'
                    printf,1,'= p5'
                    printf,1,'= p6'
                    printf,1,'= p12'
                    printf,1,'new 4 '+string(zra[z],format='(f0.2)')
                    printf,1,'puts $id "[tcloutr param 4]"'
                    for n = 0,nlo-1 do begin
                        printf,1,'new 13 '+nhval[n]
                        printf,1,'new 19 '+fs[n]
                        printf,1,'flux '+eng
                        printf,1,'puts $id "[tcloutr flux 2]"'
                    endfor
                    ;; link obscuration to torus at log NH = 22
                    printf,1,'new 13 = 1.0e-22*10.0^p7'
                    for n = nlo,nnh-1 do begin
                        printf,1,'new 7 '+nhval[n]
                        printf,1,'new 19 '+fs[n]
                        printf,1,'flux '+eng
                        printf,1,'puts $id "[tcloutr flux 2]"'
                    endfor                
                endfor
                printf,1,'close $id'
            endfor
            end
        'power': begin
            for i = 0,n_elements(engv)-1 do begin
                eng = engv[i]
                pref = prefv[i]
                if (pref eq '210') then zra = zv[0] else zra = zv
                printf,1,'set id [open spectra_'+pref+'_'+model+'_scat'+fs[0]+suffix+'.dat a]'
                for z = 0,n_elements(zra)-1 do begin
                    printf,1,'model  constant*phabs*zphabs(zphabs*cabs*cutoffpl + constant*cutoffpl)'
                    printf,1,'              1         -1          0          0      1e+10      1e+10'
                    printf,1,'           0.01         -1          0          0     100000      1e+06'
                    printf,1,'           0.01         -1          0          0     100000      1e+06'
                    printf,1,'            0.0      -0.01     -0.999     -0.999         10         10'
                    printf,1,'           0.01         -1          0          0     100000      1e+06'
                    printf,1,'= p4'
                    printf,1,'= p5'
                    printf,1,'            1.8       0.05        1.4       1.45       2.55        2.6'
                    printf,1,'            300         -1         20         50        500       2000'
                    printf,1,'              1       0.01          0          0      1e+20      1e+24'
                    printf,1,'           '+fs[0]+'       '+fs[0]+'          0          0      1e+10      1e+10'
                    printf,1,'= p8'
                    printf,1,'= p9'
                    printf,1,'= p10'
                    printf,1,'new 4 '+string(zra[z],format='(f0.2)')
                    printf,1,'puts $id "[tcloutr param 4]"'
                    for n = 0,nnh-1 do begin
                        printf,1,'new 5 '+nhval[n]
                        printf,1,'new 11 '+fs[n]
                        printf,1,'flux '+eng
                        printf,1,'puts $id "[tcloutr flux 2]"'
                    endfor
                endfor
                printf,1,'close $id'
            endfor    
            end
    endcase
    close,1
endif

;; convulate RX from model fluxes; soft and hard conversions
if keyword_set(conv) then begin
    files = file_search('spectra_*_'+model+'_scat'+fs[0]+suffix+'.dat')
    ifull = where(strmatch(files,'spectra_210*'),full0)
    ihard = where(strmatch(files,'spectra_27*'),hard0)
    isoft = where(strmatch(files,'spectra_052*'),soft0)
    if (full0 eq 0 or soft0 eq 0 or hard0 eq 0) then stop
    readcol,files[ifull],full,dfull,format='d,d'
    readcol,files[ihard],hard,dhard,format='d,d'
    readcol,files[isoft],soft,dsoft,format='d,d'
    ;; create NH array
    inh = where(dfull eq 0.,num_nh)
    ;; create restframe 2-10keV flux ratio
    full = full[inh]
    rx = alog10(full/full[0])    
    

    ;; prepare conversion arrays
    iz = where(dhard lt 0.,zlen)
    if (n_elements(zv) ne zlen) then stop
    temp_hard = dblarr(num_nh,zlen)
    temp_soft = dblarr(num_nh,zlen)
    c_hard = dblarr(num_nh,zlen)
    c_soft = dblarr(num_nh,zlen)
    for i = 0,zlen-1 do begin
        temp_hard[*,i] = hard[iz[i]+1:iz[i]+num_nh]
        temp_soft[*,i] = soft[iz[i]+1:iz[i]+num_nh]
        c_hard[*,i] = hard[iz[i]+1:iz[i]+num_nh]/full
        c_soft[*,i] = soft[iz[i]+1:iz[i]+num_nh]/full
    endfor
    hard = temp_hard
    soft = temp_soft
    ;; test reform of hard and soft fluxes (SAME)
    ;inhh = where(dhard eq 0.,nhard)
    ;inhs = where(dsoft eq 0.,nsoft)
    ;newhard = reform(hard[inhh],num_nh,zlen) 
    ;newsoft = reform(soft[inhs],num_nh,zlen) 
    ;newc_hard = newhard/rebin(full,num_nh,zlen)
    ;newc_soft = newsoft/rebin(full,num_nh,zlen)

    ;; limit NH ² 25.0
    if keyword_set(nhlim) then begin
        ilim = where(nh le 25.,ct)
        nh = nh[ilim]
        fs = fs[ilim]
        rx = rx[ilim]
        full = full[ilim]
        soft = soft[ilim,*]
        hard = hard[ilim,*]
        c_hard = c_hard[ilim,*]
        c_soft = c_soft[ilim,*]
    endif
    ;; interpolate to increase sample density
    rx_fine = [rx[0]:rx[-1]:-0.001]
    c_hard_fine = dblarr(n_elements(rx_fine),zlen)
    c_soft_fine = dblarr(n_elements(rx_fine),zlen)
    for i = 0,zlen-1 do begin
        c_hard_fine[*,i] = interpol(c_hard[*,i],rx,rx_fine)
        c_soft_fine[*,i] = interpol(c_soft[*,i],rx,rx_fine)
    endfor
    save,zv,nh,fs,full,hard,soft, $
         rx,c_hard,c_soft,rx_fine,c_hard_fine,c_soft_fine,file='rxz_'+model+'_scat'+fs[0]+suffix+'.sav'
endif


END





