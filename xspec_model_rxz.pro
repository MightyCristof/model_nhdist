PRO xspec_model_rxz, POWER = power, $
                     XCM = xcm, $
                     CONV = conv


scat = '0.008'
dscat = ''
if keyword_set(dscat) then conn = '_' else conn = ''
engv = ['0.5 2','2 7','2 10']
prefv = ['052','27','210']
if keyword_set(power) then model = 'power' else $
                           model = 'borus'
suffix = ''
nhlim = 1
howlowcanyougo = 0

;; prep NH array
if keyword_set(howlowcanyougo) then begin
    nhlo = strtrim(10.^[20.5:21.75:0.25]/1e22,2)
    suffix = '_nh20.5_25'
    nh = [20.5:25.5:0.25]
endif else begin
    nhlo = strtrim(10.^[21.:21.75:0.25]/1e22,2)
    nh = [21.:25.5:0.25]
endelse
nhhi = string([22.:25.5:0.25],format='(f5.2)')

;; redshift bins for model
zv = [0.:0.8:0.01]

;; create XSPEC model
if keyword_set(xcm) then begin
    if keyword_set(power) then begin
        nhval = strtrim(10.^nh/1e22,2)
        openw,1,'model_'+model+'_scat'+(strsplit(scat,'.',/extract))[-1]+conn+(strsplit(dscat,'.',/extract))[-1]+suffix+'.xcm'
        printf,1,'abund angr'
        printf,1,'cosmo 70 0 0.73'
        printf,1,'method leven 10 0.01'
        printf,1,'xsect vern'
        printf,1,'xset delta 0.01'
        for i = 0,n_elements(engv)-1 do begin
            eng = engv[i]
            pref = prefv[i]
            if (pref eq '210') then zra = zv[0] else zra = zv
            printf,1,'set id [open spectra_'+pref+'_'+model+'_scat'+(strsplit(scat,'.',/extract))[-1]+conn+(strsplit(dscat,'.',/extract))[-1]+suffix+'.dat a]'
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
                printf,1,'           '+scat+'       '+scat+'          0          0      1e+10      1e+10'
                printf,1,'= p8'
                printf,1,'= p9'
                printf,1,'= p10'
                printf,1,'new 4 '+string(zra[z],format='(f0.2)')
                printf,1,'puts $id "[tcloutr param 4]"'
                for n = 0,n_elements(nhval)-1 do begin
                    printf,1,'new 5 '+nhval[n]
                    printf,1,'flux '+eng
                    printf,1,'puts $id "[tcloutr flux 2]"'
                endfor
            endfor
            printf,1,'close $id'
        endfor    
    endif else begin
        openw,1,'model_'+model+'_scat'+(strsplit(scat,'.',/extract))[-1]+conn+(strsplit(dscat,'.',/extract))[-1]+suffix+'.xcm'
        printf,1,'abund angr'
        printf,1,'cosmo 70 0 0.73'
        printf,1,'method leven 10 0.01'
        printf,1,'xsect vern'
        printf,1,'xset delta 0.01'
        for i = 0,n_elements(engv)-1 do begin
            eng = engv[i]
            pref = prefv[i]
            if (pref eq '210') then zra = zv[0] else zra = zv
            printf,1,'set id [open spectra_'+pref+'_'+model+'_scat'+(strsplit(scat,'.',/extract))[-1]+conn+(strsplit(dscat,'.',/extract))[-1]+suffix+'.dat a]'
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
                printf,1,'           '+scat+'       '+scat+'          0          0      1e+10      1e+10'
                printf,1,'= p5'
                printf,1,'= p6'
                printf,1,'= p12'
                printf,1,'new 4 '+string(zra[z],format='(f0.2)')
                printf,1,'puts $id "[tcloutr param 4]"'
                for n = 0,n_elements(nhlo)-1 do begin
                    printf,1,'new 13 '+nhlo[n]
                    printf,1,'flux '+eng
                    printf,1,'puts $id "[tcloutr flux 2]"'
                endfor
                printf,1,'new 13 = 1.0e-22*10.0^p7'
                for n = 0,n_elements(nhhi)-1 do begin
                    if keyword_set(dscat) then begin
                        if (nhhi[n] eq '23.00') then printf,1,'new 19 '+dscat
                    endif
                    printf,1,'new 7 '+nhhi[n]
                    printf,1,'flux '+eng
                    printf,1,'puts $id "[tcloutr flux 2]"'
                endfor
            endfor
            printf,1,'close $id'
        endfor
    endelse
    close,1
endif

;; convulate RX from model fluxes; soft and hard conversions
if keyword_set(conv) then begin
    files = file_search('spectra_*_'+model+'_scat'+(strsplit(scat,'.',/extract))[-1]+conn+(strsplit(dscat,'.',/extract))[-1]+suffix+'.dat')
    ifull = where(strmatch(files,'*210*'),full0)
    ihard = where(strmatch(files,'*27*'),hard0)
    isoft = where(strmatch(files,'*052*'),soft0)
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
    c_hard = dblarr(num_nh,zlen)
    c_soft = dblarr(num_nh,zlen)
    for i = 0,zlen-1 do begin
        c_hard[*,i] = hard[iz[i]+1:iz[i]+num_nh]/full
        c_soft[*,i] = soft[iz[i]+1:iz[i]+num_nh]/full
    endfor
    ;; limit NH ² 25.0
    if keyword_set(nhlim) then begin
        ilim = where(nh le 25.,ct)
        nh = nh[ilim]
        rx = rx[ilim]
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
    save,zv,nh,rx,c_hard,c_soft,rx_fine,c_hard_fine,c_soft_fine,file='rxz_'+model+'_scat'+(strsplit(scat,'.',/extract))[-1]+conn+(strsplit(dscat,'.',/extract))[-1]+suffix+'.sav'
endif


END





