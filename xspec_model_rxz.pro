PRO xspec_model_rxz, SCAT = scat, $
                     XCM = xcm, $
                     CALC = calc, $
                     NHLIM = nhlim, $
                     HOWLOWCANYOUGO = howlowcanyougo


engv = ['0.5 2','2 7','2 10']
prefv = ['052','27','210']
model = 'rxz'
zv = [0.:0.8:0.01]
suffix = ''

;; create XSPEC model
if keyword_set(howlowcanyougo) then begin
    nhlo = strtrim(10.^[20.:21.75:0.25]/1e22,2)
    suffix = '_nh20'
endif else nhlo = strtrim(10.^[21.:21.75:0.25]/1e22,2)
nhhi = string([22.:25.5:0.25],format='(f5.2)')
if keyword_set(xcm) then begin
    openw,1,'model_'+model+'_scat'+(strsplit(scat,'.',/extract))[-1]+suffix+'.xcm'
    printf,1,'abund angr'
    printf,1,'cosmo 70 0 0.73'
    printf,1,'method leven 10 0.01'
    printf,1,'xsect vern'
    printf,1,'xset delta 0.01'
    for i = 0,n_elements(engv)-1 do begin
        eng = engv[i]
        pref = prefv[i]
        if (pref eq '210') then zra = zv[0] else zra = zv
        printf,1,'set id [open spectra_'+pref+'_'+model+'_scat'+(strsplit(scat,'.',/extract))[-1]+suffix+'.dat a]'
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
            printf,1,'= 1.0e-22*10.0^p7'
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
            for n = 0,n_elements(nhhi)-1 do begin
                printf,1,'new 7 '+nhhi[n]
                printf,1,'flux '+eng
                printf,1,'puts $id "[tcloutr flux 2]"'
            endfor
        endfor
        printf,1,'close $id'
    endfor
    close,1
endif

;; calculate RX from model fluxes; soft and hard conversions
if keyword_set(calc) then begin
    files = file_search('spectra_*_'+model+'_scat'+(strsplit(scat,'.',/extract))[-1]+suffix+'.dat')
    ifull = where(strmatch(files,'*210*'),full0)
    ihard = where(strmatch(files,'*27*'),hard0)
    isoft = where(strmatch(files,'*052*'),soft0)
    if (full0 eq 0 or soft0 eq 0 or hard0 eq 0) then stop
    readcol,files[ifull],full,dfull,format='d,d'
    readcol,files[ihard],hard,dhard,format='d,d'
    readcol,files[isoft],soft,dsoft,format='d,d'
    ;; create NH array
    inh = where(dfull eq 0.,num_nh)
    if keyword_set(howlowcanyougo) then nh = [20.:21.+0.25*(num_nh-1):0.25] else $
                                        nh = [21.:21.+0.25*(num_nh-1):0.25]
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
    c_hard = alog10(c_hard)
    c_soft = alog10(c_soft)
    
    save,zv,nh,rx,c_hard,c_soft,file=model+'_scat'+(strsplit(scat,'.',/extract))[-1]+suffix+'.sav'
endif


END





