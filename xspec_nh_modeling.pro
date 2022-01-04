PRO xspec_model_rxz, XCM = xcm, $
                     CONV = conv, $
                     SCAT = scat, $


model = 'borus'
gupta = 1


;; energy ranges and string prefix
engv = ['2 10'];['0.5 2','2 7','2 10']
prefv = ['210'];['052','27','210']

;; redshift bins for model
zv = [0.:0.8:0.01]

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
nhval = [nhlo,nhhi]

;; scattering fraction
readcol,'Gupta+2021_Fig1_data.csv',nhrv,fsrv,format='d,d'
;; choice of model
fsr = fsrv[0:1]  ;; exclude upper limits
nhr = nhrv[0:1]
;fsr = fsrv[2:3]  ;; MC simulations
;nhr = nhrv[2:3]
;; interpolate, convert to decimal, string
fs = string(10.^interpol(alog10(fsr),nhr,nh)/100.,format='(d6.4)')
suffix = '_gupta'



;; create XSPEC model

;; XSPEC script setup
openw,1,'model_'+model+'_scat'+fs[0]+suffix+'.xcm'
printf,1,'abund angr'
printf,1,'cosmo 70 0 0.73'
printf,1,'method leven 10 0.01'
printf,1,'xsect vern'
printf,1,'xset delta 0.01'


for i = 0,n_elements(engv)-1 do begin
    eng = engv[i]
    pref = prefv[i]



for energy
for photon index
for scattering fraction
for compton reflection





END



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
        ;;  Model Model Component  Parameter  Unit     Value
        ;;   par  comp
        ;;     1    1   constant   factor              1.00000      frozen
        ;;     2    2   phabs      nH         10^22    1.00000E-02  frozen
        ;;     3    3   zphabs     nH         10^22    1.00000E-02  frozen
        ;;     4    3   zphabs     Redshift            0.0          frozen
        ;;     5    4   borus02    PhoIndex            1.80000      +/-  0.0
        ;;     6    4   borus02    Ecut       keV      300.000      frozen
        ;;     7    4   borus02    logNHtor            22.0000      +/-  0.0
        ;;     8    4   borus02    CFtor               0.500000     +/-  0.0
        ;;     9    4   borus02    cos(thInc)          5.00000E-02  frozen
        ;;    10    4   borus02    A_Fe       A_Sun    1.00000      frozen
        ;;    11    4   borus02    z                   0.0          = p4
        ;;    12    4   borus02    norm                1.00000      +/-  0.0
        ;;    13    5   zphabs     nH         10^22    1.00000      +/-  0.0
        ;;    14    5   zphabs     Redshift            0.0          = p4
        ;;    15    6   cabs       nH         10^22    1.00000      = p13
        ;;    16    7   cutoffpl   PhoIndex            1.80000      = p5
        ;;    17    7   cutoffpl   HighECut   keV      300.000      = p6
        ;;    18    7   cutoffpl   norm                1.00000      = p12
        ;;    19    8   constant   factor              0.177200     +/-  0.0
        ;;    20    9   cutoffpl   PhoIndex            1.80000      = p5
        ;;    21    9   cutoffpl   HighECut   keV      300.000      = p6
        ;;    22    9   cutoffpl   norm                1.00000      = p12    
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

close,1
endif






