PRO xspec_model_zra, XCM = xcm, $
                     GRID = grid


if keyword_set(xcm) then begin
    zra = [0.:0.8:0.01]

    openw,1,'run_model_nhdist.xcm' 
    printf,1,'abund angr'
    printf,1,'cosmo 70 0 0.73'
    printf,1,'method leven 10 0.01'
    printf,1,'xsect vern'
    printf,1,'xset delta 0.01'
    printf,1,'set id [open fx210_model_nhdist.dat a]'
    for i = 0,n_elements(zra)-1 do begin
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
        printf,1,'          0.005      0.005          0          0      1e+10      1e+10'
        printf,1,'= p5'
        printf,1,'= p6'
        printf,1,'= p12'
        printf,1,'new 4 '+string(zra[i],format='(f0.2)')
        printf,1,'puts $id "[tcloutr param 4]"'
        printf,1,'new 13 1.00000e-01'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
        printf,1,'new 13 1.77828e-01'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
        printf,1,'new 13 3.16228e-01'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
        printf,1,'new 13 5.62341e-01'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
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
        printf,1,'          0.005      0.005          0          0      1e+10      1e+10'
        printf,1,'= p5'
        printf,1,'= p6'
        printf,1,'= p12'
        printf,1,'new 4 '+string(zra[i],format='(f0.2)')
        printf,1,'new 7 22.0'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
        printf,1,'new 7 22.25'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
        printf,1,'new 7 22.5'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
        printf,1,'new 7 22.75'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
        printf,1,'new 7 23.0'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
        printf,1,'new 7 23.25'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
        printf,1,'new 7 23.5'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
        printf,1,'new 7 23.75'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
        printf,1,'new 7 24.0'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
        printf,1,'new 7 24.25'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
        printf,1,'new 7 24.5'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
        printf,1,'new 7 24.75'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
        printf,1,'new 7 25.0'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
        printf,1,'new 7 25.25'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
        printf,1,'new 7 25.5'
        printf,1,'flux 2 10'
        printf,1,'puts $id "[tcloutr flux 2]"'
    endfor
    printf,1,'close $id'
    close,1
endif


if keyword_set(grid) then begin
    readcol,'fx210_model_nhdist.dat',col,dcol,format='d,d'
    ibreak = where(dcol lt 0.,zlen)
    rx_grid = dblarr(width(ibreak)-1,zlen)
    row_z = col[ibreak]
    num_nh = width(ibreak)-1
    col_nh = [21.:21.+0.25*(num_nh-1):0.25]
    ibreak = [ibreak,ibreak[-1]+width(ibreak)]
    for i = 0,zlen-1 do begin
        readcol,'fx210_model_nhdist.dat',fx,format='d',skipline=ibreak[i]+1,numline=num_nh,/silent
        rx_grid[*,i] = alog10(fx/fx[0])
    endfor
    save,rx_grid,col_nh,row_z,file='rx_grid.sav'
endif




END





