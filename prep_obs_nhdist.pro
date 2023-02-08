PRO prep_obs_nhdist


common _data

file = 'data_prep/'+['Lansbury+2015_Fig13_data_cha.csv','Lansbury+2015_Fig13_data_nst.csv', $
                     'Ananna+2019_Fig10_data_lowL.csv','Ananna+2019_Fig10_data_hiL.csv', $
                     'Masoura+2021_Fig2_data.csv']
par = ['NH_LAN_CHA','NH_LAN_NST','NH_ANA_LOX','NH_ANA_HIX','NH_MAS']
vars = []

for f = 0,n_elements(file)-1 do begin
    readcol,file[f],x,y,format='d,d',delim=',',/silent
    xh = []
    yh = []

    for i = 0,n_elements(x)-2 do begin
        w = width(x[i:i+1])
        if (w eq 0) then continue
        xh = [xh,x[i]+w/2.]
        yh = [yh,y[i]]
    endfor
    w = width(xh,/med)
    xh = [xh[0]-w,xh,xh[-1]+w]
    yh = [0.,yh,0.]
    ;yh = yh/total(yh)

    re = execute(par[f]+' = soa2aos({xh:xh,yh:yh})')
    vars = [vars,par[f]]
endfor

;; add NH=20-21 based on 20% type11 sources at LX==3.6 from Merloni+2014
funob = 0.25
xh = nh_lan_nst.xh
yh = nh_lan_nst.yh
nm = rnd(total(yh[where(xh gt 21.)])/(1.-funob),1)
yh[where(xh gt 20. and xh lt 21.)] = funob*nm
;yh = yh/total(yh)
nh_lan_cor = soa2aos({xh:xh,yh:yh})
vars = [vars,'NH_LAN_COR']

;; add Ricci+2015
; # L_14-195<43.7
; nH_ricci_2015_high_lum=np.array([0.23869,0.08225,0.29148,0.38698])
; #L_14-195>43.7
; nH_ricci_2015_low_lum=np.array([0.38401,0.11192,0.20247,0.30126])
readcol,'data_prep/Ricci+2015_Fig4_data_lowL.csv',xh,yh,format='d,d',delim=',',/silent
;yh = yh/total(yh)
nh_ric_lox = soa2aos({xh:xh,yh:yh})
readcol,'data_prep/Ricci+2015_Fig4_data_hiL.csv',xh,yh,format='d,d',delim=',',/silent
;yh = yh/total(yh)
nh_ric_hix = soa2aos({xh:xh,yh:yh})
xh = nh_ric_hix.xh
yh = (nh_ric_lox.yh+nh_ric_hix.yh)/2.
;yh = yh/total(yh)
yh[where(xh lt 21.)] = 0.
nh_ric_avg = soa2aos({xh:xh,yh:yh})
vars = [vars,'NH_RIC_LOX','NH_RIC_HIX','NH_RIC_AVG']

;; add Ricci+2017 intrinsic
readcol,'data_prep/Ricci+2017_Fig23_data.csv',xh,yh,format='d,d',delim=',',/silent
nh_ric_int = soa2aos({xh:xh,yh:yh})
;; remove NH<21
yh[where(xh lt 21.)] = 0.
nh_ric_adj = soa2aos({xh:xh,yh:yh})
vars = [vars,'NH_RIC_INT','NH_RIC_ADJ']


var_str = strjoin(vars,',')
re = execute('save,'+var_str+',file="observed_nhdist.sav"')

END


