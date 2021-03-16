PRO prep_nhdist_obs


common _data

file = 'data_prep/'+['Lansbury+2015_Fig13_data_cha.csv','Lansbury+2015_Fig13_data_nst.csv', $
                     'Ananna+2019_Fig10_data_lowL.csv','Ananna+2019_Fig10_data_hiL.csv', $
                     'Masoura+2021_Fig2_data.csv']
par = ['nh_lan_cha','nh_lan_nst','nh_ana_lo','nh_ana_hi','nh_mas']
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
    w = width(xh)
    xh = [xh[0]-w,xh,xh[-1]+w]
    yh = [0.,yh,0.]

    re = execute(par[f]+' = soa2aos({xh:xh,yh:yh})')
    vars = [vars,par[f]]
endfor

;; add NH=20-21 based on 25% unobscured fraction of Merloni+2014
nm_nst = total(nh_lan_nst[where(nh_lan_nst.xh lt 24.)].yh)
nh_xh = nh_lan_nst.xh
nh_yh = nh_lan_nst.yh/nm_nst
nh_yh[where(nh_xh gt 20. and nh_xh lt 21.)] += total(nh_yh)*0.25
nh_lan = soa2aos({xh:nh_xh,yh:nh_yh})
vars = [vars,'nh_lan']

var_str = strjoin(vars,',')
re = execute('save,'+var_str+',file="nhdist_obs.sav"')

END


