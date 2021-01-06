PRO prep_nhdist_obs


file = 'data_prep/'+['Lansbury+2015_Fig13_data_cha.csv','Lansbury+2015_Fig13_data_nst.csv','Ananna+2019_Fig10_data_lowL.csv','Ananna+2019_Fig10_data_hiL.csv']
par = ['nh_lan_cha','nh_lan_nst','nh_ana_lo','nh_ana_hi']
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

var_str = strjoin(vars,',')
re = execute('save,'+var_str+',file="nhdist_obs.sav"')

END


