PRO set_nhobs, NHOBS = nhobs
                    

common _data
common _nhdist

nhobs = strupcase(nhobs)
;; set observed NH distribution for sampling
case nhobs of 
    'NH_LAN_NST': nh_obs = nh_lan_nst
    'NH_LAN_CHA': nh_obs = nh_lan_cha
    'NH_LAN_COR': nh_obs = nh_lan_cor
    'NH_RIC_AVG': nh_obs = nh_ric_avg
    'NH_RIC_ADJ': nh_obs = nh_ric_adj
    'NH_RIC_INT': nh_obs = nh_ric_int
    else: message, 'INVALID NH_OBS.'
endcase

sav_vars = ['NH_OBS']
sav_inds = []

sav_str = strjoin([sav_vars,sav_inds],',')
re = execute('save,'+sav_str+',/compress,file="select_nhobs.sav"')


END

