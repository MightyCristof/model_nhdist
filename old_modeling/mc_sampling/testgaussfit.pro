restore,'fakedata.sav'
ploterror,t,r,rerr
guessp = [900.,2.,1.,1000.]
;guessp = [900.,-1.2,1.,1000.]

parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                       limits:[0.D,0]}, 4)
parinfo[*].value = guessp
parinfo[0].limited(0) = 1
parinfo[0].limits(0) = 0.d

fa = {X:t, Y:r, ERR:rerr}

p = mpfit('mygaussian', guessp, functargs=fa)
plt = plot(t,p[3] + p[0] * exp(-0.5 * ((t-p[1])/p[2])^2.),color='red',thick=5)

print,'a0 = ',p[0]
print,'a1 = ',p[1]
print,'a2 = ',p[2]
print,'a3 = ',p[3]

stop
end