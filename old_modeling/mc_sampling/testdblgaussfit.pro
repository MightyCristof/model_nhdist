restore,'lldata.sav'



lerr = lly*0.01
guessp = [1.,.1,stddev(lly),0.,.2,-0.5,stddev(lly),0.]

parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                       limits:[0.D,0]}, 8)
parinfo[*].value = guessp
parinfo[0].limited(0) = 1
parinfo[0].limits(0) = 0.d
parinfo[4].limited(0) = 1
parinfo[4].limits(0) = 0.d

fa = {X:llx, Y:lly, ERR:lerr}

p = mpfit('dblgaussian', guessp, functargs=fa)
plt = plot(llx,p[3] + p[0] * exp(-0.5 * ((llx-p[1])/p[2])^2.),color='red',thick=5)
plt = plot(llx,p[7] + p[4] * exp(-0.5 * ((llx-p[5])/p[6])^2.),color='blue',thick=5,/ov)
plt = plot(llx,lly,'o',/ov)

print,'a0 = ',p[0]
print,'a1 = ',p[1]
print,'a2 = ',p[2]
print,'a3 = ',p[3]
print,'a0 -2 = ',p[4]
print,'a1 -2 = ',p[5]
print,'a2 -2 = ',p[6]
print,'a3 -2 = ',p[7]

stop
end



