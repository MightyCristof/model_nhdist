FUNCTION dblgaussian, p, x=x, y=y, err=err
	; p = 1st[a0, a1, a2, a3]], 2nd[a4, a5, a6, a7]
	model = (p[3] + p[0] * exp(-0.5 * ((x-p[1])/p[2])^2.)) + (p[7] + p[4] * exp(-0.5 * ((x-p[5])/p[6])^2.))
	return, (y-model)/err
end


