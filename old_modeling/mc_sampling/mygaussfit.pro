function mygaussian, p, x=x, y=y, err=err

	model = p[3] + p[0] * exp(-0.5 * ((x-p[1])/p[2])^2.)

	return,(y-model)/err

end

