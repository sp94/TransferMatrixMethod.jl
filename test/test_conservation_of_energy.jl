# Check there is no reflection when impedance is unchanged
g = Geometry([0, 10, 0], fill(1+rand(),3), fill(1+rand(),3))
res = tmm(g, 1, pi/3*rand(), pi/2*rand())
@test isapprox(res.R, 0, atol=1e-15)
@test isapprox(res.T, 1)
@test isapprox(res.A, 0, atol=1e-15)

# Check there is no absorbance in a lossless device
g = Geometry([0, 10, 0], [1, 2, 1], fill(1+rand(),3))
res = tmm(g, 1, pi/3*rand(), pi/2*rand())
@test isapprox(res.A, 0, atol=1e-15)

# Check there is absorbance when imag(ep) > 0
g = Geometry([0, 10, 0], [1, 2+1im, 1], fill(1+rand(),3))
res = tmm(g, 1, pi/3*rand(), pi/2*rand())
@test res.A > 0

# Check there is gain when imag(ep) < 0
g = Geometry([0, 10, 0], [1, 2-1im, 1], fill(1+rand(),3))
res = tmm(g, 1, pi/3*rand(), pi/2*rand())
@test res.A < 0
