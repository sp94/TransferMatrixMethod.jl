n1, n2 = 3, 4

g = Geometry([0,0], [n1^2,n2^2], [1,1])

function fresnel(pol,n1,n2,θ)
    θt = asin(n1/n2*sin(θ))
    @assert isapprox(n1*sin(θ), n2*sin(θt))
    if pol == "s"
        r_s = (n1*cos(θ)-n2*cos(θt)) / (n1*cos(θ)+n2*cos(θt))
        t_s = 2n1*cos(θ) / (n1*cos(θ)+n2*cos(θt))
        @assert isapprox(t_s, r_s+1)
        return r_s, t_s
    elseif pol == "p"
        r_p = (n2*cos(θ)-n1*cos(θt)) / (n2*cos(θ)+n1*cos(θt))
        t_p = 2n1*cos(θ) / (n2*cos(θ)+n1*cos(θt))
        @assert isapprox(n2/n1*t_p, r_p+1)
        return r_p, t_p
    else
        throw("Polarisation must be 's' or 'p'.")
    end
end

# p polarisation
for θ in (0:5:85)*pi/180
    res = tmm(g, 1, θ, 0)
    r, t = res.ref.p, res.trn.p
    @test isapprox(res.ref.s, 0, atol=1e-15)
    @test isapprox(res.trn.s, 0, atol=1e-15)
    r_, t_ = fresnel("p", n1, n2, θ)
    @test isapprox(r, r_)
    @test isapprox(t, t_)
end

# s polarisation
for θ in (0:5:85)*pi/180
    res = tmm(g, 1, θ, pi/2)
    r, t = res.ref.s, res.trn.s
    @test isapprox(res.ref.p, 0, atol=1e-15)
    @test isapprox(res.trn.p, 0, atol=1e-15)
    r_, t_ = fresnel("s", n1, n2, θ)
    @test isapprox(r, r_)
    @test isapprox(t, t_)
end