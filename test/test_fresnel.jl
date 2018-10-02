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
        # Note: compared to Wikipedia, we have not reversed the sign of r_p
        r_p = (n1*cos(θt)-n2*cos(θ)) / (n2*cos(θ)+n1*cos(θt))
        t_p = 2n1*cos(θ) / (n2*cos(θ)+n1*cos(θt))
        @assert isapprox(n2/n1*t_p, 1-r_p)
        return r_p, t_p
    else
        throw("Polarisation must be 's' or 'p'.")
    end
end

# p polarisation
for θ in (0:5:85)*pi/180
    s = Source(n1, 1, θ, 0)

    θt = asin(n1/n2*sin(θ))
    st = Source(n2, 1, θt, 0)

    c1m, c2p = tmm(g, s)
    # It would be better to do something like dot(c1m,[s.Ex,s.Ey])/dot([s.Ex,s.Ey],[s.Ex,s.Ey])
    # Then we could merge the tests for p and s polarisation
    # (Or, calculate r and t in the tmm function!)
    r, t = c1m[1]/s.Ex, c2p[1]/st.Ex
    r_, t_ = fresnel("p", n1, n2, θ)
    @test isapprox(r, r_)
    @test isapprox(t, t_)
end

# s polarisation
for θ in (0:5:85)*pi/180
    s = Source(n1, 1, θ, pi/2)

    θt = asin(n1/n2*sin(θ))
    st = Source(n2, 1, θt, pi/2)

    c1m, c2p = tmm(g, s)
    r, t = c1m[2]/s.Ey, c2p[2]/s.Ey
    r_, t_ = fresnel("s", n1, n2, θ)
    @test isapprox(r, r_)
    @test isapprox(t, t_)
end