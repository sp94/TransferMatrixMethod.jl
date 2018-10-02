epr_list = 1 .+ rand(ComplexF64, 2)
mur_list = 1 .+ rand(ComplexF64, 2)

@show epr_list, mur_list

for P3_list in ([0,0], [0,1/2], [1/2,0])
    g = Geometry([0,0], epr_list, mur_list, P3_list)
    s = Source(n1, 1, pi/4, pi/4)
    c1m, c2p = tmm(g, s)
    E1x, E1y, H1x, H1y = W(g, s, 1) * vcat([s.Ex,s.Ey],c1m)
    E2x, E2y, H2x, H2y = W(g, s, 2) * vcat(c2p,[0,0])
    ΔP3 = P3_list[2] - P3_list[1]
    @test isapprox(E1x, E2x)
    @test isapprox(E1y, E2y)
    @test isapprox(H1x-2im*α*ΔP3*E1x, H2x)
    @test isapprox(H1y-2im*α*ΔP3*E1y, H2y)
end