α = TransferMatrixMethod.α

d_list = [0, 0]
epr_list = 1 .+ rand(ComplexF64, 2)
mur_list = 1 .+ rand(ComplexF64, 2)

for P3_list in ([0,0], [0,1/2], [1/2,0])
    g = Geometry(d_list, epr_list, mur_list, P3_list)
    res = tmm(g, k0, pi/3*rand(), pi/2*rand())
    E1x = res.inc.Ex + res.ref.Ex
    E1y = res.inc.Ey + res.ref.Ey
    E1z = res.inc.Ez + res.ref.Ez
    H1x = res.inc.Hx + res.ref.Hx
    H1y = res.inc.Hy + res.ref.Hy
    H1z = res.inc.Hz + res.ref.Hz
    E2x = res.trn.Ex
    E2y = res.trn.Ey
    E2z = res.trn.Ez
    H2x = res.trn.Hx
    H2y = res.trn.Hy
    H2z = res.trn.Hz
    ΔP3 = P3_list[2] - P3_list[1]
    @test isapprox(E1x, E2x)
    @test isapprox(E1y, E2y)
    @test isapprox(epr_list[1]*E1z-2im*α*ΔP3*mur_list[1]*H1z, epr_list[2]*E2z) # this BC was not programmed into the TMM at all :)
    @test isapprox(H1x-2im*α*ΔP3*E1x, H2x)
    @test isapprox(H1y-2im*α*ΔP3*E1y, H2y)
    @test isapprox(mur_list[1]*H1z, mur_list[2]*H2z)
end
