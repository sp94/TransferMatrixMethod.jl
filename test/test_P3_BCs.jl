α = TransferMatrixMethod.α

k0 = 2pi/1
d_list = vcat(0, rand(10), 0)
epr_list = 1 .+ rand(ComplexF64, length(d_list))
mur_list = 1 .+ rand(ComplexF64, length(d_list))
P3_list = rand([0,1/2], length(d_list))

g = Geometry(d_list, epr_list, mur_list, P3_list)
res = tmm(g, k0, pi/3*rand(), pi/2*rand())

z_list = cumsum(d_list)
for (layer,z) in enumerate(z_list[1:end-1])
    xs = -5:5
    zs = [z-1e-15, z+1e-15]
    E, H = get_fields(res, xs, zs)
    E1x = E[:,1,1]
    E1y = E[:,1,2]
    E1z = E[:,1,3]
    H1x = H[:,1,1]
    H1y = H[:,1,2]
    H1z = H[:,1,3]
    E2x = E[:,2,1]
    E2y = E[:,2,2]
    E2z = E[:,2,3]
    H2x = H[:,2,1]
    H2y = H[:,2,2]
    H2z = H[:,2,3]
    ΔP3 = P3_list[layer+1] - P3_list[layer]
    @test isapprox(E1x, E2x)
    @test isapprox(E1y, E2y)
    @test isapprox(epr_list[layer].*E1z.-2im*α*ΔP3*mur_list[layer].*H1z, epr_list[layer+1].*E2z) # this BC was not programmed into the TMM at all :)
    @test isapprox(H1x-2im*α*ΔP3*E1x, H2x)
    @test isapprox(H1y-2im*α*ΔP3*E1y, H2y)
    @test isapprox(mur_list[layer]*H1z, mur_list[layer+1]*H2z)
end
