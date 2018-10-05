# Reproduction of figure 7 from "Long-range surface polaritons in ultra-thin films of silicon", Giannini et al, Optics Express 2008

using PyPlot
using TransferMatrixMethod

k0 = 2pi / 375 # 1/nanometers
d_list = [0, 345, 13, 0] # nanometers
epr_list = [1.67^2, 1.48^2, 10.8+19.9im, 1.48^2]
mur_list = ones(length(epr_list))
g = Geometry(d_list, epr_list, mur_list)
θs = (60:0.1:68)*pi/180
Rp = [tmm(g,k0,θ,0).R for θ in θs]
Rs = [tmm(g,k0,θ,pi/2).R for θ in θs]

plot(θs, Rp, label="p")
plot(θs, Rs, label="s")
ylim(0, 1)
xlabel(L"θ"*" [degrees]")
ylabel("Reflectivity")
legend()
