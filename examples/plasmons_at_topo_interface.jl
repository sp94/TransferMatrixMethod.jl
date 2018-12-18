using TransferMatrixMethod
using PyPlot
using LinearAlgebra

# 300 micrometer source
k0 = 2pi / 300 # 1/um

# Define geometry with
## prism region (silicon)
## spacer region (eg benzocyclobutene, 200um)
## thin film of TI (bismuth selenide, 8um)
## and substrate (silicon)
g1 = Geometry([0,200,8,0], [3.5^2,2.5,-3.4+36.7im,2.5], [1,1,1,1], [0,0,1/2,0])
# and also define the same geometry without the prism
g2 = Geometry([0,8,0], [2.65,-3.4+36.7im,2.5], [1,1,1], [0,1/2,0])

# Find absorption peak
θ_list = (0:0.001:40)*pi/180
A_list = [tmm(g1,k0,θ,0).A for θ in θ_list]
θmax = θ_list[argmax(A_list)]

# Near-field H for p-polarisation
figure()
title(L"H_y/H_0")
xlabel(L"x"*" [μm]")
ylabel(L"y"*" [μm]")
res = tmm(g1, k0, θmax, 0)
xs = -300:10:300
zs = -100:10:500
E, H = get_fields(res, xs, zs)
imshow(real(H[:,:,2])', extent=[xs[1],xs[end],zs[end],zs[1]], interpolation="bicubic", vmin=-10, vmax=10);
for z in cumsum(g1.d_list)
    axhline(z, color="white")
end
colorbar()

# Reflection, transmission, and absorption vs angle of incidence
figure()
θ_list = (25:0.001:30)*pi/180
axvline(θmax/pi*180, color="lightgrey", ls="-")
plot(θ_list/pi*180, [tmm(g1,k0,θ,0).R for θ in θ_list], label="R")
plot(θ_list/pi*180, [tmm(g1,k0,θ,0).T for θ in θ_list], "--", label="T")
plot(θ_list/pi*180, [tmm(g1,k0,θ,0).A for θ in θ_list], ":", label="A")
xlim(26, 28)
xticks([26,27,28])
ylim(0,1)
legend(loc="center right", handlelength=1.5)
xlabel("θ [deg]")
ylabel("Efficiency")

# Amplitude of transmitted s-polarised wave with and without prism
figure()
θ_list = (20:0.001:30)*pi/180
axvline(θmax/pi*180, color="lightgrey")
plot(θ_list/pi*180, [abs(tmm(g1,k0,θ,0).trn.s)*1000 for θ in θ_list], label="With\nprism")
plot(θ_list/pi*180, [abs(tmm(g2,k0,θ,0).trn.s)*1000 for θ in θ_list], "--", label="Without\nprism")
xlim(25, 28)
legend(loc="upper left", handlelength=1.5)
xlabel("θ [deg]")
ylabel(L"|s_\mathrm{trn}|\times10^{3}")
ylim(0,2.2)
