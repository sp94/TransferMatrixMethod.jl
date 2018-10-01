module TransferMatrixMethod

using LinearAlgebra

α = 1/137

P(kx_,ky_,epr,mur) = [
    kx_*ky_          epr*mur-kx_*kx_;
    ky_*ky_-epr*mur         -kx_*ky_
] / epr
Q(kx_,ky_,epr,mur) = [
    kx_*ky_          epr*mur-kx_*kx_;
    ky_*ky_-epr*mur         -kx_*ky_
] / mur
function W(g, s, i)
    epr = g.epr_list[i]
    mur = g.mur_list[i]
    kz = sqrt(epr*mur*s.k0^2 - s.kx^2 - s.ky^2)
    V = -1im * Q(s.kx/s.k0, s.ky/s.k0, g.epr_list[i], g.mur_list[i]) / (kz/s.k0)
    return [Matrix{ComplexF64}(I,size(V)) Matrix{ComplexF64}(I,size(V)); V -V]
end

struct Geometry
    d_list::Array{Real,1}
    epr_list::Array{ComplexF64,1}
    mur_list::Array{ComplexF64,1}
end

struct Source
    k0::Real
    kx::Real
    ky::Real
    Ex::ComplexF64
    Ey::ComplexF64
end

function Source(n, k0, θ, ϕ)
    kx = n*k0 * sin(θ)
    ky = 0
    Ex = cos(θ)cos(ϕ)
    Ey = sin(ϕ)
    return Source(k0, kx, ky, Ex, Ey)
end

struct Solution
end

function tmm(g, s)
    T = I
    for i in 1:length(g.d_list)-1
        d = g.d_list[i]
        epr1, epr2 = g.epr_list[i:i+1]
        mur1, mur2 = g.mur_list[i:i+1]
        W1 = W(g, s, i)
        W2 = W(g, s, i+1)
        kz = sqrt(epr1*mur1*s.k0^2 - s.kx^2 - s.ky^2)
        Tϕ = diagm(0=>[exp(+1im*kz*d), exp(+1im*kz*d), exp(-1im*kz*d), exp(-1im*kz*d)])
        #TΔ(ΔP3) = diagm(0=>[1,1,1,1], -2=>[-2im*α*ΔP3,-2im*α*ΔP3])
        T = W2 \ W1 * Tϕ * T
    end
    T11 = T[1:2,1:2]
    T12 = T[1:2,3:4]
    T21 = T[3:4,1:2]
    T22 = T[3:4,3:4]
    c1m = -T22\T21*[s.Ex,s.Ey]
    c2p = T11*[s.Ex,s.Ey] + T12*c1m
    return c1m, c2p
end

export tmm, Geometry, Source

end # module
