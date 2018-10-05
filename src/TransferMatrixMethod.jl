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
function W(epr, mur, k0, kx, ky, kz)
    V = -1im * Q(kx/k0, ky/k0, epr, mur) / (kz/k0)
    U = Matrix{ComplexF64}(I,size(V))
    return [+U +U; +V -V]
end

struct Geometry
    d_list::Array{Real,1}
    epr_list::Array{ComplexF64,1}
    mur_list::Array{ComplexF64,1}
    P3_list::Array{Real,1}
    function Geometry(d_list, epr_list, mur_list, P3_list)
        @assert d_list[1] == 0 && d_list[end] == 0
        @assert all((P3_list.==0) .| (P3_list.==1/2))
        return new(d_list, epr_list, mur_list, P3_list)
    end
end

function Geometry(d_list, epr_list, mur_list)
    return Geometry(d_list, epr_list, mur_list, zero(d_list))
end

struct PlaneWave # to do: put (kx,ky,kz) etc into a single array, k
    # Wave vector
    k0::Real
    kx::ComplexF64
    ky::ComplexF64
    kz::ComplexF64
    # Electric field
    Ex::ComplexF64
    Ey::ComplexF64
    Ez::ComplexF64
    # (Normalised) magnetic field
    Hx::ComplexF64
    Hy::ComplexF64
    Hz::ComplexF64
    # P and S polarisation vectors
    Px::ComplexF64
    Py::ComplexF64
    Pz::ComplexF64
    Sx::ComplexF64
    Sy::ComplexF64
    Sz::ComplexF64
    # Projections onto P and S vectors
    p::ComplexF64
    s::ComplexF64
    function PlaneWave(epr, mur, k0, kx, ky, kz, Ex, Ey, Ez)
        @assert ky == 0 # for now, until we generalise P and S polarisation vectors
        k_norm = [kx,ky,kz] 
        Hx, Hy, Hz = 1im * cross([kx,ky,kz],[Ex,Ey,Ez]) ./ (k0*mur)
        Sx, Sy, Sz = [0, 1, 0]
        Px, Py, Pz = cross([Sx,Sy,Sz],[kx,ky,kz]) ./ sqrt(epr*mur*k0^2)
        p = Ex*Px + Ey*Py + Ez*Pz
        s = Ex*Sx + Ey*Sy + Ez*Sz
        return new(k0, kx, ky, kz, Ex, Ey, Ez, Hx, Hy, Hz, Px, Py, Pz, Sx, Sy, Sz, p, s)
    end
end

function IncidentPlaneWave(epr, mur, k0, θ, ϕ)
    n = sqrt(epr*mur)
    kx = n*k0*sin(θ)
    ky = 0
    kz = n*k0*cos(θ)
    Ex =  cos(θ)cos(ϕ)
    Ey =        sin(ϕ)
    Ez = -sin(θ)cos(ϕ)
    return PlaneWave(epr, mur, k0, kx, ky, kz, Ex, Ey, Ez)
end
# Alternate definition:
# function IncidentPlaneWave(epr, mur, k0, θ, ϕ)
#     kx = sqrt(epr*mur)*k0 * sin(θ)
#     Ex = cos(θ)cos(ϕ)
#     Ey = sin(ϕ)
#     inc = TransmittedPlaneWave(epr, mur, k0, kx, 0, Ex, Ey)
#     @assert isapprox(inc.Ez, -sin(θ)cos(ϕ))
#     return inc
# end

function TransmittedPlaneWave(epr, mur, k0, kx, ky, Ex, Ey)
    kz = +sqrt(epr*mur*k0^2 - kx^2 - ky^2)
    Ez = -(kx*Ex + ky*Ey) / kz
    return PlaneWave(epr, mur, k0, kx, ky, kz, Ex, Ey, Ez)
end

function ReflectedPlaneWave(epr, mur, k0, kx, ky, Ex, Ey)
    kz = -sqrt(epr*mur*k0^2 - kx^2 - ky^2)
    Ez = -(kx*Ex + ky*Ey) / kz
    return PlaneWave(epr, mur, k0, kx, ky, kz, Ex, Ey, Ez)
end

struct Result
    g::Geometry
    inc::PlaneWave
    ref::PlaneWave
    trn::PlaneWave
    R::Real
    T::Real
    A::Real
    function Result(g, inc, ref, trn)
        @assert isapprox(norm([inc.Ex,inc.Ey,inc.Ez]), 1)
        R = norm([ref.Ex,ref.Ey,ref.Ez])^2
        T = norm([trn.Ex,trn.Ey,trn.Ez])^2
        T *= real(trn.kz/g.mur_list[end])
        T /= real(inc.kz/g.mur_list[1])
        A = 1 - R - T
        return new(g, inc, ref, trn, R, T, A)
    end
end

function T(g, inc, i)
    d = g.d_list[i]
    epr1, epr2 = g.epr_list[i:i+1]
    mur1, mur2 = g.mur_list[i:i+1]
    k1z = TransmittedPlaneWave(epr1, mur1, inc.k0, inc.kx, inc.ky, NaN, NaN).kz
    k2z = TransmittedPlaneWave(epr2, mur2, inc.k0, inc.kx, inc.ky, NaN, NaN).kz
    W1 = W(epr1, mur1, inc.k0, inc.kx, inc.ky, k1z)
    W2 = W(epr2, mur2, inc.k0, inc.kx, inc.ky, k2z)
    ϕ = exp(1im*k1z*d)
    Tϕ = diagm(0=>[ϕ, ϕ, 1/ϕ, 1/ϕ])
    ΔP3 = g.P3_list[i+1] - g.P3_list[i]
    TΔ = W2 \ diagm(0=>[1,1,1,1], -2=>[-2im*α*ΔP3,-2im*α*ΔP3]) * W1
    return TΔ * Tϕ
end

function tmm(g, k0, θ, ϕ)
    inc = IncidentPlaneWave(g.epr_list[1], g.mur_list[1], k0, θ, ϕ)
    # Build global transfer matrix
    Tg = I
    for i in 1:length(g.d_list)-1
        Tg = T(g,inc,i) * Tg
    end
    # Solve transfer matrix equation
    T11 = Tg[1:2,1:2]
    T12 = Tg[1:2,3:4]
    T21 = Tg[3:4,1:2]
    T22 = Tg[3:4,3:4]
    c1p = [inc.Ex, inc.Ey]
    c1m = -T22\T21*c1p
    c2p = T11*c1p + T12*c1m
    c2m = [0, 0]
    # Post-processing
    ref = ReflectedPlaneWave(
            g.epr_list[1], g.mur_list[1], k0, inc.kx, inc.ky, c1m...)
    trn = TransmittedPlaneWave(
            g.epr_list[end], g.mur_list[end], k0, inc.kx, inc.ky, c2p...)
    return Result(g, inc, ref, trn)
end

function get_fields(res::Result, xs::AbstractArray{<:Real,1}, zs::AbstractArray{<:Real,1})
    g = res.g
    inc = res.inc
    # Calculate the forward and backward traveling waves
    fws = Array{PlaneWave}(undef, length(g.d_list))
    bws = Array{PlaneWave}(undef, length(g.d_list))
    fws[1] = res.inc
    bws[1] = res.ref
    for i in 1:length(g.d_list)-1
        c = T(res.g,res.inc,i) * [fws[i].Ex, fws[i].Ey, bws[i].Ex, bws[i].Ey]
        epr = res.g.epr_list[i+1]
        mur = res.g.mur_list[i+1]
        fws[i+1] = TransmittedPlaneWave(res.g.epr_list[i+1], res.g.mur_list[i+1],
            res.inc.k0, res.inc.kx, res.inc.ky, c[1], c[2])
        bws[i+1] = ReflectedPlaneWave(res.g.epr_list[i+1], res.g.mur_list[i+1],
            res.inc.k0, res.inc.kx, res.inc.ky, c[3], c[4])
    end
    for w in vcat(fws, bws)
        if any(isnan.([w.Ex,w.Ey,w.Ez,w.kx,w.ky,w.kz]))
            error()
        end
    end
    # Calculate E and H
    E = zeros(ComplexF64, length(xs), length(zs), 3)
    H = zeros(ComplexF64, length(xs), length(zs), 3)
    z_list = cumsum(res.g.d_list)
    layer = 1
    @assert issorted(zs)
    for (iz,z) in enumerate(zs), (ix,x) in enumerate(xs)
        while layer < length(z_list) && z >= z_list[layer]
            layer += 1
        end
        Δx = x - xs[1]
        Δz = z - vcat(0,z_list)[layer]
        @assert layer == 1  || Δz >= 0
        @assert layer == length(res.g.d_list) || Δz <= res.g.d_list[layer]
        fw, bw = fws[layer], bws[layer]
        E[ix,iz,:] .+= [fw.Ex,fw.Ey,fw.Ez].*exp(1im*fw.kx*Δx+1im*fw.kz*Δz)
        E[ix,iz,:] .+= [bw.Ex,bw.Ey,bw.Ez].*exp(1im*bw.kx*Δx+1im*bw.kz*Δz)
        H[ix,iz,:] .+= [fw.Hx,fw.Hy,fw.Hz].*exp(1im*fw.kx*Δx+1im*fw.kz*Δz)
        H[ix,iz,:] .+= [bw.Hx,bw.Hy,bw.Hz].*exp(1im*bw.kx*Δx+1im*bw.kz*Δz)
    end
    return E, H
end

export tmm, Geometry, PlaneWave, TransmittedPlaneWave, ReflectedPlaneWave, get_fields

end # module
