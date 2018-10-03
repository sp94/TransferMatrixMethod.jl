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

function tmm(g, k0, θ, ϕ)
    inc = IncidentPlaneWave(g.epr_list[1], g.mur_list[1], k0, θ, ϕ)
    # Build transfer matrix
    T = I
    for i in 1:length(g.d_list)-1
        d = g.d_list[i]
        epr1, epr2 = g.epr_list[i:i+1]
        mur1, mur2 = g.mur_list[i:i+1]
        ΔP3 = g.P3_list[i+1] - g.P3_list[i]
        k1z = TransmittedPlaneWave(epr1, mur1, k0, inc.kx, inc.ky, NaN, NaN).kz
        k2z = TransmittedPlaneWave(epr2, mur2, k0, inc.kx, inc.ky, NaN, NaN).kz
        W1 = W(epr1, mur1, k0, inc.kx, inc.ky, k1z)
        W2 = W(epr2, mur2, k0, inc.kx, inc.ky, k2z)
        Tϕ = diagm(0=>[exp(+1im*k1z*d), exp(+1im*k1z*d), exp(-1im*k1z*d), exp(-1im*k1z*d)])
        TΔ = W2 \ diagm(0=>[1,1,1,1], -2=>[-2im*α*ΔP3,-2im*α*ΔP3]) * W1
        T = TΔ * Tϕ * T
    end
    # Solve transfer matrix equation
    T11 = T[1:2,1:2]
    T12 = T[1:2,3:4]
    T21 = T[3:4,1:2]
    T22 = T[3:4,3:4]
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

export tmm, Geometry, PlaneWave, TransmittedPlaneWave, ReflectedPlaneWave

# temporary, for testing BCs
# ultimately, the tmm function should return Ex,Ey,Hx,Hy
# without the need for the user to use W
export W, α

end # module
