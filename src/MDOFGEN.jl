using Base: nothing_sentinel
using LinearAlgebra
include("./HARMONIC.jl")

# * Define Structs
"""
NONLINEARITY

Struct storing a nonlinearity.

# Fields
- `type::Symbol`: Type of nonlinearity. One of:
  1. :Inst (for instantaneous nonlinearity)
  2. :Hyst (for hysteretic nonlinearity)
  3. :Fdom (for nonlinearity defined in frequency domain)

     Out of these, only the first 2 can be used for transient simulations.
- `func::Function`: Function handle for evaluating nonlinearity. Has signature
  - `(t, u, udot)` for :Inst
  - `(t, u, s, up, sp)` for :Hyst, where `s` is some internal state
- `L::Matrix{Float64}`: Selection Matrix
- `Lf`: Force shape matrix (can be empty)
"""
struct NONLINEARITY
    type::Symbol	# Type of Nonlinearity, either
    			# :Inst for instantaneous nl, or
                        # :Hyst for hysteretic nl, or
                        # :Fdom for frequency domain definition.

    func::Function	# Function handle to nonlinearity

    L::Matrix{Float64}  # Selection matrix
    Lf			# Force shape matrix (can be empty)

    NONLINEARITY(type, func, L) = new(type, func, L, []);
    NONLINEARITY(type, func, L, Lf) = new(type, func, L, Lf);
end

"""
MDOFGEN

Struct storing the generic MDOF system.

# Fields
- `M::Matrix{Float64}`: Mass matrix
- `C::Matrix{Float64}`: Damping Matrix
- `K::Matrix{Float64}`: Stiffness Matrix
- `L`: Displacement Transform Matrix (can be empty)
- `Ndofs::Int64`: Number of DOFs
- `NLTs::Vector{NONLINEARITY}`: Vector of nonlinearities present (see `NONLINEARITY`)
"""
struct MDOFGEN
    M::Matrix{Float64}  # Mass matrix
    C::Matrix{Float64}  # Damping Matrix
    K::Matrix{Float64}  # Stiffness Matrix
    
    L			# Displacement Transform Matrix (can be empty)

    Ndofs::Int64  # Number of DOFs

    NLTs::Vector{NONLINEARITY}
    
    # Constructors
    MDOFGEN(M, C, K) = new(M, C, K, [], size(M,1), []);
    MDOFGEN(M, C, K, L) = new(M, C, K, L, size(M,1), []);
    MDOFGEN(MDL::MDOFGEN, NLTs) = new(MDL.M, MDL.C, MDL.K, MDL.L, MDL.Ndofs, NLTs);
end

# ** Routine to Add a nonlinearity to an MDOFGEN object
function ADDNL(m::MDOFGEN, type::Symbol, func::Function, L, Lf=nothing)
    push!(m.NLTs, NONLINEARITY(type, func, L, Lf));

    m = MDOFGEN(m, m.NLTs);
end

# * Define Routines
"""
# Description

# Arguments
- `Uw`     : 
- `m`      :
- `h`      :
- `N`      :
  [optional]
- `tol`    :   
- `FNL`    : 
- `dFNLdU` : 
- `dFNLdw` : 
"""
function NLEVAL!(Uw::Union{Vector, Matrix, SparseMatrixCSC, Nothing},
    m::MDOFGEN,
    h,
    N::Int64; tol::Float64=eps()^(4//5),
    FNL=nothing,
    dFNLdU=nothing,
    dFNLdw=nothing)
    Nhc = sum(all(h.==0, dims=2) + 2any(h.!=0, dims=2));
    t = (0:N-1)*2π/N;
    w = Uw[end];
    
    D1 = zeros(eltype(Uw), Nhc, Nhc);
    HARMONICSTIFFNESS!(D1, nothing, 0, 1.0, 0, [w], h);

    cst = reshape(AFT(float(I(Nhc)), h, N, :f2t), N,Nhc,1);
    sct = reshape(AFT(D1, h, N, :f2t), N,Nhc,1);

    if !(FNL === nothing)
        FNL[:] = zeros(eltype(Uw), m.Ndofs*Nhc);
    end
    if !(dFNLdU === nothing)
        dFNLdU[:] = zeros(eltype(Uw), m.Ndofs*Nhc, m.Ndofs*Nhc);
    end
    if !(dFNLdw === nothing)    
        dFNLdw[:] = zeros(eltype(Uw), m.Ndofs*Nhc);
    end
    for ni in 1:length(m.NLTs)
        Ndnl = size(m.NLTs[ni].L, 1);
        Unl = (m.NLTs[ni].L*reshape(Uw[1:end-1], m.Ndofs,Nhc))';  # Nhc, Ndnl

        unlt = AFT(Unl, h, N, :f2t);
        unldot = AFT(D1*Unl, h, N, :f2t);
        if m.NLTs[ni].type==:Inst
            # Instantaneous Nonlinearity
            ft, dfdu, dfdud = m.NLTs[ni].func(t, unlt, unldot);
            # (Nt,Ndnl); (Nt,Ndnl); (Nt,Ndnl) (point-wise)
            if !(FNL === nothing)
                F = AFT(ft, h, N, :t2f);
            end
            if !(dFNLdU === nothing)
                J = zeros(eltype(Uw), Ndnl * Nhc, Ndnl * Nhc);
                dFdU = reshape(AFT(reshape(dfdu .* permutedims(cst, [1 3 2]) +
                    dfdud .* permutedims(sct, [1 3 2]), N, Ndnl * Nhc),
                                   h, N, :t2f), Nhc, Ndnl, Nhc)
                for di in 1:Ndnl
                    J[di:Ndnl:end, di:Ndnl:end] = dFdU[:, di, :];
                end
            end
            if !(dFNLdw === nothing)
                dFdw = AFT(dfdud .* unldot / w, h, N, :t2f);
            end
        elseif m.NLTs[ni].type==:Hyst
            error("Hysteretics need to be implemented");
        elseif m.NLTs[ni].type==:Fdom
            error("Frequency domain need to be implemented");
        elseif m.NLTs[ni].type==:Dlag
            error("Dynamic Lagrangian needs to be implemented");
        else
            error("Unknown type ", m.NLTs[ni].type, " at the ", ni, "th Nonlinearity.");
        end

        if m.NLTs[ni].Lf === nothing  # Self adjoint forcing
            if !(FNL === nothing)
                FNL[:] += reshape(m.NLTs[ni].L' * F', Nhc * m.Ndofs, 1);
            end
            if !(dFNLdU === nothing)
                dFNLdU[:,:] += kron(I(Nhc), m.NLTs[ni].L') * J * kron(I(Nhc), m.NLTs[ni].L);
            end
            if !(dFNLdw === nothing)
                dFNLdw[:] += reshape(m.NLTs[ni].L' * dFdw', Nhc * m.Ndofs, 1);
            end
        else  # Non-self adjoint forcing
            if !(FNL === nothing)
                FNL += reshape(m.NLTs[ni].Lf * F', Nhc * m.Ndofs, 1);
            end
            if !(dFNLdU === nothing)
                dFNLdU += kron(I(Nhc), m.NLTs[ni].Lf) * J * kron(I(Nhc), m.NLTs[ni].L)
            end
            if !(dFNLdw === nothing)
                dFNLdw += reshape(m.NLTs[ni].Lf * dFdw', Nhc * m.Ndofs, 1)
            end
        end
    end
end

"""
# Description

# Arguments
- `Uw`   : 
- `m`    : 
- `Fl`   : 
- `h`    : 
- `N`    : 
- `tol`  : 
- `dRdU` : 
- `dRdw` : 
"""
function HBRESFUN!(Uw,
    m::MDOFGEN,
    Fl::Union{Vector, Matrix, SparseMatrixCSC, Nothing},
    h, N::Int64; tol::Float64=eps()^(4//5),
    R=nothing,
    dRdU=nothing,
    dRdw=nothing)
    w = Uw[end];

    E = spzeros(eltype(Uw), Nhc*m.Ndofs, Nhc*m.Ndofs);
    if !(dRdw === nothing)
        dEdw = [spzeros(eltype(Uw), Nhc*m.Ndofs, Nhc*m.Ndofs)];
    else
        dEdw = nothing;
    end
    HARMONICSTIFFNESS!(E, dEdw, m.M, m.C, m.K, [w], h);

    # Evaluate Nonlinearity in Frequency Domain
    NLEVAL!(Uw, m, h, N; tol=tol, FNL=R, dFNLdU=dRdU, dFNLdw=dRdw);
    
    # Residue
    if !(R === nothing)
        R[:] = E*Uw[1:end-1]+R-Fl;
    end
    if !(dRdU === nothing)
        dRdU[:, :] = E+dRdU;
    end
    if !(dRdw === nothing)
        dRdw[:] = dEdw[1]*Uw[1:end-1]+dRdw;
    end
    return nothing;
end

function HHTAMARCH(m::MDOFGEN, T0, T1, dt, U0, Ud0, Fex, opts=nothing)
    print("To do")
end
