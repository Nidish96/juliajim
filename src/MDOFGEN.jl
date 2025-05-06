using Base: nothing_sentinel
using LinearAlgebra
include("./HARMONIC.jl")

struct NONLINEARITY
    type::Symbol	# Type of Nonlinearity, either
    			# :Inst for instantaneous nl, or
    			# :Hyst for hysteretic nl.

    func::Function	# Function handle to nonlinearity

    L::Matrix{Float64}  # Selection matrix
    Lf			# Force shape matrix (can be empty)

    NONLINEARITY(type, func, L) = new(type, func, L, []);
    NONLINEARITY(type, func, L, Lf) = new(type, func, L, Lf);
end


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

function ADDNL(m::MDOFGEN, type::Symbol, func::Function, L, Lf=nothing)
    push!(m.NLTs, NONLINEARITY(type, func, L, Lf));

    m = MDOFGEN(m, m.NLTs);
end

function NLEVAL!(FNL::Union{Vector, Matrix, SparseMatrixCSC, Nothing},
                 dFNLdU::Union{Matrix, SparseMatrixCSC, Nothing},
                 dFNLdw::Union{Vector, Matrix, SparseMatrixCSC, Nothing},
                 m::MDOFGEN,
                 Uw::Union{Vector, Matrix, SparseMatrixCSC, Nothing}, h::Matrix{Int64},
                 N::Int64, tol::Float64=nothing)
    Nhc = sum(all(h.==0, dims=2) + 2*any(h.!=0, dims=2));
    t = (0:N-1)*2π/N;
    w = Uw[end];
    
    D1 = zeros(Nhc, Nhc);
    HARMONICSTIFFNESS!(D1, nothing, 0, 1.0, 0, [w], h);

    cst = reshape(AFT(I(Nhc), h, N, :f2t), N,Nhc,1);
    sct = reshape(AFT(D1, h, N, :f2t), N,Nhc,1);

    if !(FNL === nothing)
        FNL[:] = zeros(m.Ndofs*Nhc, 1);
    end
    if !(dFNLdU === nothing)
        dFNLdU[:] = zeros(m.Ndofs*Nhc, m.Ndofs*Nhc);
    end
    if !(dFNLdw === nothing)    
        dFNLdw[:] = zeros(m.Ndofs*Nhc, 1);
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
                J = zeros(Ndnl * Nhc, Ndnl * Nhc);
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
        elseif m.NLTs[ni].type==:Dlag
            error("Dynamic Lagrangian needs to be implemented");
        else
            error("Unknown type at the ", ni, "th Nonlinearity.");
        end

        if m.NLTs[ni].Lf === nothing  # Self adjoint forcing
            if !(FNL === nothing)
                FNL[:] += reshape(m.NLTs[ni].L' * F', Nhc * MDL.Ndofs, 1);
            end
            if !(dFNLdU === nothing)
                dFNLdU[:,:] += kron(I(Nhc), m.NLTs[ni].L') * J * kron(I(Nhc), m.NLTs[ni].L);
            end
            if !(dFNLdw === nothing)
                dFNLdw[:] += reshape(m.NLTs[ni].L' * dFdw', Nhc * MDL.Ndofs, 1);
            end
        else  # Non-self adjoint forcing
            if !(FNL === nothing)
                FNL += reshape(m.NLTs[ni].Lf * F', Nhc * MDL.Ndofs, 1);
            end
            if !(dFNLdU === nothing)
                dFNLdU += kron(I(Nhc), m.NLTs[ni].Lf) * J * kron(I(Nhc), m.NLTs[ni].L)
            end
            if !(dFNLdw === nothing)
                dFNLdw += reshape(m.NLTs[ni].Lf * dFdw', Nhc * MDL.Ndofs, 1)
            end
        end
    end
end

function HBRESFUN!(R::Union{Vector, Matrix, SparseMatrixCSC, Nothing},
                   dRdU::Union{Matrix, SparseMatrixCSC, Nothing},
                   dRdw::Union{Vector, Matrix, SparseMatrixCSC, Nothing},
                   m::MDOFGEN,
                   Uw::Union{Vector, Matrix, SparseMatrixCSC, Nothing},
                   Fl::Union{Vector, Matrix, SparseMatrixCSC, Nothing},
                   h::Matrix{Int64}, N::Int64, tol::Float64=nothing)
    w = Uw[end];

    E = spzeros(Nhc*MDL.Ndofs, Nhc*MDL.Ndofs);
    if !(dRdw === nothing)
        dEdw = [spzeros(Nhc*MDL.Ndofs, Nhc*MDL.Ndofs)];
    else
        dEdw = nothing;
    end
    HARMONICSTIFFNESS!(E, dEdw, m.M, m.C, m.K, [w], h);

    # Evaluate Nonlinearity in Frequency Domain
    NLEVAL!(R, dRdU, dRdw, m, Uw, h, N, tol);
    
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
end

function HHTAMARCH(m::MDOFGEN, T0, T1, dt, U0, Ud0, Fex, opts=nothing)
    print("To do")
end
