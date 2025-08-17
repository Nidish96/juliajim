using NonlinearSolve
using LinearAlgebra
using ProgressMeter
using Markdown
using Printf
using Infiltrator

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
  - `(t, u, up, sp)` for :Hyst, where `s` is some internal state (`sp` is `s` at `t-Δt`)
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
# ** Frequency Domain
"""
# Description
Evaluate the nonlinearities (and Jacobian) for one period and return the harmonics

# Arguments
- `Uw`         : 
- `m::MDOFGEN` :
- `h`          :
- `N::Int64`   :
  [optional]
- `tol`        :   
- `FNL`        : 
- `dFNLdU`     : 
- `dFNLdw`     : 
"""
function NLEVAL!(Uw::Union{Vector, Matrix, SparseMatrixCSC, Nothing},
    m::MDOFGEN, h, N::Int64;
    tol::Float64=eps()^(4//5), FNL=nothing, dFNLdU=nothing, dFNLdw=nothing)
    
    Nhc = sum(all(h.==0, dims=2) + 2any(h.!=0, dims=2));
    t = (0:N-1)*2π/N;
    w = Uw[end];
    eltp = eltype(Uw);
    
    D1 = zeros(eltp, Nhc, Nhc);
    HARMONICSTIFFNESS!(D1, nothing, 0, 1.0, 0, [w], h);

    cst = AFT(float(I(Nhc)), h, N, :f2t);
    sct = AFT(D1, h, N, :f2t);

    if !(FNL === nothing)
        FNL[:] = zeros(eltp, m.Ndofs*Nhc);
    end
    if !(dFNLdU === nothing)
        dFNLdU[:] = zeros(eltp, m.Ndofs*Nhc, m.Ndofs*Nhc);
    end
    if !(dFNLdw === nothing)    
        dFNLdw[:] = zeros(eltp, m.Ndofs*Nhc);
    end
    for ni in 1:length(m.NLTs)
        Ndnl = size(m.NLTs[ni].L, 1);
        Unl = (m.NLTs[ni].L*reshape(Uw[1:end-1], m.Ndofs,Nhc))';  # Nhc, Ndnl

        unlt = AFT(Unl, h, N, :f2t); # N, Ndnl
        unldot = AFT(D1*Unl, h, N, :f2t);  # N, Ndnl
        if m.NLTs[ni].type==:Inst
            # Instantaneous Nonlinearity
            ft, dfdu, dfdud = m.NLTs[ni].func(t, unlt, unldot);
            # (N,Ndnl); (N,Ndnl)/(N,Ndnl,Ndnl); (N,Ndnl)/(N,Ndnl,Ndnl) (point-wise)
            if !(FNL === nothing)
                F = AFT(ft, h, N, :t2f);  # Nhc, Ndnl
            end
            if !(dFNLdU === nothing)
                J = zeros(eltp, Ndnl * Nhc, Ndnl * Nhc);
                if ndims(dfdu)==2
                    for di in 1:Ndnl
                        J[di:Ndnl:end, di:Ndnl:end] = AFT(dfdu[:,di].*cst +
                                                          dfdud[:,di].*sct,
                            h, N, :t2f);
                    end
                else
                    for di in 1:Ndnl
                    	for dj in 1:Ndnl
                            J[di:Ndnl:end, dj:Ndnl:end] = AFT(dfdu[:,di,dj].*cst +
                                                              dfdud[:,di,dj].*sct,
                                h, N, :t2f);
                        end
                    end
                end
            end
            if !(dFNLdw === nothing)
                dFdw = AFT(dfdud .* unldot / w, h, N, :t2f);
            end
        elseif m.NLTs[ni].type==:Hyst
            K0 = m.NLTs[ni].func(0, repeat([zeros(Ndnl)], 3)...)[2];
            ft = unlt*K0'
            dfdai = reshape(repeat(K0, 1,1,Nhc), 1, Ndnl,Ndnl,Nhc).*
                    reshape(cst, N,1,1,Nhc);

            # ft = zeros(eltp, N, Ndnl);
            # dfdai = zeros(eltp, N, Ndnl, Ndnl, Nhc);

            fprev = ft[end, :];
            its = 0;
            while its==0 || maximum(abs.(fprev-ft[end,:]))>tol
                fprev = ft[end, :];
                for (ti,tim1) in zip(1:N,circshift(1:N,1))
            	    f, dfdu, dfdup, dfdfp = m.NLTs[ni].func(t[ti], unlt[ti, :],
                        unlt[tim1, :], ft[tim1, :]);

                    ft[ti, :] .= f;
                    if !(dFNLdU === nothing)
                        dfdai[ti, :, :, :] = dfdu.*reshape(cst[ti,:], 1,1,Nhc) +
                                             dfdup.*reshape(cst[tim1,:], 1,1,Nhc) +
                                             dfdfp.*reshape(dfdai[tim1,:,:,:], 1,1,Nhc)
                    end
                end

                its += 1;
                if its>10 # Don't do more than 10 iterations
                    break;
                end
            end

            if !(FNL === nothing)
                F = AFT(ft, h, N, :t2f);
            end
            if !(dFNLdU === nothing)
                J = zeros(Ndnl*Nhc, Ndnl*Nhc);
                for di in 1:Ndnl
                    for dj in 1:Ndnl
                        J[di:Ndnl:end, dj:Ndnl:end] =
                            AFT(dfdai[:, di, dj, :], h, N, :t2f);
                    end
                end
            end
            if !(dFNLdw === nothing)
                dFdw = zeros(eltp, Nhc, Ndnl);
            end
        elseif m.NLTs[ni].type==:Fdom
            F, J, dFdw = m.NLTs[ni].func(Unl, h, N);
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
                dFNLdU[:,:] += kron(I(Nhc), m.NLTs[ni].L') * J *
                               kron(I(Nhc), m.NLTs[ni].L);
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

# *** Forced Response
"""
# Description
Residue function that can be used for Harmonic Balance forced response analysis.

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
function HBRESFUN!(Uw, m::MDOFGEN,
    Fl::Union{Vector, Matrix, SparseMatrixCSC, Nothing},
    h, N::Int64;
    tol::Float64=eps()^(4//5), R=nothing, dRdU=nothing, dRdw=nothing)
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

"""
# Description
HB residue with amplitude constraint.

# Arguments
- `Ufw`           : 
- `m::MDOFGEN`    : 
- `A::Float64`    : Amplitude level (also see `atype` below)
- `Fl`            : 
- `h`             : 
- `N::Int64`      :
- `atype::Symbol` : (default `:H1`) One of `:H1`, `:RMS`. Type of amplitude measure.
- `ashape`        : (default 1) Shape to dof for applying amplitude. If Integer, interpreted as Dof. If vector, taken as shape.
- `tol::Float64`  : (default `eps()^(4//5)`)
- `R`             : (default `nothing`)
- `dRdUf`         : (default `nothing`)
- `dRdw`          : (default `nothing`)
"""
function HBRESFUN_A!(Ufw, m::MDOFGEN, A::Float64, Fl, h, N::Int64;
    atype::Symbol=:H1, ashape::Union{Int64, Vector{Float64}}=1, 
    tol::Float64=eps()^(4//5), R=nothing, dRdUf=nothing, dRdw=nothing)
    if !(R === nothing)
        Ri = @view R[1:end-1];
    else
        Ri = nothing;
    end
    if !(dRdUf===nothing)
        dRidU = @view dRdUf[1:end-1,1:end-1];
        dRdUf[1:end-1, end] = -Fl;
    else
        dRidU = nothing;
    end
    if !(dRdw===nothing)
        dRidw = @view dRdw[1:end-1];
    else
        dRidw = nothing;
    end

    HBRESFUN!(Ufw[[1:end-2; end]], m, Ufw[end-1]*Fl, h, N;
        tol=tol, R=Ri, dRdU=dRidU, dRdw = dRidw)


    _,_, zinds,rinds,iinds = HINDS(1, h);
    if isa(ashape, Integer)
        dof = ashape
        if dof>0
            ashape = zeros(m.Ndofs)
            ashape[dof] = 1.0;
        else
            ashape = ones(m.Ndofs);
        end
    end

    hweights = zeros(eltype(Ufw), Nhc);
    if (atype==:H1)
        hweights[[rinds[1],iinds[1]]] .= 1.0;
    elseif (atype==:RMS)
        hweights .+= 0.5;
        hweights[zinds] *= 2;
    else
        error("Unknown amplitude measure ", atype)
    end
    ahweights = kron(hweights, ashape);

    if !(R===nothing)
        R[end] = Ufw[1:end-2]'*(ahweights.*Ufw[1:end-2]) - A^2;
    end
    if !(dRdUf===nothing)
        dRdUf[end, 1:end-1] = 2*(ahweights.*Ufw[1:end-2])';
    end

    return nothing;
end

# *** EPMC
"""
# Description

# Arguments
- Uwxa                : 
- m::MDOFGEN          :
- Fl                  : (Nd*Nhc) Forcing vector. Zero harmonics are used as static loads. Harmonic portions are used to provide phase constraint.
- h                   : 
- N::Int64            : 
- tol::Float64        : (default eps()^(4//5))
- atype::Symbol       : (default :H1)
- ashape::Union{Int64 : 
- Vector{Float64}}    : (default 1)
- R                   : (default nothing)
- dRdUwx              : (default nothing)
- dRda                : (default nothing)
"""
function EPMCRESFUN!(Uwxa, m::MDOFGEN, Fl, h, N::Int64;
    atype::Symbol=:H1,
    tol::Float64=eps()^(4//5), R=nothing, dRdUwx=nothing, dRda=nothing)

    Nhc = sum(all(h.==0, dims=2) + 2any(h.!=0, dims=2));

    # Preamble
    eltp = eltype(Uwxa);
    Ri = nothing;
    dRidU = nothing;
    dRidw = nothing;
    if !(R===nothing)
        Ri = @view R[1:end-2]
    end
    if !(dRdUwx===nothing)
        dRidU = @view dRdUwx[1:end-2, 1:end-2];
        dRidw = @view dRdUwx[1:end-2, end-1];
    end
    if (!(dRda===nothing) && (dRdUwx===nothing))
        dRidU = zeros(eltp, m.Ndofs*Nhc, m.Ndofs*Nhc);
        dRidw = zeros(eltp, m.Ndofs*Nhc);
    end

    # Amplitude
    la = last(Uwxa);
    A = 10^la;
    dAdla = A*log(10);

    ahweights = nothing
    Asc = A*ones(Nhc);  # DoF scaling vector
    dAscdA = ones(Nhc);
    begin
        _,_, zinds,rinds,iinds = HINDS(1, h);

        hweights = zeros(eltype(Uwxa), Nhc);
        if (atype==:H1)
            hweights[[rinds[1],iinds[1]]] .= 1.0;
        elseif (atype==:RMS)
            hweights .+= 0.5;
            hweights[zinds] *= 2;
        else
            error("Unknown amplitude measure ", atype)
        end
        ahweights = kron(diagm(hweights), m.M);

        Asc[zinds] .= 1.0;
        dAscdA[zinds] .= 0.0;

        Asc = kron(Asc, ones(m.Ndofs));
        dAscdA = kron(dAscdA, ones(m.Ndofs));
    end

    xi = Uwxa[end-1];
    w = Uwxa[end-2];
    
    # Linear harmonic stiffness
    E = spzeros(eltp, Nhc*m.Ndofs, Nhc*m.Ndofs);
    if !(dRdw === nothing)
        dEdw = [spzeros(eltp, Nhc*m.Ndofs, Nhc*m.Ndofs)];
    else
        dEdw = nothing;
    end
    HARMONICSTIFFNESS!(E, dEdw, m.M, m.C-xi*m.M, m.K, w, h);
    dEdw = dEdw[1];
    dEdxi, _ = HARMONICSTIFFNESS(zeros(size(m.M)), -m.M, zeros(size(m.M)), w, h);

    # Evaluate Nonlinear Forces
    NLEVAL!([Asc.*Uwxa[1:end-3];w], m, h, N; tol=tol, FNL=Ri, dFNLdU=dRidU, dFNLdw=dRidw);

    # Load vector
    _,_,zinds,rinds,iinds = HINDS(m.Ndofs, h);
    Fstat = zeros(eltp, m.Ndofs*Nhc);
    Fstat[zinds] = Fl[zinds];
    Fl = Fl-Fstat;

    # Build Residue
    if !(R===nothing)
        R[:] = [E*(Asc.*Uwxa[1:end-3])+Ri - Fstat;
                Uwxa[1:end-3]'*ahweights*Uwxa[1:end-3]-1.;
                Fl'Uwxa[1:end-3]];
    end
    if !(dRda===nothing)
        dRda[:] = [(E+dRidU)*(dAscdA.*Uwxa[1:end-3])*dAdla; 0; 0];
    end
    if !(dRdUwx===nothing)
        dRdUwx[:,:] = [(E+dRidU)*diagm(Asc) dEdw*(Asc.*Uwxa[1:end-3])+dRidw dEdxi*(Asc.*Uwxa[1:end-3]);
                       2Uwxa[1:end-3]'ahweights 0 0; Fl' 0 0];
    end
    
    return nothing;
end

# ** Time Domain
"""
# Description

# Arguments
- m::MDOFGEN : 
- t          : 
- U          : 
- Ud         : 
- tp         : (default nothing)
- Up         : (default nothing)
- Udp        : (default nothing)
- Sp         : (default nothing)
"""
function NLFORCE(m::MDOFGEN, t, U, Ud; tp=nothing, Up=nothing, Udp=nothing, Sp=nothing)
    eltp = eltype(U);
    if Up===nothing
        Up = zeros(eltp, m.Ndofs);
    end
    if Udp===nothing
        Udp = zeros(eltp, m.Ndofs);
    end
    if Sp===nothing
        Sp = Vector{Vector{eltp}}(undef, length(m.NLTs));
        for ni in 1:length(m.NLTs)
            if (m.NLTs[ni].type==:Hyst)
                Sp[ni] = [0.];
            end
        end
    end
    S = copy(Sp);


    F = zeros(eltp, m.Ndofs)
    dFdU = zeros(eltp, m.Ndofs, m.Ndofs);
    dFdUd = zeros(eltp, m.Ndofs, m.Ndofs);

    for (ni, nl) in enumerate(m.NLTs)  # loop over nonlinearities
        if nl.type==:Inst
            f, dfdu, dfdud = nl.func(t, nl.L*U, nl.L*Ud);
        elseif nl.type==:Hyst
            f, dfdu, _, _ = nl.func(t, nl.L*U, nl.L*Up, Sp[ni]);
            S[ni] = f;
            dfdud = zeros(eltype(dfdu), size(dfdu));
        else
            error("Nonlinearity type ", nl.type, " has not been implemented for transients");
        end

        if ndims(dfdu)==1
            dfdu = diagm(dfdu);
        end
        if ndims(dfdud)==1
            dfdud = diagm(dfdud);
        end

        # Distribute the forces
        if nl.Lf === nothing
            F += nl.L'f;
            dFdU += nl.L'dfdu*nl.L;
            dFdUd += nl.L'dfdud*nl.L;
        else
            F += nl.Lf'f;
            dFdU += nl.Lf'dfdu*nl.L;
            dFdUd += nl.Lf'dfdud*nl.L;
        end
    end

    return F, dFdU, dFdUd, S;
end

"""
# Description

# Arguments
- m::MDOFGEN : 
- T0         : 
- T1         : 
- dt         : 
- U0         : 
- Ud0        : 
- Fex        : 
"""
function NEWMARKMARCH(m::MDOFGEN, T0, T1, dt, U0, Ud0, Fex;
    beta::Float64=1.0/4, gamma::Float64=1.0/2,
    verbosity::Int=1, S0=nothing)
    
    eltp = eltype(U0);
    U0 = U0[:];
    Ud0 = Ud0[:];

    # Utility matrices
    Z1 = m.M + gamma*dt*m.C + beta*dt^2*m.K;
    Z2 = m.M - (1-gamma)*dt*m.C - (0.5-beta)*dt^2*m.K;
    Z3 = dt*m.K;

    # Setup storage
    t_ = i->T0+(i-1)*dt;
    N = floor(Int, (T1-T0)/dt+1);

    U = zeros(eltp, m.Ndofs, N);
    Ud = zeros(eltp, m.Ndofs, N);
    Udd = zeros(eltp, m.Ndofs, N);
    S = Matrix{Vector{eltp}}(undef, length(m.NLTs), N);
    if !(S0===nothing)
        S[:, 1] = S0;
    end
    if !(U0===nothing)
        U[:, 1] = U0;
    end
    if !(Ud0===nothing)
        Ud[:, 1] = Ud0;
    end

    # Initialization
    Fnl, dFnldu, dFnldud, S[:,1] = NLFORCE(m, T0, U0, Ud0);
    Udd[:, 1] = m.M\(Fex(T0)-m.C*Ud0-m.K*U0-Fnl);
    # Gradients if necessary
    dUdd0 = -m.M\[m.K+dFnldu m.C+dFnldud];  # Derivative of accel wrt initial U0,Ud0

    inds = @. !isfinite(Udd[:,1]);
    if any(inds)
        dUdd0[inds, :] .= 0.0;
        Udd[inds, 1] .= 0.0;
    end

    PHI = diagm(ones(eltp, 2m.Ndofs)); # State Transition Matrix    

    dRdudd = zeros(eltp, m.Ndofs, m.Ndofs); 
    dRd0 = zeros(eltp, m.Ndofs, 2m.Ndofs);  # residue deriv. wrt prev step [U0 Ud0]
    step = 0;
    try
        @showprogress dt=1/verbosity desc="Marching in time..." for i in 2:N
            step = i;
            # Explicit Predictor
            Udd[:, i] = Udd[:, i-1];

            # Iterative Corrector
            if isempty(m.NLTs)  # Fully linear model
                Udd[:, i] = Z1\ (Z2*Udd[:, i-1] - Z3*Ud[:, i-1] +
                                 (Fex(t_(i))-Fex(t_(i-1))) );
                dUdd1 = Z1\ (Z2*dUdd0 - kron([0. 1.], Z3));
            else
                resfun! = (udd; R=nothing,dRdudd=nothing, dRd0=nothing) -> begin
                    FnlP, dFnlPdu, dFnlPdud, _ = NLFORCE(m, t_(i),
                        U[:, i-1] + dt*Ud[:,i-1] + dt^2*((.5-beta)*Udd[:,i-1]+beta*udd),
                        Ud[:, i-1] + dt*((1-gamma)*Udd[:,i-1]+gamma*udd);
                        tp=t_(i-1), Up=U[:,i-1], Udp=Ud[:,i-1], Sp=S[:,i-1]);
                    if !(R===nothing)
                        R[:] = Z1*udd - Z2*Udd[:, i-1] + Z3*Ud[:, i-1] +
                               (FnlP-Fnl) - (Fex(t_(i))-Fex(t_(i-1)));
                    end
                    if !(dRdudd===nothing)
                        dRdudd[:,:] = Z1 + (beta*dt^2)*dFnlPdu + (gamma*dt)*dFnlPdud;
                    end
                    if !(dRd0===nothing)
                        dRd0[:,:] = (-Z2 +dt^2*(.5-beta)*dFnlPdu+dt*(1-gamma)*dFnlPdud)*dUdd0 +
                                    kron([0. 1.], Z3) +
                                    kron([1. dt], dFnlPdu) + kron([0. 1.], dFnlPdud);
                    end
                end
                fun = NonlinearFunction((r,udd,p)->resfun!(udd; R=r),
                    jac=(J,udd,p)->resfun!(udd; dRdudd=J));
                prob = NonlinearProblem(fun, Udd[:, i]);
                sol = solve(prob);
                resfun!(sol.u, dRdudd=dRdudd, dRd0=dRd0);

                # Update
                Udd[:, i] = sol.u;
                dUdd1 = -dRdudd\dRd0;
            end

            # Update States
            Ud[:, i] = Ud[:, i-1] + dt*((1-gamma)*Udd[:, i-1]+gamma*Udd[:, i]);
            U[:, i]  = U[:, i-1] + dt*Ud[:, i-1] +
                       dt^2*((.5-beta)*Udd[:, i-1]+beta*Udd[:,i]);

            # Update State Transition
            PHI = (kron([1. dt;0. 1.], I(m.Ndofs)) +
                   [dt^2*((.5-beta)*dUdd0+beta*dUdd1);
                    dt*((1-gamma)*dUdd0+gamma*dUdd1)])*PHI;

            # Update Nonlinear Forces, States
            Fnl, dFnldu, dFnldud, S[:, i] = NLFORCE(m, t_(i), U[:,i], Ud[:,i],
                tp=t_(i-1), Up=U[:,i-1], Udp=U[:,i-1], Sp=S[:, i-1]);
            Udd[:, i] = m.M\ (-m.C*Ud[:,i]-m.K*U[:,i]-Fnl+Fex(t_(i)));
            dUdd0 = -m.M\[m.K+dFnldu m.C+dFnldud];

            if mod(i, 10)==0
                yield();
            end
        end
    catch ex
        if isa(ex, InterruptException)
            println("\nInterrupted at step $step")
        else
            println("Error ", ex, " encountered at step $step")
        end
        U = U[:, 1:step-1];
        Ud = Ud[:, 1:step-1];
        Udd = Udd[:, 1:step-1];
        S = S[:, 1:step-1];
    end
    return U, Ud, Udd, S, PHI;
end
