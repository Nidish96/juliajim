# * Preamble
using Base: nothing_sentinel
using FFTW
using LinearAlgebra
using SparseArrays
using ToeplitzMatrices

# * Custom Abstract Types
const MxTypes = Union{Float64, Matrix{Float64}, Matrix{Int64}, AbstractMatrix{Bool}};
const hTypes = Union{Int,VecOrMat{Int},UnitRange{Int}};

# * Fourier Routines ############################################################
# ** Alternating Frequency-Time Transform

"""
    AFT(yin, h::hTypes, N::Int64, dir::Symbol);

Routine to conduct Time-to-Frequency and Frequency-to-Time transforms
for the general Multi frequency case.
# Arguments
+ `yin`: (N^C, Ny) if time data, and (Nhc, Ny) if frequency data.
+ `h::hTypes`: (Nh, C) list of harmonics.
+ `N::Int64`: Number of time samples for AFT.
+ `dir::Symbol`: :t2f for time-to-frequency, and :f2t for frequency-to-time.
# Examples
+ Single Frequency Case:
```julia-repl
    N = 8;
    t = (0:N-1)*2π/N;
    yt = cos.(t) + 3sin.(2t) .+ 4;
    yt = [yt yt];

    h = collect(0:3)[:,:];
    Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));

    yf = AFT(yt, h, N, :t2f);
    YT = AFT(yf, h, N, :f2t);
```
+ 2-Frequency Case:
```julia-repl
    t1 = repeat(t, 1, N);
    t2 = repeat(t', N, 1);

    yt = cos.(t1) + 3sin.(t2) .+ 4;
    h = [0 0;1 0;0 1;1 1];
    Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));

    yf = AFT([yt[:] yt[:]], h, N, :t2f);
    YT = AFT(yf, h, N, :f2t);
```
"""
function AFT(yin::VecOrMat{Float64}, h::hTypes, N::Int, dir::Symbol)
    C = size(h, 2);
    Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));
    Ny = size(yin, 2);

    Nf = Int(N / 2) + 1;
    tN = Tuple([fill(N, C); Ny]);
    fN = Tuple([Nf; fill(N, C - 1); Ny]);
    if cmp(dir, :t2f) == 0
        yf = rfft(reshape(yin, tN), 1:C) ./ (N .^ C / 2);

        yout = zeros(Nhc, Ny);
        i0 = 0;
        if all(h[1, :] .== 0)
            yout[1, :] = yf[Tuple([fill(1, C); :])...] ./ 2;
            i0 = 1;
        end
        for i in (i0+1):size(h, 1)
            inds = Tuple([mod.(h[i, :], N) .+ 1; :]);
            yout[i0+(i-i0-1)*2+1, :] = real(yf[inds...]);
            yout[i0+(i-i0-1)*2+2, :] = -imag(yf[inds...]);
        end
    elseif cmp(dir, :f2t) == 0
        yf = zeros(ComplexF64, fN);
        i0 = 0;
        if all(h[1, :] .== 0)
            yf[Tuple([fill(1, C); :])...] = real(yin[1, :]) .* 2;
            i0 = 1;
        end
        for i in (i0+1):size(h, 1)
            inds1 = Tuple([mod.(h[i, :], N).+1; :]);
            yf[inds1...] = yin[i0+(i-i0-1)*2+1, :] - im .* yin[i0+(i-i0-1)*2+2, :];

            inds2 = [mod.(N .- h[i, :] .+ i0 .- 1, N) .+ 1; :];
            inds2[1] = inds1[1];
            inds2 = Tuple(inds2);
            if inds1 != inds2
                yf[inds2...] = yin[i0+(i-i0-1)*2+1, :] + im * yin[i0+(i-i0-1)*2+2, :];
            end
        end
        yout = reshape(irfft(yf .* (N^C / 2), N, 1:C), (N^C, Ny));
    else
        error("Unknown dir for AFT");
    end
    return yout
end

# **  Harmonic Selector

"""
    HSEL(Nhmax::Int64, ws::Any=1, hcr::Int=1)

Selection of Harmonic indices based on different criteria. First C+1 rows are
[zeros(1,C); I(C)].

Implemented criteria are,
+ `hcr=1`: Σᵢ |hᵢ| < Nhmax & ∑ᵢ h_i ≥ 0

# Arguments
+ `Nhmax::Int64`: Max harmonic for truncation.
+ `ws=1`: [opt] List of frequencies
+ `hcr::Int=1`: [opt] Truncation criterion. See above.
# Examples
```julia-repl
C = 2;
Nhmax = 3;

h = HSEL(Nhmax, 1:C)
```
"""
function HSEL(Nhmax::Int64,
              ws::Union{Nothing, Float64, Vector{Float64}, UnitRange{Int64}}=nothing,
              hcr::Int=1)
    if ws === nothing
        C = 1
    else
        C = length(ws)
    end
    hall = Tuple(fill(-Nhmax:Nhmax, C))
    hall = Iterators.product(Tuple(fill(-Nhmax:Nhmax, C))...)

    if hcr==1
        h = [zeros(Int, 1, C); I(C)]
        for hi in hall
            if sum(abs.(hi)) <= Nhmax && sum(hi) >= 0
                h = vcat(h, [i for i in hi]')
            end
        end
        h = unique(h, dims=1)
    else
        error("Unknown criterion");
    end
    h = [zeros(Int, 1, C); h[.!(vec(sum(h, dims=2)).==0 .& h[:,1].<=0), :]];
    h = Int.(sign.(h[:,1].+.5).*h);
end

# **  Evaluate Fourier Series

"""
  FSEVAL(U, h, t)

  Evaluate Fourier series at given points t


"""
function FSEVAL(h::hTypes, t::Union{StepRangeLen{Float64}, VecOrMat{Float64}},
                U::Union{Nothing,Matrix{Float64}}=nothing)
    Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));
    Nt = size(t, 1);

    _, hn0i, zinds, rinds, iinds = HINDS(1, h)

    cbass = cos.(h[hn0i,:]*t');
    sbass = sin.(h[hn0i,:]*t');

    J = zeros(Nt, Nhc);
    J[:, zinds] .= 1.0;
    J[:, [rinds; iinds]] = [cbass;sbass]';

    if U===nothing
        return J
    else
        cbass = reshape(cbass, (size(cbass)..., 1))
        sbass = reshape(sbass, (size(sbass)..., 1))
        cbass = permutedims(cbass, (1,3,2));
        sbass = permutedims(sbass, (1,3,2));

        Np = size(U, 2);
        ut = zeros(Nt, Np);
        ut = dropdims(U[zinds,:] .+ sum(U[rinds,:].*cbass + U[iinds,:].*sbass, dims=1), dims=1)';
    
        return ut, J
    end
end


# **  Harmonic Stiffness

"""
    HARMONICSTIFFNESS(M, D, K, ws::Array{Float64}, h)

    Computes the harmonic stiffness for HB simulations.

# Arguments
+ 'M, D, K': (Nd,Nd) Mass, Damping, Stiffness Matrices.
+ 'ws::Array{Float64}': (C,1) list of frequencies
+ 'h': (H, C) harmonic indices
# Outputs
+ `E`: (Nd*Nhc, Nd*Nhc) Harmonic Stiffness
+ `dEdw`: C-Vector of (Nd*Nhc, Nd*Nhc) Representing Frequency-gradients of E
# Examples
+ Single Frequency Case
```julia-repl
M = I(2);
K = [2 -1;-1 2];
D = 0.001.*K + 0.01.*M;

h = [0;1;2;3];
ws = 1;
E, dEdw = HARMONICSTIFFNESS(M, D, K, ws, h);
```
+ 2-Frequency Case
```julia-repl
Nhmax = 3;
ws = [1; π];
h = HSEL(Nhmax, ws);
E, dEdw = HARMONICSTIFFNESS(M, D, K, ws, h);
```
"""
function HARMONICSTIFFNESS(M::MxTypes,
                           D::MxTypes,
                           K::MxTypes,
                           ws::Union{Float64, Vector{Float64}},
                           h::hTypes)
    H = size(h,1);
    C = size(h,2);
    Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));

    hn = kron(h, ones(2,1));
    E1 = kron(sparse(I, H, H), [0 1;-1 0])
    if all(h[1,:].==0)
        hn = hn[2:end,:];
        E1 = E1[2:end, 2:end];
    end
    E = kron(sparse(I, Nhc, Nhc), K) + kron((hn*ws).*E1, D) - kron(spdiagm((hn*ws)[:].^2), M);
    dEdw = [kron(hn[:,i].*E1, D) -
        kron(spdiagm(2(hn*ws)[:] .* hn[:,i]), M) for i in 1:C];

    return (E, dEdw);
end

# ---------------------------------------------------------

"""
    HARMONICSTIFFNESS!(E, dEdw, M, D, K, ws, h)

    Bang version of HARMONICSTIFFNESS (preallocate E, dEdw).
    Advantage is that nothing can be passed to avoid needless computations. 

# Arguments
+ `E`: (Nd*Nhc, Nd*Nhc) Harmonic Stiffness
+ `dEdw`: C-Vector of (Nd*Nhc, Nd*Nhc) Representing Frequency-gradients of E
+ 'M, D, K': (Nd,Nd) Mass, Damping, Stiffness Matrices.
+ 'ws::Array{Float64}': (C,1) list of frequencies
+ 'h': (H, C) harmonic indices
# Examples
+ Single Frequency Case
```julia-repl
M = I(2);
K = [2 -1;-1 2];
D = 0.001.*K + 0.01.*M;

h = [0;1;2;3];
Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));
E = zeros(2*Nhc, 2*Nhc);
dEdw = zeros(2*Nhc, 2*Nhc);
ws = 1;
HARMONICSTIFFNESS!(E, dEdw, M, D, K, ws, h);
```
+ 2-Frequency Case
```julia-repl
Nhmax = 3;
ws = [1; π];
h = HSEL(Nhmax, ws);
Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));
E = zeros(2*Nhc, 2*Nhc);
dEdw = zeros(2*Nhc, 2*Nhc);
HARMONICSTIFFNESS!(E, dEdw, M, D, K, ws, h);
```
"""
function HARMONICSTIFFNESS!(E, dEdw, M, D, K, ws, h)
    H = size(h,1);
    C = size(h,2);
    Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));

    hn = kron(h, ones(2,1));
    E1 = kron(sparse(I, H, H), [0 1;-1 0])
    if all(h[1,:].==0)
        hn = hn[2:end,:];
        E1 = E1[2:end, 2:end];
    end
    if !(E === nothing)
        E[:,:] = kron(sparse(I, Nhc, Nhc), K) + kron((hn*ws).*E1, D) -
            kron(spdiagm((hn*ws)[:].^2), M);
    end
    if !(dEdw === nothing)
        tmp = hn*ws;
        dEdw[:] = [kron(hn[:,i].*E1, D) -
            kron(spdiagm(2tmp[:] .* hn[:,i]), M) for i in 1:C];
    end
end

# **  Harmonic Indices

"""
# Description

# Arguments
- Ndofs::Int64     : 
- h::hTypes :
# Outputs
- zinds
- hinds
- rinds0
- rinds
- iinds
"""
function HINDS(Ndofs::Int64, h::hTypes)
    Nh = size(h,1);
    if h[1] == 0
        rinds0 = collect(1:Ndofs);
        zinds = collect(1:Ndofs);
        hinds = collect((Ndofs+1):(Ndofs*Nh));

        rinds = vec(collect(1:Ndofs) .+ 2Ndofs*(0:Nh-2)' .+ Ndofs);
        iinds = rinds .+ Ndofs;
    else
        rinds0 = [];
        zinds = [];
        hinds = collect(1:Ndofs*Nh);

        rinds = vec(collect(1:Ndofs) .+ 2Ndofs*(0:Nh-1)');
        iinds = rinds .+ Ndofs;
    end

    return (zinds, hinds, rinds0, rinds, iinds);
end

# **  Differentiation Matrix

"""
  DFOUR(h::hTypes, ws=nothing)

  Returns the Fourier differentiation matrix (useful for HB representations).

# Arguments
+ h::hTypes : (H, C) Harmonic indices
+ ws::Array{Float64} : (C, 1) (optional) list of frequencies
# Outputs
+ `D` : (Nhc, Nhc) Harmonic differentiation matrix
+ `dDdw` : (Nhc, Nhc) Harmonic differentiation matrix (derivative wrt ws).
"""
function DFOUR(h::hTypes, ws=nothing)
    H = size(h, 1);
    C = size(h, 2);
    if ws === nothing
        ws = ones(C,1);
    end

    hn = kron(h, ones(2, 1));
    E1 = kron(sparse(I, H, H), [0 1; -1 0]);
    if all(h[1, :] .== 0)
        hn = hn[2:end, :];
        E1 = E1[2:end, 2:end];
    end
    D = (hn*ws).*E1;
    dDdw = [hn[:, i] .* E1 for i in 1:C];
    return (D, dDdw);
end

"""
  DFOUR!(h::hTypes, D, dDdw=nothing, ws=nothing)

  Bang version of DFOUR.

# Arguments
+ h::hTypes : (H, C) Harmonic indices
+ `D` : (Nhc, Nhc) Harmonic differentiation matrix
+ `dDdw` : (Nhc, Nhc) Harmonic differentiation matrix (derivative wrt ws).
+ ws::Array{Float64} : (C, 1) (optional) list of frequencies
"""
function DFOUR!(h::hTypes, D, dDdw=nothing, ws=nothing)
    H = size(h, 1)
    C = size(h, 2)
    if ws === nothing
        ws = ones(C,1);
    end

    hn = kron(h, ones(2, 1))
    E1 = kron(sparse(I, H, H), [0 1; -1 0])
    if all(h[1, :] .== 0)
        hn = hn[2:end, :]
        E1 = E1[2:end, 2:end]
    end
    if !(D === nothing)
        D[:, :] = (hn*ws).*E1;
    end
    if ~(dDdw === nothing)
        dDdw[:] = [hn[:, i] .* E1 for i in 1:C];
    end
end

# **  Fourier Product Matrix

"""
# Description

# Arguments
- U    : 
- h    : 
- Hmax : (default nothing)
- D    : (default nothing)
- Lb   : (default nothing)
"""
function PRODMAT_FOUR(U, h::hTypes, Hmax=nothing, D=nothing, L=nothing)
    if Hmax === nothing
        Hmax = maximum(h);
    end
    hfull = HSEL(Hmax);

    hfri = vec((h.-1)'.*2 .+ [1; 2]);
    hfri = hfri[hfri.>=0];

    Lb = I(2*Hmax+1);
    Lb = Lb[hfri.+1, :];

    U = Lb'*U;
    h = hfull;

    # Compute D
    Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));
    a0 = U[1];
    as = [0; U[2:2:end]];
    bs = [0; U[3:2:end]];

    zind = 1;
    cinds = 2:2:Nhc;
    sinds = 3:2:Nhc;

    zh = zeros(length(h)-1,1);

    Df = zeros(Nhc, Nhc);
    Df[zind, [zind; cinds; sinds]] = [a0 as[2:end]'/2 bs[2:end]'/2];

    Df[cinds, zind] = as[2:end];
    tmp = (Toeplitz(as,as)+Hankel(as, vec([as[end]; zh])))/2;
    Df[cinds, cinds] = tmp[2:end, 2:end];
    tmp = -(Toeplitz(bs,-bs)-Hankel(bs, vec([bs[end]; zh])))/2;
    Df[cinds, sinds] = tmp[2:end, 2:end];

    Df[sinds, zind] = bs[2:end];
    tmp = (Toeplitz(bs,-bs)+Hankel(bs, vec([bs[end]; zh])))/2;
    Df[sinds, cinds] = tmp[2:end, 2:end];
    tmp = (Toeplitz(as,as)-Hankel(as, vec([as[end]; zh])))/2;
    Df[sinds, sinds] = tmp[2:end, 2:end];

    Df[2:end, 2:end] = Df[2:end, 2:end] .+ I(Nhc-1).*a0;

    #
    if (D === nothing)
        return (Lb*Df*Lb');
    elseif (L === nothing)
        D[:, :] = Lb*Df*Lb';
    else
        L[:, :] = Lb;
        D[:, :] = Df;
    end
end

# * Chebyshev Routines ##########################################################

# ** Alternating Chebyshev-"Time" Transform
"""
  ACT(yin, h::hTypes, N::Int64, dir::Symbol)

  Routine for the alternating Chebyshev Transform. Uses AFT internally. 

# Arguments
- yin              : 
- h::hTypes : 
- N::Int64         : 
- dir::Symbol      : 
"""
function ACT(yin::VecOrMat{Float64}, h::hTypes, N::Int64, dir::Symbol)
    Nhc = sum(all(h.==0, dims=2) + 2*any(h.!=0, dims=2));
    # Chebyshev only uses cosines
    L = I(Nhc);
    if all(h[1, :] .== 0)
        L = L[[1; 2:2:end],:];
    else
        L = L[2:2:end, :];
    end
    Ny = size(yin, 2);
    C = size(h,2);
    
    if cmp(dir, :t2f) == 0
        Nt = 2N-2;

        rinds = Tuple([fill([1:N; N-1:-1:2], C); :]);
        rNs = Tuple([fill(N, C); Ny]);
        yin = reshape(reshape(yin, rNs)[rinds...], :, Ny);
        yout = AFT(yin, h, Nt, dir);
        yout = L*yout;
    elseif cmp(dir, :f2t) == 0
        Nt = 2(N-1);
        yout = AFT(L'*yin, h, Nt, dir);
        
        tN = Tuple([fill(Nt, C); Ny]);
        tNs = Tuple([fill(1:N, C); :]);
        yout = reshape(reshape(yout, tN)[tNs...], :, Ny);        
    else
        error("Unknown dir for ACT");
    end
    if size(yout,2)==1
        yout = vec(yout);
    end
    return yout;
end

# **  Chebyshev Differentiation Matrix

"""
  DCHEB(h::hTypes)

  Returns the Chebyshev Differentiation matrix.

# Arguments
- h::hTypes : 
"""
function DCHEB(h::hTypes)
    H = size(h, 1);
    C = size(h, 2);
    if C!=1
        error("DCHEB Currently only implemented for C=1.");
    end

    D = zeros(H, H);
    for hi in h[h.!=0]
        if mod(hi, 2) == 0  # even
            D[2:2:hi, 1+hi] .= 2*hi;
        else  # odd
            D[1:2:hi, 1+hi] .= 2*hi;
            D[1, 1+hi] = D[1, 1+hi] .- hi;
        end
    end
    return D;
end

# **  Chebyshev Product Matrix

"""
# Description

# Arguments
- U    : 
- h    : 
- Hmax : (default nothing)
- D    : (default nonothing)
- L    : (default nonothing)
"""
function PRODMAT_CHEB(U, h, Hmax=nothing, D=nothing, L=nothing)
    if Hmax === nothing
        Hmax = maximum(h)
    end

    hfull = collect(0:Hmax);
    Lb = I(Hmax+1);
    Lb = Lb[1 .+ h, :];

    U = Lb'*U;
    h = hfull;

    ## Compute D
    Nh = length(h);
    a0 = U[1];
    as = [0; U[2:end]];

    zind = 1;
    cinds = 2:Nh;

    Df = zeros(Nh, Nh);
    Df[zind, [zind; cinds]] = [a0 as[2:end]'/2];
    Df[cinds, zind] = as[2:end];
    tmp = (Toeplitz(as,as)+Hankel(as,vec([as[end]; zeros(Nh-1,1)])))/2;
    Df[cinds, cinds] = tmp[2:end, 2:end] + I(Nh-1).*a0;

    ##
    if D === nothing
        return Lb*Df*Lb';
    elseif L === nothing
        D[:, :] = Lb*Df*Lb;
    else
        D[:, :] = Df;
        L[:, :] = Lb;
    end
end
