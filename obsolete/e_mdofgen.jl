using Revise
using LinearAlgebra
using SparseArrays
using NLsolve
include("../ROUTINES/MDOFGEN.jl")

M = I(2);
K = [2 -1;-1 2];
C = 0.01*M+0.001*K;

MDL = MDOFGEN(M, C, K);

# Nonlinearity
β = 0.1;
fnl = (t,u,ud)->return β.*u.^3, 3β.*u.^2, zeros(size(u));;
L = [0.0 1.0];

MDL = ADDNL(MDL, :Inst, fnl, L);

# * Trial
h = HSEL(3, 1);
N = 128;
t = (0:N-1)*2π/N;

Nhc = sum(all(h.==0, dims=2) + 2*any(h.!=0, dims=2));

_, _, rinds0, rinds, iinds = HINDS(MDL.Ndofs, h);
Fl = zeros(Nhc*MDL.Ndofs,1);
Fl[rinds[1]] = 1.0;

E = spzeros(Nhc*MDL.Ndofs, Nhc*MDL.Ndofs);
dEdw = [spzeros(Nhc*MDL.Ndofs, Nhc*MDL.Ndofs)];

Wst = 0.6;
HARMONICSTIFFNESS!(E, dEdw, MDL.M, MDL.C, MDL.K, [Wst], h);

Uw0 = [E\Fl; Wst];

# *
FNL = spzeros(MDL.Ndofs*Nhc,1);
dFNLdU = spzeros(MDL.Ndofs*Nhc,MDL.Ndofs*Nhc);
dFNLdw = spzeros(MDL.Ndofs*Nhc,1);
NLEVAL!(FNL, dFNLdU, dFNLdw, MDL, Uw0, h, N, 1e-6);

# *
R = spzeros(MDL.Ndofs*Nhc,1);
dRdU = spzeros(MDL.Ndofs*Nhc,MDL.Ndofs*Nhc);
dRdw = spzeros(MDL.Ndofs*Nhc,1);
HBRESFUN!(R, dRdU, dRdw, MDL, Uw0, Fl, h, N, 1e-6);

HBRESFUN!(R, dRdU, nothing, MDL, Uw0, Fl, h, N, 1e-6);

res = nlsolve(only_fj!((F,J,x)->HBRESFUN!(F, J, nothing, MDL, [x;Wst], Fl, h, N, 1e-6)),
              Uw0[1:end-1], show_trace=true);
