# ```@meta
# CurrentModule = juliajim
# ```

# ## * [Example D3: Branch Switching Using Normal Forms](@id ex_d3)

# This example is a direct continuation of [Example D2](@ref ex_d2). In that example we did a branch switching analysis by arbitrarily choosing the amplitude of perturbation. Here, we use the [`NORMALFORMFIT`](@ref juliajim.NOLDYN.NORMALFORMFIT) function to fit the system to a cubic nonlinear "normal form" (through a similarity transformation eliminating certain terms).

# The normal form can easily be interpreted to provide an estimate for the bifurcated branch amplitudes. It can be seen that the normal form approach provides a near exact estimate of the bifurcated branch amplitude in the end of this example.

# ## Preamble: Load Packages
using GLMakie
using LinearAlgebra
using SparseArrays
using NonlinearSolve
using DSP
using Debugger

using Revise
using juliajim.HARMONIC
using juliajim.CONTINUATION
using juliajim.MDOFUTILS
using juliajim.NLDYN

# ## System Setup
# Here we setup the system with parameters.

pars = (c=0.02, k=4., F=1., eta=0.1);

mdl = MDOFGEN(1.0, -pars.c, pars.k);
fnl = (t,u,ud) -> return pars.eta.*u.^2 .*ud, 2pars.eta.*u.*ud, pars.eta.*u.^2;
mdl = ADDNL(mdl, :Inst, fnl, 1.0);

# ## Setup HB

h = 0:5;
N = 128;
t = range(0, 2π, N+1)[1:N];

Nhc = NHC(h);

inds0, hinds, zinds, rinds, iinds = HINDS(mdl.Ndofs, h);
Fl = zeros(Nhc*mdl.Ndofs, 1);
Fl[rinds[1]] = pars.F;

Om0 = 0.1;
Om1 = 2.0;
E, dEdw = HARMONICSTIFFNESS(mdl.M, mdl.C, mdl.K, Om0, h);

Uw0 = [E\Fl; Om0];

# ## Setup HB Residue, Get First Point
fun = NonlinearFunction((r,u,p)->HBRESFUN!([u;p], mdl, Fl, h, N; R=r),
    jac=(J,u,p)->HBRESFUN!([u;p], mdl, Fl, h, N; dRdU=J),
    paramjac=(Jp,u,p)->HBRESFUN!([u;p], mdl, Fl, h, N; dRdw=Jp));

prob = NonlinearProblem(fun, Uw0[1:end-1], Uw0[end]);
sol = solve(prob, show_trace=Val(true));

# ## Continuation
Om0 = 0.1;
Om1 = 4.0;
dOm = 0.2;
cpars = (parm=:arclength, nmax=100, save_jacs=true);

sols, _, _, _, _ = CONTINUATE(Uw0[1:end-1], fun, [Om0, Om1], dOm; cpars...);

# Obtain Harmonics
uh = zeros(Complex, maximum(h)+1, length(sols));
uh[[inds0; hinds], :] = hcat([[up[zinds];up[rinds,:]+1im*up[iinds,:]] for up in sols.up]...);
Oms = [up[end] for up in sols.up];

# ### Stability Analysis

E0, _ = HARMONICSTIFFNESS(0., -2mdl.M, 0, 1, h);
E0 = collect(E0);
stab = zeros(length(Oms));
for (i,(J,Om)) in enumerate(zip(sols.J, Oms))
    evs = eigvals(J[2:3,2:3], Om*E0[2:3,2:3]);
    stab[i] = sum(real(evs).>=0);
end

# ## Bifurcation Analysis

# ### Setup QPHB for Branch Switching
Nhmax = 3;
hq = HSEL(Nhmax, [1.,1.]);

inds0q,hindsq, zindsq,rindsq,iindsq = HINDS(1, hq);
Nhcq = NHC(hq);
Nq = 32;

h01 = findall(vec(all(hq.==0, dims=2)));
hq1s = setdiff(findall(hq[:,2].==0), h01).-length(h01);  # w1 alone
hq2s = setdiff(findall(hq[:,1].==0), h01).-length(h01);  # w2 alone

Flq = zeros(Nhcq);
Flq[rindsq[hq2s[1]]] = pars.F;

# ### Detection
bifis = findall(stab[1:end-1].!=stab[2:end]);
bifis[2] += 1;
dxis = [-1, 1];

bi = 2;
bifi = bifis[bi];
eVals, eVecs = eigen(sols.J[bifi][2:3, 2:3], Oms[bifi]*collect(E0[2:3,2:3]));
eVecsC = eVecs[1,:]-1im.*eVecs[2,:];  # Complexify

ei = argmax(abs.(eVecsC));  # Only one should "survive" for the Hopf bifurcation 

sig = imag(eVals[ei]);
Wself = Oms[bifi] + sig;  # Self Excited Frequency
Pvec = normalize([real(eVecsC[ei]), -imag(eVecsC[ei])]);  # Perturbation Vector

# Setup Initial Guess
Uhq0 = zeros(Nhcq);
Uhq0[[rindsq[hq2s]; iindsq[hq2s]]] = sols.up[bifi][[rinds[1:length(hq2s)];
                                                    iinds[1:length(hq2s)]]];

# ### Fit the HB Residue to a Local Normal form

Fl1 = Fl[[rinds[1], iinds[1]]];  # Single harmonic excitation
uh1 = sols.up[bifi][[rinds[1], iinds[1]]];  # Bifurcated Solution
vh1 = [normalize(real(eVecs[:,ei])) normalize(imag(eVecs[:,ei]))];  # Projector

uh0 = sols.up[bifi-dxis[bi]][[rinds[1], iinds[1]]];  # Previous Solution
d_amp = norm(vh1'*(uh1-uh0));

E0f, _ = HARMONICSTIFFNESS(0, -2mdl.M, 0, Oms[bifi], 1)
xyfun = (u) -> vh1\ (E0f\HBRESFUN!([uh1+vh1*u; Oms[bifi]], mdl, Fl1, 1, N));

A, B, C, H = NORMALFORMFIT(xyfun, d_amp/10, 100);
λ = sum(diag(A))/2;
α = sum(C[[1,4,5,8]])/4;

pVal, pVec = eigen(A);
vi = argmax(abs.(pVec[1,:]+1im*pVec[2,:]));
Pvec = vh1*pVec[:,vi];
Pvec = normalize([real(pVec[1]+1im*pVec[2]), imag(pVec[1]+1im*pVec[2])]);

# Perturbation Vector
Phq0 = zeros(Nhcq);
Phq0[[rindsq[hq1s[1]]; iindsq[hq1s[1]]]] = Pvec;
# Perturbation Amplitude
qamp = √(-λ/α);  

# ### Apply phase constraint and converge with deflation
cL = I(Nhcq+2)[:, setdiff(1:Nhcq+2, iindsq[hq1s[1]])];
uC = cL'*[Uhq0; Wself; Oms[bifi]];

funq = NonlinearFunction((r,u,p)-> QPHBRESFUN!([u;p], mdl, Flq, hq, Nq; R=r,cL=cL,U0=uC),
    jac=(J,u,p)->QPHBRESFUN!([u;p], mdl, Flq, hq, Nq; dRdU=J,cL=cL,U0=uC),
    paramjac=(Jp,u,p)->QPHBRESFUN!([u;p], mdl, Flq, hq, Nq; dRdw=Jp,cL=cL,U0=uC));

u0 = cL'*[Uhq0+qamp*Phq0; Wself; Oms[bifi]];
probq = NonlinearProblem(funq, u0[1:end-1], u0[end]);
solq = solve(probq, show_trace=Val(true));

# We can now verify that the converged solution does indeed have an amplitude that's almost exactly equal to the normal form prediction.
[cL[2,:]'*[solq.u; Oms[bifi]], qamp]
