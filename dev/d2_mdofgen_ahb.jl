using GLMakie
using LinearAlgebra
using SparseArrays
using ForwardDiff
using NonlinearSolve
using DSP

using Revise
using juliajim.HARMONIC
using juliajim.CONTINUATION
using juliajim.MDOFUTILS

# * System Setup
M = I(2);
K = [2 -1;-1 2];
C = 0.01*M+0.001*K;

mdl = MDOFGEN(M, C, K);

# Nonlinearity
kt = 1.0;
fs = 1.0;
fnl = (t,u,up,fp)-> if all(abs.(fp+kt*(u-up)).<fs)
    return fp+kt*(u-up), kt*ones(size(u)), -kt*ones(size(u)), ones(size(u)); else
        return fs*sign.(fp+kt*(u-up)), zeros(size(u)), zeros(size(u)), zeros(size(u));
end
L = [0.0 1.0];

mdl = ADDNL(mdl, :Hyst, fnl, L);

# * Trial
h = HSEL(3, 1.0);
h = 1:2:5;
# h = 0:5;
N = 256;
t = (0:N-1)*2π/N;

Nhc = sum(all(h.==0, dims=2) + 2any(h.!=0, dims=2));

_, _, zinds, rinds, iinds = HINDS(mdl.Ndofs, h)
Fl = zeros(Nhc*mdl.Ndofs);
Fl[rinds[1]] = 1.0;

Wst = 0.6;
E, dEdw = HARMONICSTIFFNESS(mdl.M, mdl.C, mdl.K, Wst, h);
U0 = E\Fl;

# * HB Residue
R = zeros(mdl.Ndofs*Nhc+1);
dRdUf = zeros(mdl.Ndofs*Nhc+1, mdl.Ndofs*Nhc+1);
dRdw = zeros(mdl.Ndofs*Nhc+1);
Amp = 1e0;
dof = 1;

HBRESFUN_A!([U0; 1.0; Wst], mdl, Amp, Fl, h, N; R=R, dRdUf=dRdUf, dRdw=dRdw)

fun = NonlinearFunction((r,uf,p)->HBRESFUN_A!([uf;p], mdl, Amp, Fl, h, N; R=r),
    jac=(J,uf,p)->HBRESFUN_A!([uf;p], mdl, Amp, Fl, h, N; dRdUf=J),
    paramjac=(Jp,uf,p)->HBRESFUN_A!([uf;p], mdl, Amp, Fl, h, N; dRdw=Jp));

prob = NonlinearProblem(fun, [U0;1.0], Wst);
sol = solve(prob, show_trace=Val(true));

# * Continuation
Om0 = 0.1;
Om1 = 3;
dOm = 0.1;

cpars = (parm=:arclength, nmax=400, Dsc=:none, itopt=4);
sols, its, dss, xis, Dsc = CONTINUATE([U0;1.0], fun, [Om0, Om1], dOm; cpars...);

uh = zeros(Complex, maximum(h)+1, length(sols), 2);
for i in 1:2
    uh[h.+1, :, i] = hcat([[s.up[zinds[i:2:end]];
                            s.up[rinds[i:2:end]]+1im*s.up[iinds[i:2:end]]]
                          for s in sols]...);
end
Oms = [s.up[end] for s in sols];
Fs = [s.up[end-1] for s in sols];

# Plot
his = [1, 3, 5];

set_theme!(theme_black())
fsz = 24;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr) && isdefined(Main, :GLMakie)
   scr = GLMakie.Screen();
end

for i in eachindex(his[his.<=maximum(h)])
    ax = Axis(fig[1, i], xlabel=L"Excitation Frequency $\Omega$",
              ylabel=L"$H_%$(his[i])$ Response (m/N)", yscale=log10);
    scatterlines!(ax, Oms, abs.(uh[his[i].+1, :, 1])./Fs, label="x1");
    scatterlines!(ax, Oms, abs.(uh[his[i].+1, :, 2])./Fs, label="x2");

    ax = Axis(fig[2, i], xlabel=L"Excitation Frequency $\Omega$",
              ylabel=L"$H_%$(his[i])$ Phase (rad)");
    scatterlines!(ax, Oms, unwrap(angle.(uh[his[i].+1, :, 1])), label="x1");
    scatterlines!(ax, Oms, unwrap(angle.(uh[his[i].+1, :, 2])), label="x2");
end

if isdefined(Main, :GLMakie)
   display(scr, fig);
else
    fig
end

