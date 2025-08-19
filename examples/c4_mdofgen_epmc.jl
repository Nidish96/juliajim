# # [Example C4: Nonlinear Normal Modes of a 2DoF Hysteretic Oscillator using EPMC](@id ex_c4)
# ## Preamble: Load Packages
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

# ## System Setup
M = I(2);
K = [2 -1;-1 2];
C = 0.01*M+0.001*K;

mdl = MDOFGEN(M, C, K);

## Nonlinearity
kt = 1.0;
fs = 1.0;
fnl = (t,u,up,fp)-> if all(abs.(fp+kt*(u-up)).<fs)
    return fp+kt*(u-up), kt*ones(size(u)), -kt*ones(size(u)), ones(size(u)); else
        return fs*sign.(fp+kt*(u-up)), zeros(size(u)), zeros(size(u)), zeros(size(u));
end
L = [0.0 1.0];

mdl = ADDNL(mdl, :Hyst, fnl, L);

# ## Linear Modal Analysis
Knl0 = L'kt*L;
K0 = mdl.K+Knl0;
D, V = eigen(K0, mdl.M)
W0s = sqrt.(D);
Xis = diag(V'mdl.C*V)

mi = 1;  # Mode of interest

# ## Trial
h = HSEL(3, 1.0);
h = 1:2:3;
h = 0:5;
N = 256;
t = (0:N-1)*2π/N;

Nhc = sum(all(h.==0, dims=2) + 2any(h.!=0, dims=2));

_, _, zinds, rinds, iinds = HINDS(mdl.Ndofs, h)
Fl = zeros(Nhc*mdl.Ndofs);
Fl[rinds[1]] = 1.0;

E, dEdw = HARMONICSTIFFNESS(mdl.M, mdl.C, mdl.K, W0s[mi], h);

U0 = zeros(mdl.Ndofs*Nhc);
U0[iinds[1:mdl.Ndofs]] = V[:, mi];

# ## EPMC Residue
R = zeros(mdl.Ndofs*Nhc+2);
dRdU = zeros(mdl.Ndofs*Nhc+2, mdl.Ndofs*Nhc+2);
dRda = zeros(mdl.Ndofs*Nhc+2);

EPMCRESFUN!([U0; W0s[mi]; Xis[mi]; -2], mdl, Fl, h, N; R=R, dRdUwx=dRdU, dRda=dRda)

# ## 
A0 = -2.0;

U0 = zeros(mdl.Ndofs*Nhc);
## U0[rinds[1:mdl.Ndofs]] = (10^A0)*V[:, mi];
U0[iinds[1:mdl.Ndofs]] = V[:, mi];

fun = NonlinearFunction((r,uwx,p)->EPMCRESFUN!([uwx;p], mdl, Fl, h, N; R=r),
    jac=(J,uwx,p)->EPMCRESFUN!([uwx;p], mdl, Fl, h, N; dRdUwx=J),
    paramjac=(Jp,uwx,p)->EPMCRESFUN!([uwx;p], mdl, Fl, h, N; dRda=Jp));

prob = NonlinearProblem(fun, [U0;W0s[mi];Xis[mi]], A0);
sol = solve(prob, show_trace=Val(true), maxiters=100);

# ## Continuation
A0 = -1.;
A1 = 3.;
da = 0.05;

cpars = (parm=:arclength, nmax=1000, Dsc=:none, itopt=4, dpbnds=[da/5, 2da]);
sols, its, dss, xis, Dsc = CONTINUATE([U0;W0s[mi];Xis[mi]], fun, [A0, A1], da; cpars...);

uh = zeros(Complex, maximum(h)+1, length(sols), 2);
for i in 1:2
    uh[h.+1, :, i] = hcat([[s.up[zinds[i:2:end]];
                            s.up[rinds[i:2:end]]+1im*s.up[iinds[i:2:end]]]
                          for s in sols]...);
end
Oms = [s.up[end-2] for s in sols];
Zts = [s.up[end-1] for s in sols]./2Oms;
As = [10^s.up[end] for s in sols];

# ## Plot

his = [1, 3, 5];

set_theme!(theme_latexfonts());
fsz = 20;
fig = Figure(fontsize=fsz, size=(1000, 600));
if !isdefined(Main, :scr) && isdefined(Main, :GLMakie) #src
   scr = GLMakie.Screen(); #src
end #src

axs = [];
for i in eachindex(his[his.<=maximum(h)])
    ax = Axis(fig[1, i],
        ylabel=L"$H_%$(his[i])$ Response (m)", xscale=log10, yscale=log10);
    scatterlines!(ax, As, abs.(uh[his[i].+1, :, 1]), label="x1");
    scatterlines!(ax, As, abs.(uh[his[i].+1, :, 2]), label="x2");
    push!(axs, ax)

    ax = Axis(fig[2, i], xlabel=L"Modal Amplitude $a$",
        ylabel=L"$H_%$(his[i])$ Phase (rad)", xscale=log10);
    scatterlines!(ax, As, unwrap(angle.(uh[his[i].+1, :, 1])), label="x1");
    scatterlines!(ax, As, unwrap(angle.(uh[his[i].+1, :, 2])), label="x2");
    push!(axs, ax)
end
ax = Axis(fig[1, length(his)+1],
    ylabel=L"Natural Frequency $\omega_n$ (rad/s)", xscale=log10);
scatterlines!(ax, As, Oms);
push!(axs, ax)

ax = Axis(fig[2, length(his)+1], xlabel=L"Modal Amplitude $a$",
    ylabel="Damping Factor (%)", xscale=log10);
scatterlines!(ax, As, 100Zts);
push!(axs, ax)

linkxaxes!(axs...)

if isdefined(Main, :GLMakie) #src
   display(scr, fig); #src
else #src
    fig
end #src
