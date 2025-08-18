# # [Example d0](@id ex_d0)
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
β = 0.1;
fnl = (t,u,ud)->return β.*u.^3, 3β.*u.^2, zeros(size(u));;
L = [0.0 1.0];

mdl = ADDNL(mdl, :Inst, fnl, L);

# ## Trial
h = HSEL(3, 1.0);
## h = 1:2:5;
h = 0:5;
N = 128;
t = (0:N-1)*2π/N;

Nhc = sum(all(h.==0, dims=2) + 2any(h.!=0, dims=2));

_, _, zinds, rinds, iinds = HINDS(mdl.Ndofs, h)
Fl = zeros(Nhc*mdl.Ndofs, 1);
Fl[rinds[1]] = 1.0;

Wst = 0.6;
E, dEdw = HARMONICSTIFFNESS(mdl.M, mdl.C, mdl.K, Wst, h);

Uw0 = [E\Fl; Wst];

# ## Nonlinear Evaluation
FNL = zeros(mdl.Ndofs*Nhc);
dFNLdU = zeros(mdl.Ndofs*Nhc, mdl.Ndofs*Nhc);
dFNLdw = zeros(mdl.Ndofs*Nhc);
NLEVAL!(Uw0, mdl, h, N; FNL=FNL, dFNLdU=dFNLdU, dFNLdw=dFNLdw)

# ## HB Residue
R = zeros(mdl.Ndofs*Nhc);
dRdU = zeros(mdl.Ndofs*Nhc, mdl.Ndofs*Nhc);
dRdw = zeros(mdl.Ndofs*Nhc);
HBRESFUN!(Uw0, mdl, Fl, h, N; R=R, dRdU=dRdU, dRdw=dRdw)

fun = NonlinearFunction((r,u,p)->HBRESFUN!([u;p], mdl, Fl, h, N; R=r),
    jac=(J,u,p)->HBRESFUN!([u;p], mdl, Fl, h, N; dRdU=J),
    paramjac=(Jp,u,p)->HBRESFUN!([u;p], mdl, Fl, h, N; dRdw=Jp));

prob = NonlinearProblem(fun, Uw0[1:end-1], Uw0[end]);
sol = solve(prob, show_trace=Val(true));

# ## Continuation
Om0 = 0.1;
Om1 = 3;
dOm = 0.01;

cpars = (parm=:arclength, nmax=2000, Dsc=:auto, minDsc=1e-5);
sols, its, dss, xis, Dsc = CONTINUATE(Uw0[1:end-1], fun, [Om0, Om1], dOm; cpars...);

uh = zeros(Complex, maximum(h)+1, length(sols), 2);
for i in 1:2
    uh[h.+1, :, i] = hcat([[s.up[zinds[i:2:end]];
                            s.up[rinds[i:2:end]]+1im*s.up[iinds[i:2:end]]]
                          for s in sols]...);
end
Oms = [s.up[end] for s in sols];

# ## Plot

his = [1, 3, 5];

fsz = 24;
fig = Figure(fontsize=fsz, size=(1000, 600));
if !isdefined(Main, :scr) && isdefined(Main, :GLMakie) #src
    scr = GLMakie.Screen(); #src
end #src

for i in eachindex(his[his.<=maximum(h)])
    ax = Axis(fig[1, i],
        ylabel=L"$H_%$(his[i])$ Response (m)", yscale=log10);
    lines!(ax, Oms, abs.(uh[his[i].+1, :, 1]), label="x1");
    lines!(ax, Oms, abs.(uh[his[i].+1, :, 2]), label="x2");

    ax = Axis(fig[2, i], xlabel=L"Excitation Frequency $\Omega$",
        ylabel=L"$H_%$(his[i])$ Phase (rad)");
    lines!(ax, Oms, unwrap(angle.(uh[his[i].+1, :, 1])), label="x1");
    lines!(ax, Oms, unwrap(angle.(uh[his[i].+1, :, 2])), label="x2");
end

if isdefined(Main, :GLMakie) #src
    display(scr, fig); #src
else #src
    display(fig) #src
    fig
end #src
