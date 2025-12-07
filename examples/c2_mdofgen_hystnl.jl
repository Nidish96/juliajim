# ```@meta
# CurrentModule = juliajim
# ```

# # [Example C2: Forced Response of a 2DoF Oscillator with Hysteretic Nonlinearity](@id ex_c2)

# Much like Examples [B1](@ref ex_b1) and [B2](@ref ex_b2), this builds off of [Example C1](@ref ex_c1) to conduct the forced response analysis of a 2DoF oscillator with a hysteretic support. 

# The dynamical system that will be studied here is:
# ```math
# \underbrace{\begin{bmatrix} 1&0\\0&1 \end{bmatrix}}_{\mx{M}} \begin{bmatrix} \ddot{x_1}\\ \ddot{x_2} \end{bmatrix} +
# \left( 0.01 \mx{M} + 0.001 \mx{K}\right) \begin{bmatrix} \dot{x_1}\\ \dot{x_2} \end{bmatrix} + \underbrace{\begin{bmatrix} 2&-1\\-1&2 \end{bmatrix}}_{\mx{K}} \begin{bmatrix} x_1\\x_2 \end{bmatrix} + \begin{bmatrix} 0\\ f_{fr}(x_2) \end{bmatrix} = \begin{bmatrix} 1\\0 \end{bmatrix}\cos\Omega t.
# ```

# The steps followed in this file are the same as in [Example C1](@ref ex_c1), except for the fact that we will now specify that the nonlinearity type is `:Hyst`:
# 1. [First](@ref exc2_setup) we specify the system under study in terms of its "MCK" matrices and the nonlinearities it contains.
# 2. [Next](@ref exc2_hbsetup) we setup the parameters for Harmonic Balance and setup the harmonic excitation vector.
# 3. And that's it! We're ready to [compute the forced responses](@ref exc2_cont) and [visualize the results](@ref exc2_plot).

# ## Preamble: Load Packages
using GLMakie
using LinearAlgebra
using SparseArrays
using ForwardDiff
using NonlinearSolve
using DSP

using Revise #src
using juliajim.HARMONIC
using juliajim.CONTINUATION
using juliajim.MDOFUTILS

# ## [System Setup](@id exc2_setup)
# We provide the MCK matrices just as before.

M = collect(1.0I(2));
K = [2. -1.;-1. 2.];
C = 0.01*M+0.001*K;

mdl = MDOFGEN(M, C, K);

# For the nonlinearity setup the signature of the nonlinear force is different (read the documentation of [`ADDNL`](@ref juliajim.MDOFUTILS.ADDNL)), now conducting incremental evaluation and making explicit reference to the displacement and force at the previous step(s). This is added to the `MDOFGEN` model the same way as before, using the [`ADDNL`](@ref juliajim.MDOFUTILS.ADDNL) routine.

## Nonlinearity
kt = 1.0;
fs = 1.0;
fnl = (t,u,up,fp)-> if all(abs.(fp+kt*(u-up)).<fs)
    return fp+kt*(u-up), kt*ones(size(u)), -kt*ones(size(u)), ones(size(u)); else
        return fs*sign.(fp+kt*(u-up)), zeros(size(u)), zeros(size(u)), zeros(size(u));
end
L = [0.0 1.0];

mdl = ADDNL(mdl, :Hyst, fnl, L);

# ## [Setup HB](@id exc2_hbsetup)
# We setup the HB, test nonlinear force evaluation through [`NLEVAL!`](@ref juliajim.MDOFUTILS.NLEVAL!), and test the HB residue.

h = 0:5;
N = 256;
t = (0:N-1)*2π/N;

Nhc = sum(all(h.==0, dims=2) + 2any(h.!=0, dims=2));

_, _, zinds, rinds, iinds = HINDS(mdl.Ndofs, h)
Fl = zeros(Nhc*mdl.Ndofs);
Fl[rinds[1]] = 1.0;

Wst = 0.6;
E, dEdw = HARMONICSTIFFNESS(mdl.M, mdl.C, mdl.K, Wst, h);

Uw0 = [E\Fl; Wst];

# ### [Evaluate the Nonlinear Forces](@id exc2_nleval)
FNL = zeros(mdl.Ndofs*Nhc);
dFNLdU = zeros(mdl.Ndofs*Nhc, mdl.Ndofs*Nhc);
dFNLdw = zeros(mdl.Ndofs*Nhc);
NLEVAL!(2Uw0, mdl, h, N; FNL=FNL, dFNLdU=dFNLdU, dFNLdw=dFNLdw)

# ### [Test HB Residue](@id exc2_testhb)
R = zeros(mdl.Ndofs*Nhc);
dRdU = zeros(mdl.Ndofs*Nhc, mdl.Ndofs*Nhc);
dRdw = zeros(mdl.Ndofs*Nhc);
HBRESFUN!(Uw0, mdl, Fl, h, N; R=R, dRdU=dRdU, dRdw=dRdw)

Famp = 1.0;

fun = NonlinearFunction((r,u,p)->HBRESFUN!([u;p], mdl, Famp*Fl, h, N; R=r),
    jac=(J,u,p)->HBRESFUN!([u;p], mdl, Famp*Fl, h, N; dRdU=J),
    paramjac=(Jp,u,p)->HBRESFUN!([u;p], mdl, Famp*Fl, h, N; dRdw=Jp));

prob = NonlinearProblem(fun, Uw0[1:end-1], Uw0[end]);
sol = solve(prob, show_trace=Val(true));

# ## [Continuation](@id exc2_cont)
# Also like before, we do the continuation to compute the full forced response curve.

# Note that the code in the last two sections ([Setup HB](@ref exc2_hbsetup) and [Continuation](@ref exc2_cont)) are identical to their counterparts in [Example C1](@ref ex_c1). This is what [`juliajim`](@ref) is all about - the analysis is abstracted enough that the code should look nearly the same except for the setup of the nonlinearity.

Om0 = 0.1;
Om1 = 3;
dOm = 0.2;
cpars = (parm=:arclength, nmax=1000, Dsc=:none);

sols, its, dss, xis, Dsc = CONTINUATE(Uw0[1:end-1], fun, [Om0, Om1], dOm; cpars...);

uh = zeros(Complex, maximum(h)+1, length(sols), 2);
for i in 1:2
    uh[h.+1, :, i] = hcat([[s.up[zinds[i:2:end]];
                            s.up[rinds[i:2:end]]+1im*s.up[iinds[i:2:end]]]
                          for s in sols]...);
end
Oms = [s.up[end] for s in sols];

# ## [Plotting](@id exc2_plot)
# Just like before, we visualize the results. Notice, just like in [`Example B2`](@ref ex_b2) that the non-smooth behavior of the nonlinearity is very striking.

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
    fig
end #src
