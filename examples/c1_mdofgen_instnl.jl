# ```@meta
# CurrentModule = juliajim
# ```

# # [Example C1: Forced Response of a 2DoF Oscillator with Instantaneous Nonlinearity](@id ex_c1)

# This example is designed to expose the core features of the [`MDOFUTILS`](@ref) module. As already noted in [Example b2](@ref exb2_outro), much of the setup while doing these analyses are very repetitive and common. [`MDOFUTILS`](@ref) abstracts much of this so the user can focus on the dynamics.

# The dynamical system that will be studied here is:
# ```math
# \underbrace{\begin{bmatrix} 1&0\\0&1 \end{bmatrix}}_{\mx{M}} \begin{bmatrix} \ddot{x_1}\\ \ddot{x_2} \end{bmatrix} +
# \left( 0.01 \mx{M} + 0.001 \mx{K}\right) \begin{bmatrix} \dot{x_1}\\ \dot{x_2} \end{bmatrix} + \underbrace{\begin{bmatrix} 2&-1\\-1&2 \end{bmatrix}}_{\mx{K}} \begin{bmatrix} x_1\\x_2 \end{bmatrix} + \begin{bmatrix} 0\\ \beta x_2^3 \end{bmatrix} = \begin{bmatrix} 1\\0 \end{bmatrix}\cos\Omega t.
# ```

# The steps followed in this file are:
# 1. [First](@ref exc1_setup) we specify the system under study in terms of its "MCK" matrices and the nonlinearities it contains.
# 2. [Next](@ref exc1_hbsetup) we setup the parameters for Harmonic Balance and setup the harmonic excitation vector.
# 3. And that's it! We're ready to [compute the forced responses](@ref exc1_cont) and [visualize the results](@ref exc1_plot).

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

# ## [System Setup](@id exc1_setup)
# We specify the linear "MCK" matrices of the system here. Also check out the documentation for the [`MDOFGEN`](@ref juliajim.MDOFUTILS.MDOFGEN) struct. 

M = I(2);
K = [2 -1;-1 2];
C = 0.01*M+0.001*K;

mdl = MDOFGEN(M, C, K);

# Next we specify the nonlinearity using a function handle. This function handle must take 3 input arguments (time, displacement, velocity) and output 3 results (force, displacement derivative, velocity derivative). Although the nonlinearity is a simple spring here, this same routine also supports more complicated nonlinearities involving multiple degrees of freedom. In such cases, `u,ud` must be an \(N\times d\) matrix where \(N\) is the number of time samples and \(d\) is the number of nonlinear degrees-of-freedom. It should also return $d$ forces (for adjoint forcing).

# Along with this, we also specify a "selector matrix" \(\mx{L}\) such that the nonlinear DoFs are \(\mx{L}\vc{u}\) when \(\vc{u}\) is the vector of all the DoFs.

## Nonlinearity
β = 0.1;
fnl = (t,u,ud)->return β.*u.^3, 3β.*u.^2, zeros(size(u));;
L = [0.0 1.0];

# We "add" the nonlinearity to the `MDOFGEN` object (`mdl` here) using the [`ADDNL`](@ref MDOFUTILS.ADDNL) routine. The second argument here specifies that the nonlinearity can be evaluated "instantaneously", like a cubic spring (UNlike a Jenkins element, as in [`Example b2`](@ref ex_b2). [`ADDNL`](@ref MDOFUTILS.ADDNL) also supports non-self adjoint forcing - check its documentation for this.
mdl = ADDNL(mdl, :Inst, fnl, L);

# ## [Setup HB](@id exc1_hbsetup)
# Next we setup the harmonics necessary and the number of time samples to use per period for the Alternating Frequency Time (AFT) step. We setup the excitation vector in the frequency domain by using the indices returned by [`HINDS`](@ref juliajim.HARMONIC.HINDS) - we set `Fl` to be non-zero at the first cosine harmonic index, which will correspond to the first harmonic of the first DoF. 

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

# ### [Evaluate the Nonlinear Forces](@id exc1_nleval)
# The [`NLEVAL!`](@ref juliajim.MDOFUTILS.NLEVAL!) routine is defined in [`MDOFUTILS`](@ref). This uses the appropriate methodology to evaluate each nonlinearity that has been added to the given `MDOFGEN` model. The function is written in the in-place evaluation for efficiency. The routine directly just returns the nonlinear force harmonics (`FNL` here).
FNL = zeros(mdl.Ndofs*Nhc);
dFNLdU = zeros(mdl.Ndofs*Nhc, mdl.Ndofs*Nhc);
dFNLdw = zeros(mdl.Ndofs*Nhc);
NLEVAL!(Uw0, mdl, h, N; FNL=FNL, dFNLdU=dFNLdU, dFNLdw=dFNLdw);
FNL

# ### [Test HB Residue](@id exc1_testhb)
# [`HBRESFUN!`](@ref juliajim.MDOFUTILS.HBRESFUN!) implements the Harmonic Balance forced response residue. This can be used in the same way that the `RESFUN!` functions were used in [`Example b1`](@ref ex_b1) and [`Example b2`](@ref ex_b2) for forced response evaluation.

R = zeros(mdl.Ndofs*Nhc);
dRdU = zeros(mdl.Ndofs*Nhc, mdl.Ndofs*Nhc);
dRdw = zeros(mdl.Ndofs*Nhc);
HBRESFUN!(Uw0, mdl, Fl, h, N; R=R, dRdU=dRdU, dRdw=dRdw)

fun = NonlinearFunction((r,u,p)->HBRESFUN!([u;p], mdl, Fl, h, N; R=r),
    jac=(J,u,p)->HBRESFUN!([u;p], mdl, Fl, h, N; dRdU=J),
    paramjac=(Jp,u,p)->HBRESFUN!([u;p], mdl, Fl, h, N; dRdw=Jp));

prob = NonlinearProblem(fun, Uw0[1:end-1], Uw0[end]);
sol = solve(prob, show_trace=Val(true));

# ## [Continuation](@id exc1_cont)
# [`HBRESFUN!`](@ref juliajim.MDOFUTILS.HBRESFUN!) can also be used in tandem with [`CONTINUATE`](@ref juliajim.CONTINUATION.CONTINUATE) for getting the full forced response. 

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

# ## [Plotting](@id exc1_plot)
# Now we plot the response harmonics. It can be seen that two primary and 3 secondary (super harmonic) resonances have been picked up. The nonlinearity level is quite strong, as evidenced by the non-trivial shapes of the response curves.

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
