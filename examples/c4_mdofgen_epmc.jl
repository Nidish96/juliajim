# ```@meta
# CurrentModule = juliajim
# ```

# # [Example C4: Nonlinear Normal Modes of a 2DoF Hysteretic Oscillator using the Extended Periodic Motion Concept](@id ex_c4)

# !!! note "The Extended Periodic Motion Concept"
#     Nonlinear normal modes is a very succinct methodology for describing the dynamics of nonlinear systems. It leads to descriptions of the dynamics in terms of amplitude-dependent natural frequency, daming factor, and mode-shapes for each mode. The [Extended Periodic Motion Concept](https://doi.org/10.1016/j.compstruc.2015.03.008) is a very popular methodology for capturing these numerically. The idea here is simple - it iteratively looks for what would be a linear self excitation (in terms of a negative viscous factor) that must be supplied to a system so that it has a periodic solution as an invariant set in its unforced state.

#     Mathematically the system we will be solving is:
#     ```math
#     \mx{M} \ddot{\vc{u}} + \mx{C} \dot{\vc{u}} + \mx{K} \vc{u} + \vc{f_{nl}} = \xi \mx{M} \dot{\vc{u}}
#     ```
#     where $\xi$ controls the strength of the self-excitation. Since the system is always solved by the trivial solution, an amplitude constraint is often imposed. In order to simplify calculations later on, the first harmonic amplitude of the degrees of freedom is constrained as:
#     ```math
#     \vc{U}_{c1}^T \mx{M} \vc{U}_{c1} + \vc{U}_{s1}^T \mx{M} \vc{U}_{s1} = a^2
#     ```
#     where $\vc{U}_{c1} = \frac{\Omega}{\pi} \int_0^{2\pi/\Omega} \vc{u} \cos\Omega t dt$ and $\vc{U}_{s1} = \frac{\Omega}{\pi} \int_0^{2\pi/\Omega} \vc{u} \sin\Omega t dt$. \(a\) is the user-specified amplitude, which serves as the continuation parameter in this context. In addition to this, there is still a phase indeterminacy present in the system. So we choose an arbitrary DoF (which is expected to be non-zero in the mode under consideration), and set its cosine harmonic to be zero (other constraints are also possible). With this, we will have a square system (same number of unknowns as equations) and it can be solved to obtain the amplitude-dependent nonlinear modal characteristics of the system.

# The [`EPMCRESFUN!`](@ref juliajim.MDOFUTILS.EPMCRESFUN!) routine in [`MDOFUTILS`](@ref) implements the residue function for EPMC that can be used in tandem with [`CONTINUATION`](@ref) to obtain nonlinear modal characteristics. This example will illustrate this on the frictional 2DoF oscillator under consideration in Examples [C2](@ref ex_c2) and [C3](@ref ex_c3). The steps followed are:
# 1. [First](@ref exc4_setup) the system is setup (along with the nonlinearity, see descriptions in [Example C2](@ref exc2_setup).
# 2. [Next](@ref exc4_linmodal) we conduct linear modal analysis on the linearized system. This is important to set the initial guess for the EPMC run next.
# 3. [Now](@ref exc4_hbsetup) we setup the Harmonic Balance parameters and then [test out the EPMC residue routine](@ref exc4_epmcres) using the linearized modes estimated above.
# 4. [Finally](@ref exc4_cont) we conduct the contuation over the modal amplitude \(a\) and then [visualize](@ref exc4_plot) the results.

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

# ## [System Setup](@id exc4_setup)
# MCK matrices and nonlinearity setup: identical to [Example C2](@ref exc2_setup).

M = collect(1.0I(2));
K = [2. -1.;-1. 2.];
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

# ## [Linear Modal Analysis](@id exc4_linmodal)
# Remember that in the limit of small amplitude, the Jenkins element acts like a linear spring. So we add this additional stiffness to the linearized stiffness matrix and conduct the modal analysis (eigen decomposition). The mode shapes, natural frequencies, and damping are saved since these will be used next as the initial guess for EPMC.

Knl0 = L'kt*L;
K0 = mdl.K+Knl0;
D, V = eigen(K0, mdl.M)
W0s = sqrt.(D);
Xis = diag(V'mdl.C*V)

mi = 1;  # Mode of interest

# ## [Setup HB](@id exc4_hbsetup)
# Setting up parameters for Harmonic Balance: identical to [Example C2](@ref exc2_setup).

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

# ## [Test EPMC Residue](@id exc4_epmcres)
# The ordering of unknowns that will be sent to the [`EPMCRESFUN!`](@ref juliajim.MDOFUTILS.EPMCRESFUN!) routine is \(\begin{bmatrix} \vc{u}^T&\omega&\xi&\log_{10}(a) \end{bmatrix}^T\), where the last parameter is the log of the amplitude (base 10). The log-scaling is done on the amplitude to help with numerical conditioning issues.

# In the following we setup everything and test it to ensure that the residue can be used to converge to the solution.

R = zeros(mdl.Ndofs*Nhc+2);
dRdU = zeros(mdl.Ndofs*Nhc+2, mdl.Ndofs*Nhc+2);
dRda = zeros(mdl.Ndofs*Nhc+2);

EPMCRESFUN!([U0; W0s[mi]; Xis[mi]; -2], mdl, Fl, h, N; R=R, dRdUwx=dRdU, dRda=dRda)

A0 = -2.0;

U0 = zeros(mdl.Ndofs*Nhc);
U0[iinds[1:mdl.Ndofs]] = V[:, mi];

fun = NonlinearFunction((r,uwx,p)->EPMCRESFUN!([uwx;p], mdl, Fl, h, N; R=r),
    jac=(J,uwx,p)->EPMCRESFUN!([uwx;p], mdl, Fl, h, N; dRdUwx=J),
    paramjac=(Jp,uwx,p)->EPMCRESFUN!([uwx;p], mdl, Fl, h, N; dRda=Jp));

prob = NonlinearProblem(fun, [U0;W0s[mi];Xis[mi]], A0);
sol = solve(prob, show_trace=Val(true), maxiters=100);

# ## [Continuation](@id exc4_cont)
# Now we are ready to conduct continuation to obtain the amplitude dependent resonance backbone for the system. The invocation of [`CONTINUATE`](@ref juliajim.CONTINUATION.CONTINUATE) is exactly as encountered in the earlier examples.

A0 = -1.;
A1 = 3.;
da = 0.05;
cpars = (nmax=1000, Dsc=:auto, angopt=deg2rad(30));

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

# ## [Plotting](@id exc4_plot)
# Finally we visualize the results in terms of harmonics as well as in terms of natural frequency and (equivalent) damping factor characteristic curves. 

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
    if i==1
        axislegend(ax)
    end
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












