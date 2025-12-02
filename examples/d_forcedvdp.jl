# ```@meta
# CurrentModule = juliajim
# ```

# # [Example D: Forced Response of a Van der Pol Oscillator](@id ex_d)

# This example is designed to demonstrate stability analysis and branch switching across a (secondary) Hopf bifurcation that occurs in a forced Van der Pol oscillator.

# The dynamical system that will be studied here is:
# ```math
# \ddot{x} - c \dot{x} + k x + \eta x^2 \dot{x} = F\cos\Omega t.
# ```
# This is an SDoF nonlinear oscillator with negative viscous damping. The nonlinearity provides the dissipation that saturates the response. Classically known as the Van der Pol oscillator, this is a useful system for learning about limit cycles. In the considered version, we also have an external periodic excitation.

# The system will have a stable periodic response close to resonance, which will lose its stability as we move away from it.

# ## Preamble: Load Packages
using GLMakie
using LinearAlgebra
using SparseArrays
using NonlinearSolve
using DSP

using Revise
using juliajim.HARMONIC
using juliajim.CONTINUATION
using juliajim.MDOFUTILS

# ## [System Setup](@id exd_setup)
# Here we setup the system with parameters.

pars = (c=0.02, k=4., F=1., eta=0.1);

mdl = MDOFGEN(1.0, -pars.c, pars.k);
fnl = (t,u,ud) -> return pars.eta.*u.^2 .*ud, 2pars.eta.*u.*ud, pars.eta.*u.^2;
mdl = ADDNL(mdl, :Inst, fnl, 1.0);

# ## [Setup HB](@id exd_hbsetup)

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

# ## [Continuation](@id exd_cont)
# Just like in other examples we call the [`CONTINUATE`](@ref juliajim.CONTINUATION.CONTINUATE) utility to obtain the periodic forced response.
Om0 = 0.1;
Om1 = 4.0;
dOm = 0.2;
cpars = (parm=:arclength, nmax=100, save_jacs=true);

sols, its, dss, xis, Dsc = CONTINUATE(Uw0[1:end-1], fun, [Om0, Om1], dOm; cpars...);

# Obtain Harmonics
uh = zeros(Complex, maximum(h)+1, length(sols));
uh[[inds0; hinds], :] = hcat([[up[zinds];up[rinds,:]+1im*up[iinds,:]] for up in sols.up]...);
Oms = [up[end] for up in sols.up];

# ### Stability AnalysisWe now use an averaging formulation to obtain the stability coefficients based on just the first harmonics. This will work if the response is dominantly single harmonic.

E0, _ = HARMONICSTIFFNESS(0., mdl.M, 0, 1, h[h.!=0]);
stab = zeros(length(Oms));
for (i,(J,Om)) in enumerate(zip(sols.J, Oms))
    # evs = eigvals(J[2:end,2:end], -Om*collect(E0));  # Multiharmonic
    evs = eigvals(J[2:3,2:3], -Om*collect(E0[1:2,1:2]));
    stab[i] = sum(real(evs).>=0);
end

# The variable `stab` will store the number of unstable eigenvalues that have been detected. One may interpret these as the Floquet exponents of the system. 

# ##  [Plotting](@id exd_plot)
set_theme!(theme_latexfonts())
fsz = 24;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr) && Makie.current_backend()==GLMakie
   scr = GLMakie.Screen();
end

ax = Axis(fig[1, 1], xlabel=L"Excitation Frequency $\Omega$", ylabel="Response");
scatterlines!(ax, Oms./(stab.==0), [norm(u) for u in eachcol(uh)], label="Stable")
scatterlines!(ax, Oms./(stab.==2), [norm(u) for u in eachcol(uh)], label="Unstable")

axislegend(ax)
xlims!(Om0, Om1)
if Makie.current_backend()==GLMakie
   display(scr, fig);
else
   display(fig)
end
