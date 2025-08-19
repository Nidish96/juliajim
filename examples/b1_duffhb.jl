# ```@meta
# CurrentModule = juliajim
# ```

# # [Example B1: Numerical Continuation on a Duffing Oscillator](@id ex_b1)

# This example exposes the core functionality of [`juliajim.CONTINUATION`](@ref) through a very minimal example. The governing equations for the oscillator are taken as
# ```math
# \ddot{x} + 2\zeta_0\omega_0 \dot{x} + \omega_0^2 x + \alpha x^3 = F\cos\Omega t.
# ```

# The steps that will be followed in this file are as follows:
# 1. [First](@ref exb1_res) we define a "harmonic residue" function for the Duffing oscillator in a way that returns the residue, its Jacobian with respect to the vector of harmonics (see [`Example a`](@ref ex_a)), and also with respect to the excitation frequency.
# 2. [Next](@ref exb1_setup) we define the parameters for the problem and setup the necessary variables for the Harmonic Balance (and AFT).
# 3. [Then](@ref exb1_cont) we conduct the actual continuation after setting up parameters for this.
# 4. [Finally](@ref exb1_plot) the results are plotted in terms of the different harmonics present.

# ## Preamble: Load packages
using GLMakie
using LaTeXStrings
using LinearAlgebra
using NonlinearSolve
using ForwardDiff
using DSP
using Infiltrator #src

using Revise #src
using juliajim.HARMONIC
using juliajim.CONTINUATION

# ## [Define Residue Function](@id exb1_res)

# We setup the Harmonic Balance residue function such that it takes as input the list of harmonics and the excitation frequency as a single vector. This is convenient for bordered continuation strategies. 

# The computation of the residue itself is as follows:
# 1. Estimate the linear harmonic stiffness
# 2. Transform the given harmonic coefficients into the time domain (for one cycle) using [`AFT`](@ref juliajim.HARMONIC.AFT), and use that to compute the nonlinearity in one oscillation cycle.
# 3. Transform the computed force into the Frequency domain by computing its harmonics using [`AFT`](@ref juliajim.HARMONIC.AFT).
# 4. Assemble the overall harmonic balance equation.

# Jacobian computations also follow similar steps. The function is implemented "in-place" in order to be efficient - the Jacobian(s) are computed only when requested.   
function RESFUN!(Uw, Fl, pars, h, Nt; R=nothing, dRdU=nothing, dRdw=nothing)
    (; z0, w0, al, F) = pars;

    Om = Uw[end];
    Nhc = sum((h.==0)+2(h.!=0));

    ## Linear Portion
    E, dEdw = HARMONICSTIFFNESS(1.0, 2z0*w0, w0^2, Om, h);

    ## AFT For nonlinear force
    ut  = AFT(Uw[1:end-1], h, Nt, :f2t);
        
    ## Construct Residue
    if !(R === nothing)
        ft  = al*ut.^3;
        Fnl = AFT(ft, h, Nt, :t2f);
        
        R[:] = E*Uw[1:end-1] + Fnl - Fl*F;
    end
    if !(dRdU === nothing)
        cst = AFT(Float64.(I(Nhc)), h, Nt, :f2t);
        dfdat = (3al*ut.^2) .* cst;
        Jnl    = AFT(dfdat, h, Nt, :t2f);
        
        dRdU[:, :] = E + Jnl;
    end
    if !(dRdw === nothing)
        dRdw[:] = dEdw*Uw[1:end-1];
    end
    return nothing;
end

# ## [Setup](@id exb1_setup)

# We now setup the parameters of the system and forcing vector for Harmonic balance. 

  ξ = 1e0;  # Scaling fudge factor - change to test the robustness of the continuation
pars = (z0 = 0.5e-2, w0 = 2., al = 0.1*ξ^2, F = 0.1/ξ);

h = (0:5);  # Choose list of harmonics
## h = 1:2:5;  # Also possible
Om = 0.1;

Nhc = sum((h.==0)+2(h.!=0));
Nt = 2^9;

Fl = zeros(Nhc);
_,_,zinds,rinds,iinds = HINDS(1, h)
Fl[rinds[1]] = 1.0;
# Now `Fl` is the harmonic forcing vector. 
# Next, we define a `NonlinearFunction` object (from [`NonlinearSolve.jl`]) that will be used by `juliajim` for the continuation. Note that for this, the unknowns are the harmonics and the parameter (has to be scalar) is the excitation frequency, with which the continuation procedure will be carried out.

# We also use the same `RESFUN!` function defined above to specify the `jac` (jacobian w.r.t. the unknowns) and paramjac (jacobian w.r.t. the frequency). If these are left empty, `ForwardDiff` will be used to estimate them. 

fun = NonlinearFunction((r,u,p)->RESFUN!([u;p],Fl,pars,h,Nt;R=r),
                        jac=(J,u,p)->RESFUN!([u;p],Fl,pars,h,Nt;dRdU=J),
                        paramjac=(Jp,u,p)->RESFUN!([u;p],Fl,pars,h,Nt;dRdw=Jp));

# Now we are ready!

# ### Check by solving for a Single point

# We set the initial guess as the solution of the linearized oscillator and solve the problem using `NonlinearSolve.jl`. The following should work properly and converge within 3-4 iterations.

E = zeros(Nhc, Nhc);
HARMONICSTIFFNESS!(E, nothing, 1.0, 2pars.z0*pars.w0, pars.w0^2, Om, h);
U0 = E\ (Fl*pars.F);

prob = NonlinearProblem(fun, U0, Om);
sol = solve(prob, show_trace=Val(true));

# Inspecting the solution would show that the solution has only odd-harmonic content, as expected from a simple system with a cubic spring:
sol.u	   

# ## [Continuation](@id exb1_cont)
# Now we are ready to conduct the actual continuation. First things first, we setup the starting and ending frequencies (`Om0`, `Om1`) and the starting step length `dOm`. **In `juliajim` the units of the step length are always in the same units as the continuation parameter, the excitation frequency in this case.**

Om0 = 2pars.w0;
Om1 = 0.02pars.w0;
dOm = 0.1pars.w0;
cpars = (parm=:arclength, nmax=2000);

# We also set two parameters relevant for continuation (see the documentation of [`CONTINUATE`](@ref juliajim.CONTINUATION.CONTINUATE) for an exhaustive list of parameters with descriptions). Here we are specifying that we would use the arclength parameterization and will do the computation for a maximum of 2000 points.
# We provide the initial guess the same as before (of the linearized problem) and can reuse the same `NonlinearFunction` object from before.

HARMONICSTIFFNESS!(E, nothing, 1.0, 2pars.z0*pars.w0, pars.w0^2, Om0, h);
U0 = E\ (Fl*pars.F);

sols, its, dss, xis, Dsc = CONTINUATE(U0, fun, [Om0, Om1], dOm; cpars...);

# CONTINUATE returns, firstly, a vector of [`myNLSoln`](@ref juliajim.CONTINUATION.myNLSoln) objects. These store the solution vector `up`, local tanget `dupds`, and any Jacobians if required. 

# In order to help with plotting and other postprocessing, it is sometimes helpful to express the harmonics in the complex notation. So we do that here before moving on to plotting.
uh = zeros(Complex, maximum(h)+1, length(sols));
uh[h.+1, :] = hcat([[s.up[zinds]; s.up[rinds]+1im*s.up[iinds]] for s in sols]...);
Oms = [s.up[end] for s in sols];

# ## [Plotting](@id exb1_plot)

# Now we choose a few harmonics and plot them using [`GLMakie`](https://docs.makie.org/stable/explanations/backends/glmakie.html) (I highly recommend Makie for interativity as well as for its rich feature-set). The output should show the primary resonance (visible in all harmonics) as well as secondary superharmonic resonances (visible in the higher harmonics).

his = [1, 3, 5];

fsz = 24;
fig = Figure(fontsize=fsz, size = (1000, 600));
if !isdefined(Main, :scr) && isdefined(Main, :GLMakie) #src
   scr = GLMakie.Screen(); #src
end #src

ax1s = [];
ax2s = [];
for i in eachindex(his[his.<=maximum(h)])
    ax = Axis(fig[1, i], ylabel=L"$H_%$(his[i])$ Response (m)", yscale=log10);
    scatterlines!(ax, Oms, abs.(uh[his[i].+1, :]));
    push!(ax1s, ax);

    ax = Axis(fig[2, i], xlabel=L"Excitation Frequency $\Omega$",
              ylabel=L"$H_%$(his[i])$ Phase (rad)");
    scatterlines!(ax, Oms, unwrap(angle.(uh[his[i].+1, :])));
    push!(ax2s, ax);
end
linkxaxes!(vcat(ax1s, ax2s)...)

if isdefined(Main, :GLMakie) #src
   display(scr, fig); #src
else #src
    display(fig) #src
    fig
end #src

# ## Outro

# This file shows how to use the routines in [`HARMONIC`](@ref) to construct a residue and then use this function with [`CONTINUATION`](@ref) to conduct numerical continuation to compute its forced response.
