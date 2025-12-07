# ```@meta
# CurrentModule = juliajim
# ```

# # [Example C5: Transient Response of a 2DoF Hysteretic Oscillator using A Newmark Scheme](@id ex_c5)

# In this example we show how the [`MDOFUTILS`](@ref) suite can also be used for conducting transient analysis. A natural question may be why not just use SciML's [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/). The answer lies in hysteretic nonlinearities implemented in the incremental formalism. It is not trivial to use existing methods for integrating these systems without additional complexity (computational as well as conceptual).

# However, classical second order methods (central difference, linear acceleration, Newmark-\(\beta\), Wilson-\(\theta\), HHT-\(\alpha\), etc.) **naturally lend themselves to hysterestic systems**! This is because these integrators require the difference in force between the previous and the next time instants, which is the most convenient way of specifying an incremental hysteretic nonlinearity.

# Here we analyze the same 2DoF oscillator as in the last few examples:
# ```math
# \underbrace{\begin{bmatrix} 1&0\\0&1 \end{bmatrix}}_{\mx{M}} \begin{bmatrix} \ddot{x_1}\\ \ddot{x_2} \end{bmatrix} +
# \left( 0.01 \mx{M} + 0.001 \mx{K}\right) \begin{bmatrix} \dot{x_1}\\ \dot{x_2} \end{bmatrix} + \underbrace{\begin{bmatrix} 2&-1\\-1&2 \end{bmatrix}}_{\mx{K}} \begin{bmatrix} x_1\\x_2 \end{bmatrix} + \begin{bmatrix} 0\\ f_{fr}(x_2) \end{bmatrix} = \begin{bmatrix} 1\\0 \end{bmatrix}\cos\Omega t.
# ```

# The steps followed in the example are:
# 1. [First](@ref exc5_setup) we setup the system and the nonlinearity (and also check the nonlinear force evaluation).
# 2. [Then](@ref exc5_trans) we setup the initial conditions and transient analysis parameters and do the actual analysis.
# 3. [Finally](@ref exc5_trans) we visualize the results.

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

# ## [System Setup](@id exc5_setup)

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
Fv = [1.0, 0.0];

mdl = ADDNL(mdl, :Hyst, fnl, L);

# ### Nonlinear Force Evaluation

# Here we test the nonlinear force evaluation using the [`NLFORCE`](@ref juliajim.MDOFUTILS.NLFORCE) routine, which is used internally to compute the nonlinear forces in time domain. This function requires, at minimum, the time and displacement-velocity states. But it also accepts the states at the previous instant for the computation (very much relevant for hysteretic nonlinearities). If these are not provided, they are all defaulted to zero. 
# Apart from the nonlinear force and its Jacobians, the routine also outputs `S`, which is a Vector as long as the number of nonlinearities present. Each element of the vector is a vector of the forces of the nonlinearity. These are relevant for hysteretic nonlinearities where explicit history dependence is very important.

Fnl, dFnldU, dFnldUd, S = NLFORCE(mdl, 0, zeros(2), zeros(2))

# ## [Transient March](@id exc5_trans)

# Now we are ready to conduct the transient simulation. We set trivial initial conditions and provide periodic excitation at a point well below the first resonance (which is at \(1\) rad/s).

U0 = zeros(2);
Ud0 = zeros(2);
fsamp = 10;
Om = 0.6;  # rad/s
Famp = 4;
Fex = t-> Famp*Fv*cos(Om*t);
Tmax = 2π/Om*100;

t = 0:1/fsamp:Tmax;

U, Ud, Udd, S = NEWMARKMARCH(mdl, 0, Tmax, 1/fsamp, U0, Ud0, Fex, verbosity=100)

# ## [Plotting](@id exc5_plot)
# Now we plot all the results. We plot them in terms of the states (displacements and velocities from `U`, `Ud`) as well as the slider internal forces from `S`. A strongly nonlinear response can be observed with the steady-state zoom-in (to the last cycle) showing highly non-smooth behavior. The displacement and velocity plots also show a significant presence of higher harmonics, underlying the same.

set_theme!(theme_latexfonts())
fsz = 20;
fig = Figure(fontsize=fsz, size=(1000, 600));
if !isdefined(Main, :scr) && Makie.current_backend()==GLMakie #src
   scr = GLMakie.Screen(); #src
end #src

ax1 = Axis(fig[1, 1:2], ylabel="Displacement (m)", title="Complete Response");
for i in 1:2
    lines!(ax1, t, U[i,:])
end
ax10 = Axis(fig[1, 3], title="Steady State");
for i in 1:2
    scatterlines!(ax10, t, U[i,:], label="x$(i)")
end
axislegend(ax10, position=:ct)
xlims!(ax10, Tmax-2π/Om, Tmax);
ylims!(ax10, -maximum(abs.(U)), maximum(abs.(U)));

ax2 = Axis(fig[2, 1:2], ylabel="Velocity (m/s)");
for i in 1:2
    lines!(ax2, t, Ud[i,:])
end
ax20 = Axis(fig[2, 3]);
for i in 1:2
    scatterlines!(ax20, t, Ud[i,:])
end
xlims!(ax20, Tmax-2π/Om, Tmax);
ylims!(ax20, -maximum(abs.(Ud)), maximum(abs.(Ud)));

ax3 = Axis(fig[3, 1:2], xlabel="Time (s)", ylabel="Friction force (N)")
lines!(ax3, t, getindex.(S, 1)[:])
ax30 = Axis(fig[3, 3], xlabel="Time (s)")
scatterlines!(ax30, t, getindex.(S, 1)[:])
xlims!(ax30, Tmax-2π/Om, Tmax);

linkxaxes!(ax1, ax2, ax3);

if Makie.current_backend()==GLMakie #src
   display(scr, fig); #src
else #src
   fig
end #src
