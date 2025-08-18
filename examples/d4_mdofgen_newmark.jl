# # [Example d4](@id ex_d4)
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
Fv = [1.0, 0.0];

mdl = ADDNL(mdl, :Hyst, fnl, L);

# ## Nonlinear Force Evaluation
Fnl, dFnldU, dFnldUd, S = NLFORCE(mdl, 0, zeros(2), zeros(2))

# ## Transient March
U0 = zeros(2);
Ud0 = zeros(2);
fsamp = 10;
Om = 0.6;
Famp = 4;
Fex = t-> Famp*Fv*cos(Om*t);
Tmax = 2π/Om*100;

t = 0:1/fsamp:Tmax;

U, Ud, Udd, S = NEWMARKMARCH(mdl, 0, Tmax, 1/fsamp, U0, Ud0, Fex, verbosity=100)

# ## Plots
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
