using NonlinearSolve
using NonlinearSolve: ForwardDiff
using LinearAlgebra
using Printf
using GLMakie
using LaTeXStrings

using Revise
includet("../src/CONTINUATION.jl")

# * Duffing
function duffresfun!(uOm, p; du=nothing, J=nothing, Jp=nothing)
    (; z0, w0, al, F) = p;
    A0 = uOm[1];
    b0 = uOm[2];
    Om = uOm[3];

    if du !== nothing
        du[:] = [-z0*w0*A0-F/2w0*sin(b0), -(Om-w0)+3al/8w0*A0^2-F/2w0/A0*cos(b0)];
    end
    if J !== nothing
        J[:,:] = [-z0*w0 -F/2w0*cos(b0); 3al/4w0*A0+F/2w0/A0^2*cos(b0) F/2w0/A0*sin(b0)];
    end
    if Jp !== nothing
    	Jp[:] = [0., -1.0];
    end
    return nothing;
end

# pars = (z0 = 0.5e-2, w0 = 2., al = 0.1, F = 0.1);
ξ = 1e3;
pars = (z0 = 0.5e-2, w0 = 2., al = 0.1*ξ^2, F = 0.1/ξ);
Om = 0.1;
funduff = NonlinearFunction((du,u,Om)->duffresfun!([u;Om],pars;du=du),
                            jac=(J,u,Om)->duffresfun!([u;Om],pars;J=J),
                            paramjac=(JOm,u,Om)->duffresfun!([u;Om],pars;Jp=JOm));

Alin = pars.F/(pars.w0^2-Om^2+2im*pars.z0*pars.w0*Om);
Ab0 = [abs(Alin), angle(Alin)];
dprob = NonlinearProblem(funduff, Ab0, Om);

sol = solve(dprob, show_trace=Val(true));

# * Routinized Version
Om0 = 0.85pars.w0;
Om1 = 1.15pars.w0;

Alin = pars.F/(pars.w0^2-Om0^2+2im*pars.z0*pars.w0*Om0);
Ab0 = [abs(Alin), angle(Alin)];

dOm = 0.05;
sols, its, dss, xis = CONTINUATE(Ab0, funduff, [Om0, Om1], dOm; parm=:arclength);

# ** Plots in 2D
fsz = 18;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr) && isdefined(Main, :GLMakie)
   scr = GLMakie.Screen();
end

ax = Axis(fig[1:3, 1], ylabel="Response (m)");
scatterlines!(ax, [s.up[end] for s in sols], [s.up[1] for s in sols]);

ax = Axis(fig[4:6, 1], xlabel="Excitation Frequency (rad/s)", ylabel="Response Phase (rad)");
scatterlines!(ax, [s.up[end] for s in sols], [s.up[2] for s in sols]);

ax = Axis(fig[1:2, 2], ylabel="Response Tangent (m)");
scatterlines!(ax, [s.up[end] for s in sols], [s.dupds[1] for s in sols]);

ax = Axis(fig[3:4, 2], ylabel="Phase Tangent (rad)");
scatterlines!(ax, [s.up[end] for s in sols], [s.dupds[2] for s in sols]);

ax = Axis(fig[5:6, 2], xlabel="Excitation Frequency (rad/s)", ylabel="Frequency Tangent (rad/s)");
scatterlines!(ax, [s.up[end] for s in sols], [s.dupds[3] for s in sols]);

ax = Axis(fig[1:2, 3], ylabel="Response Secant (m)");
scatterlines!(ax, [s.up[end] for s in sols[1:end-1]], [s.up[1] for s in (sols[2:end].-sols[1:end-1])]);

ax = Axis(fig[3:4, 3], ylabel="Phase Secant (rad)");
scatterlines!(ax, [s.up[end] for s in sols[1:end-1]], [s.up[2] for s in (sols[2:end].-sols[1:end-1])]);

ax = Axis(fig[5:6, 3], xlabel="Excitation Frequency (rad/s)", ylabel="Frequency Secant (rad/s)");
scatterlines!(ax, [s.up[end] for s in sols[1:end-1]], [s.up[3] for s in (sols[2:end].-sols[1:end-1])]);

if isdefined(Main, :GLMakie)
    display(scr, fig);
else
    fig   
end
