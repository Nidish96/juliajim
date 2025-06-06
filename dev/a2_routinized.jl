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

sol = solve(dprob, show_trace=Val(true), store_trace=Val(true));

# * Routinized Version
Om0 = 0.85pars.w0;
Om1 = 1.15pars.w0;

Alin = pars.F/(pars.w0^2-Om0^2+2im*pars.z0*pars.w0*Om0);
Ab0 = [abs(Alin), angle(Alin)];

funduff = NonlinearFunction((du,u,Om)->duffresfun!([u;Om],pars;du=du));
funduffJ = NonlinearFunction((du,u,Om)->duffresfun!([u;Om],pars;du=du),
                            jac=(J,u,Om)->duffresfun!([u;Om],pars;J=J),
                            paramjac=(JOm,u,Om)->duffresfun!([u;Om],pars;Jp=JOm));

dOm = 0.05pars.w0;
sols, its, dss, xis, Dsc = CONTINUATE(Ab0, funduff, [Om0, Om1], dOm);

# Plots in 2D
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

# * Analytical Solution

Om1_f(A0) = -(sqrt(pars.F^2-4A0^2*pars.w0^4*pars.z0^2)/(2A0*pars.w0))+pars.w0+(3A0^2*pars.al)/8pars.w0;
Om2_f(A0) = sqrt(pars.F^2-4A0^2*pars.w0^4*pars.z0^2)/(2A0*pars.w0)+pars.w0+(3A0^2*pars.al)/8pars.w0

impfun(A0,Om) = -(A0^2*pars.w0^2*pars.z0^2)-(A0*(pars.w0-Om)+(3A0^3*pars.al)/8pars.w0)^2+pars.F^2/4pars.w0^2;

isapprox.([impfun.(s.up[1], s.up[3]) for s in sols], 0.0; atol=eps())

aex = extrema([s.up[1] for s in sols]);
Na = 100;
as = exp10.(range(log10(aex[1]/10), log10(10aex[2]), Na))
b1 = [Om1_f.(Complex.(as)) as];
b2 = [Om2_f.(Complex.(as)) as];
b1 = real(b1[isreal.(b1[:,1]), :]);
b2 = real(b2[isreal.(b2[:,1]), :]);

ansol = [b1; b2[end:-1:1,:]];

fsz = 18;
fig2 = Figure(fontsize=fsz);
if !isdefined(Main, :scr2) && isdefined(Main, :GLMakie)
   scr2 = GLMakie.Screen();
end

ax = Axis(fig2[1, 1], xlabel="Excitation Frequency (rad/s)", ylabel="Response (m)");
lines!(ax, ansol[:,1], ansol[:,2], color=:blue, label="Analytical")
scatterlines!(ax, [s.up[end] for s in sols], [s.up[1] for s in sols], color=:red, label="Numerical");
xlims!(ax, Om0, Om1)

axislegend(ax);

if isdefined(Main, :GLMakie)
   display(scr2, fig2);
else
    fig2   
end
   
