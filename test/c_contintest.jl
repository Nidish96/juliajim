using NonlinearSolve
using NonlinearSolve: ForwardDiff
using LinearAlgebra
using Printf
using Test

using Revise
includet("../src/CONTINUATION.jl")

# * Duffing MMS System Setup
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

# * Continuation
Om0 = 0.85pars.w0;
Om1 = 1.15pars.w0;

Alin = pars.F/(pars.w0^2-Om0^2+2im*pars.z0*pars.w0*Om0);
Ab0 = [abs(Alin), angle(Alin)];

dOm = 0.05pars.w0;
sols, its, dss, xis, Dsc = CONTINUATE(Ab0, funduff, [Om0, Om1], dOm, verbosity=0);

# * Analytical Solution

Om1_f(A0) = -(sqrt(pars.F^2-4A0^2*pars.w0^4*pars.z0^2)/(2A0*pars.w0))+pars.w0+(3A0^2*pars.al)/8pars.w0;
Om2_f(A0) = sqrt(pars.F^2-4A0^2*pars.w0^4*pars.z0^2)/(2A0*pars.w0)+pars.w0+(3A0^2*pars.al)/8pars.w0;

impfun(A0,Om) = -(A0^2*pars.w0^2*pars.z0^2)-(A0*(pars.w0-Om)+(3A0^3*pars.al)/8pars.w0)^2+pars.F^2/4pars.w0^2;

@testset "MMS Forced Response Continuation" begin
    @test sols[end].up[end]*sign(Om1-Om0)>=Om1*sign(Om1-Om0)
    @test all(isapprox.([impfun.(s.up[1], s.up[3]) for s in sols], 0.0; atol=eps()))
end
