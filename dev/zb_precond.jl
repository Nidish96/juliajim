using NLsolve, NonlinearSolve
using LinearAlgebra
using Printf
using GLMakie
using LaTeXStrings

# * Duffing
function duffresfun!(uOm, p; du=nothing, J=nothing, JOm=nothing)
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
    if JOm !== nothing
    	JOm[:] = [0., -1.0];
    end
    return nothing;
end

# pars = (z0 = 0.5e-2, w0 = 2., al = 0.1, F = 0.1);
ξ = 1e-6;
pars = (z0 = 0.5e-2, w0 = 2., al = 0.1*ξ^2, F = 0.1/ξ);
Om = 2.0;
funduff = NonlinearFunction((du,up,Om)->duffresfun!([up;Om],pars;du=du),
                            jac=(J,up,Om)->duffresfun!([up;Om],pars;J=J));

Ab0 = [1., 0.];
Alin = pars.F/(pars.w0^2-Om^2+2im*pars.z0*pars.w0*Om);
Ab0 = [abs(Alin), angle(Alin)];
dprob = NonlinearProblem(funduff, Ab0, Om);

sol = solve(dprob, show_trace=Val(true), store_trace=Val(true));

# * Preconditioning
Dsc = [1/ξ, 1.0];
# Dsc = [1.0, 1.0];

alg = NLsolveJL(method=:newton, linsolve=(x,A,b)->begin x[:]=Dsc.*((Dsc.*A.*Dsc')\(Dsc.*b)); return x; end);

dprob = NonlinearProblem(funduff, Ab0, Om);
sol = solve(dprob, alg, show_trace=Val(true), store_trace=Val(true))
