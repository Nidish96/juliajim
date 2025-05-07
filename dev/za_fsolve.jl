# * Preamble
using NonlinearSolve
using LinearAlgebra
using UnPack

# * Learning
# ** Introductory Case
f = (u,p) -> u.*u .- p;
u0 = [1.,1.];
p = 2.;
prob = NonlinearProblem(f, u0, p);
sol = solve(prob, show_trace=Val(true));

sol.stats

# ** Also Supply Jacobian
function f!(du, u, p)
    du .= u.*u.-p;
    return nothing;
end
function J!(J,u,p)
    J[:,:] = 2Diagonal(u);
    return nothing;
end
fun = NonlinearFunction(f!; jac=J!);
prob2 = NonlinearProblem(fun, u0, p);
sol = solve(prob2, show_trace=Val(true));

sol.stats

# ** Also Supply Jacobian together
function fJ!(du, J, u, p)
    if du !== nothing
        du .= u.*u.-p;
    end
    if J !== nothing
        J[:,:] = 2Diagonal(u);
    end
    return nothing;
end
fun = NonlinearFunction((du,u,p)->fJ!(du, nothing, u, p); jac=(J,u,p)->fJ!(nothing, J, u, p));
prob3 = NonlinearProblem(fun, u0, p);
sol = solve(prob2, show_trace=Val(true));

sol.stats

# * Duffing
function duffresfun!(du, J, u, p, JOm=nothing)
    Om, z0, w0, al, F = p;
    A0 = u[1];
    b0 = u[2];

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

funduff = NonlinearFunction((du,u,p)->duffresfun!(du,nothing,u,p), jac=(J,u,p)->duffresfun!(nothing,J,u,p));
pars = (Om = 0.1, z0 = 0.1e-2, w0 = 2., al = 0.01, F = 0.1);

Ab0 = [0.01, 0.];
dprob = NonlinearProblem(funduff, Ab0, collect(pars));

sol = solve(dprob, show_trace=Val(true));

