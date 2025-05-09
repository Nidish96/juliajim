# * Preamble
using NonlinearSolve

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

pars = (z0 = 0.1e-2, w0 = 2., al = 0.01, F = 0.1);
Om = 0.1;
funduff = NonlinearFunction((du,u,Om)->duffresfun!([u;Om],pars;du=du),
                            jac=(J,u,Om)->duffresfun!([u;Om],pars;J=J));

Ab0 = [1., 0.];
dprob = NonlinearProblem(funduff, Ab0, Om);

sol = solve(dprob, show_trace=Val(true));

# * Iterator Interface
nlcache = init(dprob; show_trace=Val(true));

# step!(nlcache);

# * My Attempt
# ** Define a struct with constructor
struct myNLSoln
    u::Union{Nothing, Vector{Float64}}
    J::Union{Nothing, Matrix{Float64}}
    Jp::Union{Nothing, Vector{Float64}}
    dudp::Union{Nothing, Vector{Float64}}
end
function myNLSoln(u=nothing; J=nothing, Jp=nothing)
    if J === nothing || Jp === nothing
        dudp = nothing;
    else
        dudp = -J\Jp;
    end
    return myNLSoln(u, J, Jp, dudp);
end
function Base.show(io::IO, p::myNLSoln)
    print(io, " u = $(p.u)\n dudp = $(p.dudp)\n J = $(p.J)\n Jp = $(p.Jp)")
end
function Base.:-(v1::myNLSoln, v2::myNLSoln)
    return myNLSoln(v1.u-v2.u, v1.J-v2.J, v1.Jp-v2.Jp, v1.dudp-v2.dudp);
end

# ** Try out continuation

# Temporary Variables
J = zeros(2,2); JOm = zeros(2);

# First 2 Omega values
Om0 = 0.1; Om1 = 0.15;

# Compute first two points
prob = remake(dprob; p=Om0);
sol0 = solve(prob);
duffresfun!([sol0.u;Om0], pars;J=J,JOm=JOm);
sol0 = myNLSoln([sol0.u;Om0]; J=copy(J), Jp=copy(JOm));
prob = remake(prob; u0=sol0.u[1:2], p=Om1);
sol1 = solve(prob);
duffresfun!([sol1.u;Om1], pars;J=J, JOm=JOm);
sol1 = myNLSoln([sol1.u;Om1]; J=copy(J), Jp=copy(JOm));

dsol1 = sol1-sol0;
