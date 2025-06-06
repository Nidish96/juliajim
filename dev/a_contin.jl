using NonlinearSolve
using LinearAlgebra
using Printf
using GLMakie
using LaTeXStrings

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
ξ = 1e0;
pars = (z0 = 0.5e-2, w0 = 2., al = 0.1*ξ^2, F = 0.1/ξ);
Om = 0.1;
funduff = NonlinearFunction((du,up,Om)->duffresfun!([up;Om],pars;du=du),
                            jac=(J,up,Om)->duffresfun!([up;Om],pars;J=J));

Ab0 = [1., 0.];
dprob = NonlinearProblem(funduff, Ab0, Om);

sol = solve(dprob, show_trace=Val(true));

# ** Iterator Interface
nlcache = init(dprob; show_trace=Val(true));

# step!(nlcache);

# * Define a struct to store the solution point
struct myNLSoln
    up::Union{Nothing, Vector{Float64}}
    J::Union{Nothing, Matrix{Float64}}
    Jp::Union{Nothing, Vector{Float64}}
    dupds::Union{Nothing, Vector{Float64}}  # HAS to be a vector only!
end
function myNLSoln(up=nothing; J=nothing, Jp=nothing)
    if J === nothing || Jp === nothing
        dupds = nothing;
    else
        # dupds = -J\Jp;  # Naive, requires IFT to hold.
        dupds = nullspace([J Jp])[:,1]; # General, but have to handle bifurcations better!
        # TODO: Need to fix sign of tangent!
    end
    return myNLSoln(up, J, Jp, dupds);
end
function Base.show(io::IO, p::myNLSoln)
    print(io, " up = $(p.up)\n dupds = $(p.dupds)\n J = $(p.J)\n Jp = $(p.Jp)")
end
function Base.:-(v1::myNLSoln, v2::myNLSoln)
    return myNLSoln(v1.up-v2.up, v1.J-v2.J, v1.Jp-v2.Jp, v1.dupds-v2.dupds);
end

# * Single Step Starting
# ** Define the Extended Residue Function
function EXTRESFUN!(up, fun, sol0, ds; parm=:riks, Dsc=1, dup=nothing, Jf=nothing)

    if dup !== nothing
        du = @view dup[1:end-1]

        fun.f(du, up[1:end-1], up[end])
    end
    if Jf !== nothing
        J = @view Jf[1:end-1, 1:end-1];
        Jp = @view Jf[1:end-1, end];

        fun.jac(J, up[1:end-1], up[end]);
        fun.paramjac(Jp, up[1:end-1], up[end]);
    end

    # Extension part
    if dup !== nothing
        if parm==:riks
    	    dup[end] = (sol0.dupds)'diagm(Dsc.^2)*(up-sol0.up)-ds;
        elseif parm==:arclength
            dup[end] = (up-sol0.up)'diagm(Dsc.^2)*(up-sol0.up)-ds^2;
        end
    end
    if Jf !== nothing
        if parm==:riks
    	    Jf[end, :] = sol0.dupds'diagm(Dsc.^2);
        elseif parm==:arclength
            Jf[end, :] = 2(up-sol0.up)'diagm(Dsc.^2);
        end
    end
    return nothing;
    
end

# R = ones(3);
# Jf = ones(3,3);
# EXTRESFUN!(sol0.up, funduff, sol0, ds; dup=R, Jf=Jf);

# ** Continuation
funduff = NonlinearFunction((du,up,Om)->duffresfun!([up;Om],pars;du=du),
                            jac=(J,up,Om)->duffresfun!([up;Om],pars;J=J),
                            paramjac=(JOm,up,Om)->duffresfun!([up;Om],pars;Jp=JOm));

# Temporary Variables
R = ones(2);
J = zeros(2,2); JOm = zeros(2);

# Continuation Parameters
ds0 = 0.01;
Nopt = 6;
dsmin = ds0/5;
dsmax = ds0*5;
nmax = 5000;
parm = :arclength; # :riks, :arclength

Om0 = 0.85*pars.w0;
Om1 = 1.15*pars.w0;

Ab0 = [1., 0.];

# Storage Struct vecs
sols = myNLSoln[];
dss = Float64[];

# Starting Point
prob0 = NonlinearProblem(funduff, Ab0, Om0);
solp0 = solve(prob0);
if solp0.u[1]<0
    solp0.u[1] = -solp0.u[1];
    solp0.u[2] = mod2pi(π+solp0.u[2])-2π;
end
solp0.u[2] = mod2pi(solp0.u[2])-2π;
duffresfun!([solp0.u;Om0], pars; J=J,Jp=JOm);
sol0 = myNLSoln([solp0.u;Om0]; J=copy(J), Jp=copy(JOm));
sol0.dupds .*= sign(Om1-Om0);
push!(sols, sol0);
push!(dss, ds0);

# Setup Problem
exfun= NonlinearFunction((du,up,p)->EXTRESFUN!(up, funduff, p[1],p[2];
                                               Dsc=p[3], parm=parm, dup=du),
                         jac=(J,up,p)->EXTRESFUN!(up, funduff, p[1],p[2];
                                                  Dsc=p[3], parm=parm, Jf=J));
Ab10 = sols[end].up + dss[end]*sols[end].dupds;
prob1 = NonlinearProblem(exfun, Ab10, (sols[end], dss[end], ones(3)));

while sols[end].up[end]<Om1 && length(sols)<=nmax

    # Unknowns Scaling
    # Approach 0: No scaling
    Dsc = ones(3);
    # Dsc[1] = ξ;
    # Dsc /= norm(Dsc.*sols[end].dupds);

    # Approach 1: Setting magnitudes to unity
    # Dsc = 1.0./abs.(sols[end].up);
    # Dsc[@. !isfinite(Dsc)] .= 1.0;
    # Dsc ./= norm(Dsc.*sols[end].dupds);

    # # Approach 2: Based on secant magnitude 
    # if length(sols)==1
    #     dup = Dsc.*sols[end].dupds;
    # else
    #     dup = Dsc.*(sols[end].up-sols[end-1].up);
    #     dup /= norm(dup);
    # end
    # dup ./= abs.(sols[end].up);
    # sdup = sort(abs.(dup), rev=true);
    # knee_id = argmax(abs.(diff(sdup)))+1;
    # thresh = sdup[knee_id];

    # # Dsc = ones(length(sols[end].up));
    # Dsc[abs.(dup).<thresh] .= 0.0;
    # Dsc /= norm(Dsc.*sols[end].dupds);
    
    # Make a Step
    Ab10 = sols[end].up + dss[end]*sols[end].dupds;
    prob1_ = remake(prob1; u0=Ab10, p=(sols[end], dss[end], Dsc) );

    solp1 = solve(prob1_);
    duffresfun!(solp1.u, pars;J=J, Jp=JOm);
    push!(sols, myNLSoln(solp1.u; J=copy(J), Jp=copy(JOm)));
    sols[end].dupds .*= sign(sols[end-1].dupds'diagm(Dsc)*sols[end].dupds);

    println(@sprintf("%d. %.2f with %.4f converged in %d iterations.",
                     length(sols), sols[end].up[end], dss[end], solp1.stats.nsolve))

    # Step length adaptation
    xi = clamp(Nopt/solp1.stats.nf, 0.5, 2.0);
    push!(dss, clamp(dss[end]*xi, dsmin,dsmax));
end

for s in sols
    if s.up[1]<0
        s.up[1] = -s.up[1];
        s.up[2] = mod2pi(π+s.up[2])-2π;
    end
end

# *** Plots in 2D
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

# *** Plot in 3D
fsz = 18;
fig2 = Figure(fontsize=fsz);
if !isdefined(Main, :scr2) && isdefined(Main, :GLMakie)
   scr2 = GLMakie.Screen();
end

ax = Axis3(fig2[1, 1], xlabel=L"Frequency $\Omega$", ylabel=L"Amplitude $A_0$",
           zlabel=L"Phase $\beta_0$");
scatterlines!(ax, [s.up[3] for s in sols], [s.up[1] for s in sols],
              [s.up[2] for s in sols]);

if isdefined(Main, :GLMakie)
   display(scr2, fig2);
else
    fig2
end
