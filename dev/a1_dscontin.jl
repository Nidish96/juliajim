using NonlinearSolve
using NonlinearSolve: ForwardDiff
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

Alin = pars.F/(pars.w0^2-Om^2+2im*pars.z0*pars.w0*Om);
Ab0 = [abs(Alin), angle(Alin)];
dprob = NonlinearProblem(funduff, Ab0, Om);

sol = solve(dprob, show_trace=Val(true));

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
function EXTRESFUN!(up, fun, sol0, ds;
                    parm=:riks, Dsc::Vector{Float64}, dup=nothing, Jf=nothing)

    if dup !== nothing
        du = @view dup[1:end-1]

        fun.f(du, up[1:end-1], up[end])
    end
    if Jf !== nothing
        J = @view Jf[1:end-1, 1:end-1];
        Jp = @view Jf[1:end-1, end];

        if fun.jac !== nothing
            fun.jac(J, up[1:end-1], up[end]);
        else
            Jt = @view Jf[1:end-1, :];
            R = similar(up[1:end-1]);
            ForwardDiff.jacobian!(Jt, (du,up) -> fun.f(du,up[1:end-1], up[end]), R, up);
        end
        if fun.paramjac !== nothing
            fun.paramjac(Jp, up[1:end-1], up[end]);
        end
    end

    # Extension part
    if dup !== nothing
        if parm==:riks
            tgt = normalize(sol0.dupds./Dsc);
    	    dup[end] = tgt'*((up-sol0.up)./Dsc)-ds;
        elseif parm==:arclength
            dup[end] = norm((up-sol0.up)./Dsc)^2-ds^2;
        end
    end
    if Jf !== nothing
        if parm==:riks
            tgt = normalize(sol0.dupds./Dsc);
    	    Jf[end, :] = (tgt./Dsc)';
        elseif parm==:arclength
            Jf[end, :] = 2((up-sol0.up)./Dsc.^2)';
        end
    end
    return nothing;
    
end

# ** Continuation
funduff = NonlinearFunction((du,u,Om)->duffresfun!([u;Om],pars;du=du),
                            jac=(J,u,Om)->duffresfun!([u;Om],pars;J=J),
                            paramjac=(JOm,u,Om)->duffresfun!([u;Om],pars;Jp=JOm));

# funduff = NonlinearFunction((du,u,Om)->duffresfun!([u;Om],pars;du=du));

# Continuation Input
Om0 = 0.85pars.w0;
Om1 = 1.15pars.w0;

# Initial Guess
Alin = pars.F/(pars.w0^2-Om0^2+2im*pars.z0*pars.w0*Om0);
Ab0 = [abs(Alin), angle(Alin)];

# Continuation Input Parameters
dw0 = 0.05;  # This is in the units of Omega
itopt = :auto;
parm = :arclength;
Dsc = :auto;
nmax = 1000;
DynScale = true;

# Step Length Adaptation Parameters
dwmin = dw0/5;
dwmax = dw0*5;
nxi = 0.5;
ximin = 0.5;
ximax = 2.0;

# Continuation Routine ########################################################
Rf = similar(Ab0, length(Ab0)+1);
Jf = similar(Ab0, length(Ab0)+1,length(Ab0)+1);

R = @view Rf[1:end-1];
J = @view Jf[1:end-1,1:end-1];
Jp = @view Jf[1:end-1,end];

# Initialize Storers
sols = myNLSoln[];
dss = Float64[];
its = Int[];
xis = Float64[];

# Converge to first point
prob0 = NonlinearProblem(funduff, Ab0, Om0);
solp0 = solve(prob0, show_trace=Val(true));
if funduff.jac !== nothing
    funduff.jac(J, solp0.u, Om0);
else
    Jt = @view Jf[1:end-1, :];
    ForwardDiff.jacobian!(Jt, (R,up)->funduff.f(R,up[1:end-1],up[end]), R,
                          [solp0.u; Om0])
end
if funduff.paramjac !== nothing
    funduff.paramjac(Jp, solp0.u, Om0);
end

push!(sols, myNLSoln([solp0.u;Om0]; J=copy(J), Jp=copy(Jp)));
push!(its, solp0.stats.nsteps)

if Dsc==:auto
    Dsc = max.(abs.(sols[end].up), 100eps());
end


# Recontextualize ds0 (such that first step is as requested)
ds0 = dw0/Dsc[end]normalize(sols[end].dupds./Dsc)[end];
push!(dss,  ds0);
push!(xis,  1.0);
dsmin = dwmin/dw0*ds0;
dsmax = dwmax/dw0*ds0;

# Setup Extended Problem
exfun = NonlinearFunction((du,up,p)->EXTRESFUN!(up, funduff, p[1], p[2];
                                                parm=parm, Dsc=p[3], dup=du),
                          jac=(J,up,p)->EXTRESFUN!(up, funduff, p[1], p[2];
                                                   parm=parm, Dsc=p[3], Jf=J));


tgt = Dsc.*normalize(sols[end].dupds./Dsc);
up0 = sols[end].up + dss[end]tgt;
prob = NonlinearProblem(exfun, up0, (sols[end], dss[end], Dsc));

if itopt == :auto # Automatic itopt
    solp = solve(prob);
    itopt = solp.stats.nsteps;
end

while sols[end].up[end]*sign(Om1-Om0)<Om1*sign(Om1-Om0) && length(sols)<=nmax
    # Tangent Predictor
    global tgt = Dsc.*normalize(sols[end].dupds./Dsc);
    global up0 = sols[end].up + dss[end]tgt;

    # Constrained Corrector
    global prob = remake(prob; u0=up0, p=(sols[end], dss[end], Dsc));

    global solp = solve(prob, store_trace=Val(true));
    prob.f.jac(Jf, solp.u, prob.p);
    push!(its, solp.stats.nsteps)

    # Push to sols
    push!(sols, myNLSoln(solp.u; J=Jf[1:end-1,1:end-1], Jp=Jf[1:end-1,end]));
    # Fix New Tangent Sign
    sols[end].dupds .*= sign((tgt./Dsc)'normalize(sols[end].dupds./Dsc));

    # Print out Message
    println(@sprintf("%d. %.2f with %.4f converged in %d iterations.",
                     length(sols), sols[end].up[end], dss[end], its[end]))

    
    # Step Length Adaptation
    push!(xis, clamp((itopt/its[end])^nxi, ximin, ximax));
    push!(dss, clamp(xis[end]dss[end], dsmin, dsmax));

    if DynScale # Rescale
    	Dsc_ = copy(Dsc);

        rat = clamp.(abs.(sols[end].up)./Dsc_, 0.5, 2.0);

        dAl = 1.0;    
        global Dsc = (1-dAl)Dsc_ + dAl*rat.*Dsc_;
    end
end

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
