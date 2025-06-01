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
ξ = 1e0;
pars = (z0 = 0.5e-2, w0 = 2., al = 0.1*ξ^2, F = 0.1/ξ);
Om = 0.1;
funduff = NonlinearFunction((du,up,Om)->duffresfun!([up;Om],pars;du=du),
                            jac=(J,up,Om)->duffresfun!([up;Om],pars;J=J));

Ab0 = [1., 0.];
dprob = NonlinearProblem(funduff, Ab0, Om);

sol = solve(dprob, show_trace=Val(true), store_trace=Val(true));

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

# * Two-Step Starting
# ** Extended Residue Function 2
function EXTRESFUN2!(up, fun, sol0s; deflate=false, parm=:arclength, xi=1.0, dup=nothing, Jf=nothing)
    
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
    if deflate
        dist = norm(up-sol0s[end-1].up)^2;

        if dup === nothing
            Rp = zeros(length(up));
            R = @view Rp[1:end-1]
            fun.f(R, up[1:end-1], up[end]);
        else
            R = du;
        end

        if Jf !== nothing
            J[:,:] = J/dist - 2R*((up-sol0s[end-1].up)[1:end-1])'/dist^2;
            Jp[:] = Jp/dist - 2R*((up-sol0s[end-1].up)[end])'/dist^2;
        end
        if dup !== nothing
    	    du[:] = du/dist;
        end
    end

    # Arc Parameterization
    ur = 1.0./abs.(sols[end].up);
    ur[@. !isfinite(ur)] .= 1.0;

    # ur = ones(length(up));
    if dup !== nothing
    	if parm==:riks
            # ds = xi*(ur.*sol0s[end-1].dupds)'*(ur.*(sol0s[end].up-sol0s[end-1].up))/((ur.*sol0s[end-1].dupds)'sol0s[end-1].dupds);
            ds = xi*(ur.*sol0s[end-1].dupds)'*(ur.*(sol0s[end].up-sol0s[end-1].up));            
            dup[end] = (ur.*sol0s[end].dupds)'*(ur.*(up-sol0s[end].up)) - ds;
        elseif parm==:arclength
            ds = xi*norm(ur.*(sol0s[end].up-sol0s[end-1].up));
            dup[end] = (ur.*(up-sol0s[end].up))'*(ur.*(up-sol0s[end].up))-ds^2;
        else
            error("Unknown parm")
        end
    end
    if Jf !== nothing
        if parm==:riks
            Jf[end, :] = (ur.*sol0s[end].dupds)'*diagm(ur);
        elseif parm==:arclength
            Jf[end, :] = 2(ur.*(up-sol0s[end].up))'*diagm(ur);
        else
            error("Unknown parm")
        end
    end
end

# ** Continuation
funduff = NonlinearFunction((du,up,Om)->duffresfun!([up;Om],pars;du=du),
                            jac=(J,up,Om)->duffresfun!([up;Om],pars;J=J),
                            paramjac=(JOm,up,Om)->duffresfun!([up;Om],pars;JOm=JOm));

# Temporary Variables
J = zeros(2,2); JOm = zeros(2);
Jf = zeros(3,3);

# Continuation Parameters
parm = :riks;
Nopt = 5;
nmax = 1000;

# TODO: adaptation - need to consider alternatives.
nxi = 0.2;
xirng = (0.5, 2.0);

# Storage Struct vecs
sols = myNLSoln[];
xis = Float64[];
its = Int[];
ers = Vector{Float64}[];

# Starting Frequency Values
Om0s = [0.85, 0.9]pars.w0;
Om1 = 1.15*pars.w0;

Ab0 = [1., 0.];
for Om in Om0s
    Alin = pars.F/(pars.w0^2-Om^2+2im*pars.z0*pars.w0*Om);
    Ab0 = [abs(Alin), angle(Alin)];
        
    prob = NonlinearProblem(funduff, Ab0, Om);
    solp0 = solve(prob);
    if solp0.u[1]<0
        solp0.u[1] = -solp0.u[1];
        solp0.u[2] = mod2pi(π+solp0.u[2])-2π;
    end
    solp0.u[2] = mod2pi(solp0.u[2])-2π;
    funduff.jac(J, solp0.u, Om)
    funduff.paramjac(JOm, solp0.u, Om)

    sol0 = myNLSoln([solp0.u;Om]; J=copy(J), Jp=copy(JOm));
    sol0.dupds .*= sign((Om1-Om)*sol0.dupds[end]);
    push!(sols, sol0);
    push!(xis, 1.0);
    push!(its, solp0.stats.nsolve);
    push!(ers, [solp0.u;Om]-[Ab0;Om]);
end

# Setup Extended Problem
exfun = NonlinearFunction((du,up,p)->EXTRESFUN2!(up, funduff, p[1];
                                                 xi=p[2], parm=parm, dup=du),
                          jac=(J,up,p)->EXTRESFUN2!(up, funduff, p[1];
                                                    xi=p[2], parm=parm, Jf=J));

ur = 1.0./abs.(sols[end].up);
ur[@. !isfinite(ur)] .= 1.0;

while sols[end].up[end]<Om1 && length(sols)<=nmax
    # Secant Predictor
    sec_ = (sols[end].up-sols[end-1].up);
    sgn = sign((ur.*sols[end-1].dupds)'*(ur.*sec_))
    up0_s = sols[end].up + xis[end]sgn.*sec_;

    # tangent predictor
    ds = (ur.*sols[end-1].dupds)'*(ur.*(sols[end].up-sols[end-1].up))/((ur.*sols[end-1].dupds)'sols[end-1].dupds);
    up0_t = sols[end].up + xis[end]ds*sols[end].dupds;

    up0 = up0_t;
    
    # Correct
    prob = NonlinearProblem(exfun, up0, (sols[end-1:end], xis[end]));
    solp = solve(prob, store_trace=Val(true));
    push!(ers, solp.u-up0);
    push!(its, solp.stats.nsolve);

    # Get gradients
    exfun.jac(Jf,solp.u,(sols[end-1:end], xis[end]));
    push!(sols, myNLSoln(solp.u; J=Jf[1:end-1,1:end-1], Jp=Jf[1:end-1,end]));
    # sols[end].dupds .*= sign(sols[end-1].dupds'sols[end].dupds);

    ur0 = copy(ur);
    
    ur = 1.0./abs.(sols[end].up);
    ur[@. !isfinite(ur)] .= 1.0;
    sols[end].dupds .*= sign((ur.*sols[end-1].dupds)'*(ur.*sols[end].dupds));

    println(@sprintf("%d. %.2f converged in %d iterations.",
                     length(sols), sols[end].up[end], its[end]))

    # Step length adaptation
    xi = clamp((Nopt/solp.stats.nf)^nxi, xirng[1], xirng[2]);
    push!(xis, xi);
end

# Plot in 2D
fsz = 18;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr) && isdefined(Main, :GLMakie)
   scr = GLMakie.Screen();
end

ax = Axis(fig[1:3, 1], ylabel="Response (m)");
scatterlines!(ax, [s.up[end] for s in sols], [s.up[1] for s in sols]);

ax = Axis(fig[4:6, 1], xlabel="Excitation Frequency (rad/s)",
          ylabel="Response Phase (rad)");
scatterlines!(ax, [s.up[end] for s in sols], [s.up[2] for s in sols]);

ax = Axis(fig[1:2, 2], ylabel="Response Tangent (m)");
scatterlines!(ax, [s.up[end] for s in sols], [s.dupds[1] for s in sols]);

ax = Axis(fig[3:4, 2], ylabel="Phase Tangent (rad)");
scatterlines!(ax, [s.up[end] for s in sols], [s.dupds[2] for s in sols]);

ax = Axis(fig[5:6, 2], xlabel="Excitation Frequency (rad/s)",
          ylabel="Frequency Tangent (rad/s)");
scatterlines!(ax, [s.up[end] for s in sols], [s.dupds[3] for s in sols]);

ax = Axis(fig[1:2, 3], ylabel="Response Secant (m)");
scatterlines!(ax, [s.up[end] for s in sols[1:end-1]],
              [s.up[1] for s in (sols[2:end].-sols[1:end-1])]);

ax = Axis(fig[3:4, 3], ylabel="Phase Secant (rad)");
scatterlines!(ax, [s.up[end] for s in sols[1:end-1]],
              [s.up[2] for s in (sols[2:end].-sols[1:end-1])]);

ax = Axis(fig[5:6, 3], xlabel="Excitation Frequency (rad/s)",
          ylabel="Frequency Secant (rad/s)");
scatterlines!(ax, [s.up[end] for s in sols[1:end-1]],
              [s.up[3] for s in (sols[2:end].-sols[1:end-1])]);

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
