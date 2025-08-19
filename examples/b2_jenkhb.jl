# ```@meta
# CurrentModule = juliajim
# ```

# # [Example b2: Numerical Continuation on a Frictional Oscillator](@id ex_b2)

# The scope as well as the approach of this example are very similar to [Example b1](@ref ex_b1), with a key difference: this example involves a hysteretic nonlinearity. This poses a very important difference in how the nonlinear force and its harmonics are computed.

# The general setup for the problem is similar:
# ```math
# \ddot{x} + 2\zeta_0\omega_0\dot{x} + \omega_0^2 x + f_{nl} = F\cos\Omega t.
# ```
# The nonlinearity $f_{nl}$ is a hysteretically saturated linear spring with low amplitude stiffness $k_t$ and saturation limit ("slip" force) $f_s$.

# The steps that will be followed in this file are exactly the same as in [Example b1](@ref ex_b1):
# 1. [First](@ref exb2_res) we define a "harmonic residue" function for the Duffing oscillator in a way that returns the residue, its Jacobian with respect to the vector of harmonics (see [`Example a`](@ref ex_a)), and also with respect to the excitation frequency.
# 2. [Next](@ref exb2_setup) we define the parameters for the problem and setup the necessary variables for the Harmonic Balance (and AFT).
# 3. [Then](@ref exb2_cont) we conduct the actual continuation after setting up parameters for this.
# 4. [Finally](@ref exb2_plot) the results are plotted in terms of the different harmonics present.

# ## Preamble: Load Packages
using GLMakie
using LaTeXStrings
using LinearAlgebra
using NonlinearSolve
using ForwardDiff
using DSP
using Infiltrator #src

using Revise #src
using juliajim.HARMONIC
using juliajim.CONTINUATION

# ## [Define Residue Function](@id exb2_res)
# The structure of the residue function is identical to that in [Example 1](@ref exb1_res). A key differerence is in the computation of the nonlinear force.

# Unlike in the previous case, the restoring force of a hysteretic element can not be estimated at any instant of time independent of its history. In other words, the restoring force of a hysteretic element can be different based on its loading history. For steady state calculations, we start by assuming initially that the hysteretic element acts as a linear spring (see `ft=kt*ut` below) and then "march forward" in time. At each instant the incremental form of the friction law is used to check if saturation is violated and the restoring force is corrected appropriately. The force "shakes down" after repeating this for 2 or more cycles of the period. We use this steady-state force for the Fourier transforms involve in Harmonic Balance.
function RESFUN!(Uw, Fl, pars, h, Nt; R=nothing, dRdU=nothing, dRdw=nothing)
    (; z0, w0, kt, fs, F) = pars;

    Om = Uw[end];
    Nhc = sum((h.==0)+2(h.!=0));

    ## Linear Portion
    E, dEdw = HARMONICSTIFFNESS(1.0, 2z0*w0, w0^2, Om, h);

    ## AFT For nonlinear force
    ut  = AFT(Uw[1:end-1], h, Nt, :f2t);
        
    ## Construct Residue
    if !(R === nothing && dRdU === nothing)
        ft = kt*ut;
        if !(dRdU === nothing)
            cst = AFT(eltype(Uw).(I(Nhc)), h, Nt, :f2t);
            dfdat = kt.*cst;
        end
        
        for _ in 1:2            
            for (ti, tim1) in zip(1:Nt, circshift(1:Nt,1))
                fsp = kt*(ut[ti]-ut[tim1]) + ft[tim1];  # stick prediction
                ft[ti] = clamp(fsp, -fs, fs);
                if !(dRdU === nothing)
                    if (abs(fsp)<fs)
                        dfdat[ti, :] = kt.*(cst[ti,:]-cst[tim1,:]) + dfdat[tim1,:];
                    else
                        dfdat[ti, :] .= 0.0;
                    end
                end
            end
        end
        if !(R === nothing)
            Fnl = AFT(ft, h, Nt, :t2f);        
            R[:] = E*Uw[1:end-1] + Fnl - Fl*F;
        end
        if !(dRdU === nothing)
            Jnl    = AFT(dfdat, h, Nt, :t2f);
            dRdU[:, :] = E + Jnl;
        end
    end
    if !(dRdw === nothing)
        dRdw[:] = dEdw*Uw[1:end-1];
    end
    return nothing;
end

# ## [Setup](@id exb2_setup)
# Now we setup the specific problem by assigning numerical values to the parameters. Odd harmonics are chosen for the Harmonic Balance since we already have some intuition on the nonlinearity (it, being an "odd function" nonlinearity). As before, this is followed by setting up a `NonlinearFunction` object that will be used for the continuation.

pars = (z0 = 0.5e-2, w0 = 2., kt = 5.0, fs=1.0, F = 0.1);

## h = (0:5); # Also possible
h = 1:2:5;
Om = 0.1;

Nhc = sum((h.==0)+2(h.!=0));
Nt = 2^9;

Fl = zeros(Nhc);
_,_,zinds,rinds,iinds = HINDS(1, h)
Fl[rinds[1]] = 1.0;

fun = NonlinearFunction((r,u,p)->RESFUN!([u;p],Fl,pars,h,Nt;R=r),
    jac=(J,u,p)->RESFUN!([u;p],Fl,pars,h,Nt;dRdU=J),
    paramjac=(Jp,u,p)->RESFUN!([u;p],Fl,pars,h,Nt;dRdw=Jp));

# ### Check by solving for a Single point
# We do the single point check to ensure that everything works - no errors in the residue, the Jacobian is correct, etc. This should also converge within 3-4 iterations. 

E = zeros(Nhc, Nhc);
HARMONICSTIFFNESS!(E, nothing, 1.0, 2pars.z0*pars.w0, pars.w0^2+pars.kt, Om, h);
U0 = E\ (Fl*pars.F);

prob = NonlinearProblem(fun, U0, Om);
sol = solve(prob, show_trace=Val(true));

# ## [Continuation](@id exb2_cont)
# Finally we're ready to do the continuation so we setup the system as before: initialize continuation parameters and setup the initial guess. Then we invoke [`CONTINUATE`](@ref juliajim.CONTINUATION.CONTINUATE) to compute the solutions and then arrange the harmonics in complex notation for plotting.

Om0 = 0.02pars.w0;
Om1 = 2pars.w0;
dOm = 0.04pars.w0;
cpars = (parm=:riks, nmax=300, Dsc=:auto, minDsc=1e-2);

HARMONICSTIFFNESS!(E, nothing, 1.0, 2pars.z0*pars.w0, pars.w0^2+pars.kt, Om0, h);
U0 = E\ (Fl*pars.F);

sols, its, dss, xis, Dsc = CONTINUATE(U0, fun, [Om0, Om1], dOm; cpars...);

uh = zeros(Complex, maximum(h)+1, length(sols));
uh[h.+1, :] = hcat([[s.up[zinds]; s.up[rinds]+1im*s.up[iinds]] for s in sols]...);
Oms = [s.up[end] for s in sols];

# ## [Plotting](@id exb2_plot)
# Here again, we visualize the harmonic contents. Unlike the previous example, a very clear non-smooth behavior can be noticed here. It is only near the resonance, when amplitude is large, where the higher harmonics show up. At other frequencies, the higher harmonics are all zero. 

his = [1, 3, 5];

fsz = 24;
fig = Figure(fontsize=fsz, size=(1000, 600));
if !isdefined(Main, :scr) && isdefined(Main, :GLMakie) #src
    scr = GLMakie.Screen(); #src
end #src

ax1s = [];
ax2s = [];
for i in eachindex(his[his.<=maximum(h)])
    ax = Axis(fig[1, i],
        ylabel=L"$H_%$(his[i])$ Response (m)", yscale=log10);
    scatterlines!(ax, Oms, abs.(uh[his[i].+1, :]));
    push!(ax1s, ax)

    ax = Axis(fig[2, i], xlabel=L"Excitation Frequency $\Omega$",
        ylabel=L"$H_%$(his[i])$ Phase (rad)");
    scatterlines!(ax, Oms, unwrap(angle.(uh[his[i].+1, :])));
    push!(ax2s, ax)
end
linkxaxes!(vcat(ax1s, ax2s)...)

if isdefined(Main, :GLMakie) #src
    display(scr, fig); #src
else #src
    fig
end #src

# ## [Outro](@id exb2_outro)

# This is yet another example of using the utility functions in [`HARMONIC`](@ref) for conducting harmonic balance. Looking at [`Example b1`](@ref ex_b1) and this, it must be clear that a lot of the "setup" that has to be done for any given system is very similar. So it makes sense to have a unified interface which constructs the different residue functions for a given system with its nonlinearity. This is precisely what the [`MDOFUTILS`](@ref) module sets out to do. 
