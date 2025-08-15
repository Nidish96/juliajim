using GLMakie
using LaTeXStrings
using LinearAlgebra
using NonlinearSolve
using ForwardDiff
using DSP
using Infiltrator

using Revise
includet("../src/CONTINUATION.jl")
includet("../src/HARMONIC.jl")

# * Residue Function
function RESFUN!(Uw, Fl, pars, h, Nt; R=nothing, dRdU=nothing, dRdw=nothing)
    (; z0, w0, kt, fs, F) = pars;

    Om = Uw[end];
    Nhc = sum((h.==0)+2(h.!=0));

    # Linear Portion
    E, dEdw = HARMONICSTIFFNESS(1.0, 2z0*w0, w0^2, Om, h);

    # AFT For nonlinear force
    ut  = AFT(Uw[1:end-1], h, Nt, :f2t);
        
    # Construct Residue
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

# * Setup
pars = (z0 = 0.5e-2, w0 = 2., kt = 5.0, fs=1.0, F = 0.1);

# h = (0:5);
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

E = zeros(Nhc, Nhc);
HARMONICSTIFFNESS!(E, nothing, 1.0, 2pars.z0*pars.w0, pars.w0^2+pars.kt, Om, h);
U0 = E\ (Fl*pars.F);

prob = NonlinearProblem(fun, U0, Om);
sol = solve(prob, show_trace=Val(true))

R = zeros(Nhc);
Jf = zeros(Nhc, Nhc+1);
J = @view Jf[:, 1:Nhc];
Jp = @view Jf[:, end];

# * Continuation
Om0 = 0.02pars.w0;
Om1 = 2pars.w0;

HARMONICSTIFFNESS!(E, nothing, 1.0, 2pars.z0*pars.w0, pars.w0^2+pars.kt, Om0, h);
U0 = E\ (Fl*pars.F);

fun = NonlinearFunction((r,u,p)->RESFUN!([u;p],Fl,pars,h,Nt;R=r),
                        jac=(J,u,p)->RESFUN!([u;p],Fl,pars,h,Nt;dRdU=J),
                        paramjac=(Jp,u,p)->RESFUN!([u;p],Fl,pars,h,Nt;dRdw=Jp));

dOm = 0.04pars.w0;
# dOm = 0.1;
cpars = (parm=:riks, nmax=300, minDsc=1e-2);
sols, its, dss, xis, Dsc = CONTINUATE(U0, fun, [Om0, Om1], dOm; cpars...);

uh = zeros(Complex, maximum(h)+1, length(sols));
uh[h.+1, :] = hcat([[s.up[zinds]; s.up[rinds]+1im*s.up[iinds]] for s in sols]...);
Oms = [s.up[end] for s in sols];

# Plot

his = [1, 3, 5];

fsz = 24;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr) && isdefined(Main, :GLMakie)
   scr = GLMakie.Screen();
end

for i in eachindex(his[his.<=maximum(h)])
    ax = Axis(fig[1, i], xlabel=L"Excitation Frequency $\Omega$",
              ylabel=L"$H_%$(his[i])$ Response (m)", yscale=log10);
    scatterlines!(ax, Oms, abs.(uh[his[i].+1, :]));

    ax = Axis(fig[2, i], xlabel=L"Excitation Frequency $\Omega$",
              ylabel=L"$H_%$(his[i])$ Phase (rad)");
    scatterlines!(ax, Oms, unwrap(angle.(uh[his[i].+1, :])));
end

if isdefined(Main, :GLMakie)
   display(scr, fig);
else
    fig
end
