# ```@meta
# CurrentModule = juliajim
# ```

# # [Example C3: Response Constrained Forced Response of a 2DoF Hysteretic Oscillator](@id ex_c3)

# It is sometimes the case that it is necessary to conduct forced response simulations where the response amplitude is kept fixed. In fact, experimentally it is significantly easier to achieve response-controlled forced response measurements (in the context of, say, stepped sine testing) than force-controlled measurements. While there are several reasons for this (large displacement effects on the exciter at resonance, low SNR at anti-resonance points, etc.), we will not go into this here.

# [`MDOFUTILS`](@ref) provides [`HBRESFUN_A!`](@ref juliajim.MDOFUTILS.HBRESFUN_A!), a routine that returns the harmonic balance residue along with an amplitude constraint. The usage of this is only slightly different from [`HBRESFUN!`](@ref juliajim.MDOFUTILS.HBRESFUN_A!) (which we have already encountered in Examples [C1](@ref ex_c1) and [C2](@ref ex_c2) for forced response analysis).

# The same system as before will be studied here:
# ```math
# \underbrace{\begin{bmatrix} 1&0\\0&1 \end{bmatrix}}_{\mx{M}} \begin{bmatrix} \ddot{x_1}\\ \ddot{x_2} \end{bmatrix} +
# \left( 0.01 \mx{M} + 0.001 \mx{K}\right) \begin{bmatrix} \dot{x_1}\\ \dot{x_2} \end{bmatrix} + \underbrace{\begin{bmatrix} 2&-1\\-1&2 \end{bmatrix}}_{\mx{K}} \begin{bmatrix} x_1\\x_2 \end{bmatrix} + \begin{bmatrix} 0\\ f_{fr}(x_2) \end{bmatrix} = \begin{bmatrix} 1\\0 \end{bmatrix}\cos\Omega t.
# ```

# ## Preamble: Load Packages
using GLMakie
using LinearAlgebra
using SparseArrays
using ForwardDiff
using NonlinearSolve
using DSP

using Revise #src
using juliajim.HARMONIC
using juliajim.CONTINUATION
using juliajim.MDOFUTILS

# ## [System Setup](@id exc3_setup)
# MCK matrices and nonlinearity setup: identical to [Example C2](@ref exc2_setup).

M = [1. 0.;0. 1.];
K = [2. -1.;-1. 2.];
C = 0.01*M+0.001*K;

mdl = MDOFGEN(M, C, K);

## Nonlinearity
kt = 1.0;
fs = 1.0;
fnl = (t,u,up,fp)-> if all(abs.(fp+kt*(u-up)).<fs)
    return fp+kt*(u-up), kt*ones(size(u)), -kt*ones(size(u)), ones(size(u)); else
        return fs*sign.(fp+kt*(u-up)), zeros(size(u)), zeros(size(u)), zeros(size(u));
end
L = [0. 1.];

mdl = ADDNL(mdl, :Hyst, fnl, L);

# ## [Setup HB](@id exc3_hbsetup)

# HB setup: identical to [Example C2](@ref exc2_hbsetup).
h = 1:2:5;
## h = 0:5;
N = 256;
t = (0:N-1)*2π/N;

Nhc = sum(all(h.==0, dims=2) + 2any(h.!=0, dims=2));

_, _, zinds, rinds, iinds = HINDS(mdl.Ndofs, h)
Fl = zeros(Nhc*mdl.Ndofs);
Fl[rinds[1]] = 1.0;

Wst = 0.6;
E, dEdw = HARMONICSTIFFNESS(mdl.M, mdl.C, mdl.K, Wst, h);
U0 = E\Fl;

# ### [Test HB Residue](@id exc3_testhb)

# The residue function [`HBRESFUN_A!`](@ref juliajim.MDOFUTILS.HBRESFUN_A!) expects a vector of the form \(\begin{bmatrix} \vc{u}& F& \Omega \end{bmatrix}^T\), where \(\vc{u}\) is the vector of harmonics, \(F\) is the (scalar) excitation amplitude (also an unknown), and \(\Omega\) is the (scalar) excitation frequency. Apart from the same arguments as [`HBRESFUN!`](@ref juliajim.MDOFUTILS.HBRESFUN_A!), the required response amplitude `Amp` must also be provided here. By default, this amplitude is interpreted as the first harmonic amplitude of the first DoF. This can be customized by supplying optional keywords `atype` (one of `:H1` or `:RMS`, specifying type of amplitude measure) and `ashape` (either DoF index or vector of weights to DoFs); see the [documentation](@ref juliajim.MDOFUTILS.HBRESFUN_A!).
# Here we set this up using `NonlinearSolve.jl` and solve for a single point. 

R = zeros(mdl.Ndofs*Nhc+1);
dRdUf = zeros(mdl.Ndofs*Nhc+1, mdl.Ndofs*Nhc+1);
dRdw = zeros(mdl.Ndofs*Nhc+1);
Amp = 1e0;
dof = 1;

HBRESFUN_A!([U0; 1.0; Wst], mdl, Amp, Fl, h, N; R=R, dRdUf=dRdUf, dRdw=dRdw)

fun = NonlinearFunction((r,uf,p)->HBRESFUN_A!([uf;p], mdl, Amp, Fl, h, N; R=r),
    jac=(J,uf,p)->HBRESFUN_A!([uf;p], mdl, Amp, Fl, h, N; dRdUf=J),
    paramjac=(Jp,uf,p)->HBRESFUN_A!([uf;p], mdl, Amp, Fl, h, N; dRdw=Jp));

prob = NonlinearProblem(fun, [U0;1.0], Wst);
sol = solve(prob, show_trace=Val(true));

# ## [Continuation](@id exc3_cont)
# Now we continue as before with continuation. In fact this piece of code is exactly the same as in [Example C2](@ref ex_c2).

Om0 = 0.1;
Om1 = 3;
dOm = 0.2;

cpars = (parm=:arclength, nmax=400, Dsc=:none, itopt=4);  # Autoscaling is not working
sols, its, dss, xis, Dsc = CONTINUATE([U0;1.0], fun, [Om0, Om1], dOm; cpars...);

uh = zeros(Complex, maximum(h)+1, length(sols), 2);
for i in 1:2
    uh[h.+1, :, i] = hcat([[s.up[zinds[i:2:end]];
                            s.up[rinds[i:2:end]]+1im*s.up[iinds[i:2:end]]]
                          for s in sols]...);
end
Oms = sols.p;
Fs = [s.up[end-1] for s in sols];

# ## [Plotting](@id exc3_plot)

# We plot first just the responses, to demonstrate how the amplitude constraint has been effective. 

# ### Amplitude Responses

his = [1, 3, 5];

set_theme!(theme_latexfonts())
fsz = 24;
fig1 = Figure(fontsize=fsz, size=(1000, 600));
if !isdefined(Main, :scr1) && isdefined(Main, :GLMakie) #src
   scr1 = GLMakie.Screen(); #src
end #src

for i in eachindex(his[his.<=maximum(h)])
    ax = Axis(fig1[1, i],
        ylabel=L"$H_%$(his[i])$ Response (m)", yscale=log10);
    scatterlines!(ax, Oms, abs.(uh[his[i].+1, :, 1]), label="x1");
    scatterlines!(ax, Oms, abs.(uh[his[i].+1, :, 2]), label="x2");
    if i==1
        axislegend(ax)
    end

    ax = Axis(fig1[2, i], xlabel=L"Excitation Frequency $\Omega$",
        ylabel=L"$H_%$(his[i])$ Phase (rad)");
    scatterlines!(ax, Oms, unwrap(angle.(uh[his[i].+1, :, 1])), label="x1");
    scatterlines!(ax, Oms, unwrap(angle.(uh[his[i].+1, :, 2])), label="x2");
end

if isdefined(Main, :GLMakie) #src
   display(scr1, fig1); #src
else #src
    fig1
end #src

# It can be observed that the first harmonic of the first DoF is always fixed to the value of `Amp` ($1.0$ here), while the other DoF's amplitude changes, showing that we have a fixed amplitude forced response here. 
# ### Forced Responses

# In order to visualize the nonlinear forced response function, we plot the response over force next. Here the resonances are clearly visible, albeit very different from what was seen in [Example C2](@ref ex_c2).

his = [1, 3, 5];

set_theme!(theme_latexfonts())
fsz = 24;
fig2 = Figure(fontsize=fsz, size=(1000, 600));
if !isdefined(Main, :scr2) && isdefined(Main, :GLMakie) #src
   scr2 = GLMakie.Screen(); #src
end #src

for i in eachindex(his[his.<=maximum(h)])
    ax = Axis(fig2[1, i],
        ylabel=L"$H_%$(his[i])$ Response (m/N)", yscale=log10);
    scatterlines!(ax, Oms, abs.(uh[his[i].+1, :, 1]./Fs), label="x1");
    scatterlines!(ax, Oms, abs.(uh[his[i].+1, :, 2]./Fs), label="x2");

    ax = Axis(fig2[2, i], xlabel=L"Excitation Frequency $\Omega$",
        ylabel=L"$H_%$(his[i])$ Phase (rad)");
    scatterlines!(ax, Oms, unwrap(angle.(uh[his[i].+1, :, 1])), label="x1");
    scatterlines!(ax, Oms, unwrap(angle.(uh[his[i].+1, :, 2])), label="x2");
end

if isdefined(Main, :GLMakie) #src
   display(scr2, fig2); #src
else #src
    fig2
end #src
