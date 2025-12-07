```@meta
EditURL = "../../examples/d2_forcedvdp.jl"
```

```@meta
CurrentModule = juliajim
```

## * [Example D2: Bifurcation Analysis of the Forced Response of a Van der Pol Oscillator](@id ex_d2)

This example is designed to demonstrate stability analysis and branch switching across a (secondary) Hopf bifurcation that occurs in a forced Van der Pol oscillator.

The dynamical system that will be studied here is:
```math
\ddot{x} - c \dot{x} + k x + \eta x^2 \dot{x} = F\cos\Omega t.
```
This is an SDoF nonlinear oscillator with negative viscous damping. The nonlinearity provides the dissipation that saturates the response. Classically known as the Van der Pol oscillator, this is a useful system for learning about limit cycles. In the considered version, we also have an external periodic excitation.

The system will have a stable periodic response close to resonance, which will lose its stability as we move away from it.

## Preamble: Load Packages

````@example d2_forcedvdp
using GLMakie
using LinearAlgebra
using SparseArrays
using NonlinearSolve
using DSP
using Debugger

using Revise
using juliajim.HARMONIC
using juliajim.CONTINUATION
using juliajim.MDOFUTILS
using juliajim.NLDYN
````

## [System Setup](@id exd2_setup)
Here we setup the system with parameters.

````@example d2_forcedvdp
pars = (c=0.02, k=4., F=1., eta=0.1);

mdl = MDOFGEN(1.0, -pars.c, pars.k);
fnl = (t,u,ud) -> return pars.eta.*u.^2 .*ud, 2pars.eta.*u.*ud, pars.eta.*u.^2;
mdl = ADDNL(mdl, :Inst, fnl, 1.0);
nothing #hide
````

## [Setup HB](@id exd2_hbsetup)

````@example d2_forcedvdp
h = 0:5;
N = 128;
t = range(0, 2π, N+1)[1:N];

Nhc = NHC(h);

inds0, hinds, zinds, rinds, iinds = HINDS(mdl.Ndofs, h);
Fl = zeros(Nhc*mdl.Ndofs, 1);
Fl[rinds[1]] = pars.F;

Om0 = 0.1;
Om1 = 2.0;
E, dEdw = HARMONICSTIFFNESS(mdl.M, mdl.C, mdl.K, Om0, h);

Uw0 = [E\Fl; Om0];
nothing #hide
````

## Setup HB Residue, Get First Point

````@example d2_forcedvdp
fun = NonlinearFunction((r,u,p)->HBRESFUN!([u;p], mdl, Fl, h, N; R=r),
    jac=(J,u,p)->HBRESFUN!([u;p], mdl, Fl, h, N; dRdU=J),
    paramjac=(Jp,u,p)->HBRESFUN!([u;p], mdl, Fl, h, N; dRdw=Jp));

prob = NonlinearProblem(fun, Uw0[1:end-1], Uw0[end]);
sol = solve(prob, show_trace=Val(true));
nothing #hide
````

## [Continuation](@id exd2_cont)
Just like in other examples we call the [`CONTINUATE`](@ref juliajim.CONTINUATION.CONTINUATE) utility to obtain the periodic forced response.

````@example d2_forcedvdp
Om0 = 0.1;
Om1 = 4.0;
dOm = 0.2;
cpars = (parm=:arclength, nmax=100, save_jacs=true);

sols, _, _, _, _ = CONTINUATE(Uw0[1:end-1], fun, [Om0, Om1], dOm; cpars...);
nothing #hide
````

Obtain Harmonics

````@example d2_forcedvdp
uh = zeros(Complex, maximum(h)+1, length(sols));
uh[[inds0; hinds], :] = hcat([[up[zinds];up[rinds,:]+1im*up[iinds,:]] for up in sols.up]...);
Oms = [up[end] for up in sols.up];
nothing #hide
````

### Stability Analysis
We now use an averaging formulation to obtain the stability coefficients based on just the first harmonics. This will work if the response is dominantly single harmonic.

````@example d2_forcedvdp
E0, _ = HARMONICSTIFFNESS(0., -2mdl.M, 0, 1, h);
E0 = collect(E0);
stab = zeros(length(Oms));
for (i,(J,Om)) in enumerate(zip(sols.J, Oms))
````

evs = eigvals(J[2:end,2:end], Om*collect(E0[2:end,2:end]));  # Multiharmonic

````@example d2_forcedvdp
    evs = eigvals(J[2:3,2:3], Om*E0[2:3,2:3]);
    stab[i] = sum(real(evs).>=0);
end
````

The variable `stab` will store the number of unstable eigenvalues that have been detected. One may interpret these as the Floquet exponents of the system.

## Bifurcation Analysis

### Setup QPHB for Branch Switching

````@example d2_forcedvdp
Nhmax = 3;
hq = HSEL(Nhmax, [1.,1.]);

inds0q,hindsq, zindsq,rindsq,iindsq = HINDS(1, hq);
Nhcq = NHC(hq);
Nq = 32;

h01 = findall(vec(all(hq.==0, dims=2)));
hq1s = setdiff(findall(hq[:,2].==0), h01).-length(h01);  # w1 alone
hq2s = setdiff(findall(hq[:,1].==0), h01).-length(h01);  # w2 alone

Flq = zeros(Nhcq);
Flq[rindsq[hq2s[1]]] = pars.F;
nothing #hide
````

### Detection

````@example d2_forcedvdp
bifis = findall(stab[1:end-1].!=stab[2:end]);
bifis[2] += 1;
dxis = [-1, 1];

uhbs = [];
Omsb = [];
Wselfs = [];
nothing #hide
````

We loop over the two detected bifurcation points

````@example d2_forcedvdp
for (bi, bifi) in enumerate(bifis)
````

# Eigenanalysis to obtain unstable manifold (eigenvectors)

````@example d2_forcedvdp
    eVals, eVecs = eigen(sols.J[bifi][2:3, 2:3], Oms[bifi]*collect(E0[2:3,2:3]));
    eVecsC = eVecs[1,:]-1im.*eVecs[2,:];  # Complexify

    ei = argmax(abs.(eVecsC));  # Only one should "survive" for the Hopf bifurcation

    sig = imag(eVals[ei]);
    Wself = Oms[bifi] + sig;  # Self Excited Frequency
    eVecsC = eVecsC*exp(-1im*angle(eVecsC[ei]));  # Normalize phase
    Pvec = normalize([real(eVecsC[ei]), -imag(eVecsC[ei])]);  # Perturbation Vector
nothing #hide
````

Setup Initial Guess

````@example d2_forcedvdp
    Uhq0 = zeros(Nhcq);
    Uhq0[[rindsq[hq2s]; iindsq[hq2s]]] = sols.up[bifi][[rinds[1:length(hq2s)];
                                                        iinds[1:length(hq2s)]]];
nothing #hide
````

Perturbation Vector

````@example d2_forcedvdp
    Phq0 = zeros(Nhcq);
    Phq0[[rindsq[hq1s[1]]; iindsq[hq1s[1]]]] = Pvec;
nothing #hide
````

Perturbation amplitude

````@example d2_forcedvdp
    qamp = 1.0;
nothing #hide
````

We choose this arbitrarily for now. It is possible to use the
method of normal forms to fix this exactly, see [next example](@ref ex_d3).

### Apply phase constraint and converge with deflation

````@example d2_forcedvdp
    cL = I(Nhcq+2)[:, setdiff(1:Nhcq+2, iindsq[hq1s[1]])];
    uC = cL'*[Uhq0; Wself; Oms[bifi]];
nothing #hide
````

The matrix `cL` is defined so that the original solution vector `Uw`
(which includes harmonics and both the frequency components) can be
recovered by \(Uw = cL \hat{Uw}\).

````@example d2_forcedvdp
    funq = NonlinearFunction((r,u,p)-> QPHBRESFUN!([u;p], mdl, Flq, hq, Nq;
        R=r,cL=cL,U0=uC),
        jac=(J,u,p)->QPHBRESFUN!([u;p], mdl, Flq, hq, Nq;
            dRdU=J,cL=cL,U0=uC),
        paramjac=(Jp,u,p)->QPHBRESFUN!([u;p], mdl, Flq, hq, Nq;
            dRdw=Jp,cL=cL,U0=uC));
nothing #hide
````

The [`QPHBRESFUN!`](@ref juliajim.MDOFUTILS.QPHBRESFUN!) function
supports deflation specification through the keyword argument U0.

````@example d2_forcedvdp
    u0 = cL[1:end-1,:]'*[Uhq0+qamp*Phq0; Wself];
    probq = NonlinearProblem(funq, u0[1:end-1], Oms[bifi]);
    solq = solve(probq, show_trace=Val(true));
nothing #hide
````

Get one point before (without deflation)

````@example d2_forcedvdp
    funr = NonlinearFunction((r,u,p)-> QPHBRESFUN!([u;p], mdl, Flq, hq, Nq; R=r,cL=cL),
        jac=(J,u,p)->QPHBRESFUN!([u;p], mdl, Flq, hq, Nq; dRdU=J,cL=cL));
    probq = NonlinearProblem(funr, u0[1:end-1], Oms[bifi-dxis[bi]]);
    solprev = solve(probq, show_trace=Val(true));
nothing #hide
````

### Continue away from the bifurcation point

````@example d2_forcedvdp
    Om0b = Oms[bifi];
    Om1b = (dxis[bi]<0) ? Om0 : Om1;
    dOmb = 0.2;
    cparsb = (parm=:arclength, nmax=100, save_jacs=true);

    solsb, _, _, _, _ = CONTINUATE(solq.u, funq, [Om0b, Om1b], dOmb; cparsb...);
nothing #hide
````

Prepend previous point & expand constraint

````@example d2_forcedvdp
    sup = [cL*up for up in [[solprev.u;Oms[bifi-dxis[bi]]], solsb.up...]];
nothing #hide
````

Obtain Harmonics

````@example d2_forcedvdp
    uhq = zeros(Complex, size(hq,1), length(sup));
    uhq[[inds0q; hindsq], :] = hcat([[up[zindsq];up[rindsq,:]+1im*up[iindsq,:]]
                                    for up in sup]...);

    push!(uhbs, uhq);
    push!(Omsb, [up[end] for up in sup]);
    push!(Wselfs, [up[end-1] for up in sup]);
end
````

##  [Plotting](@id exd2_plot)

````@example d2_forcedvdp
set_theme!(theme_latexfonts())
fsz = 24;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr) && Makie.current_backend()==GLMakie
   scr = GLMakie.Screen();
end

ax = Axis(fig[1, 1], xlabel=L"Excitation Frequency $\Omega$", ylabel="Response");
scatterlines!(ax, Oms./(stab.==0), [norm(u) for u in eachcol(uh)], label="Stable")
scatterlines!(ax, Oms./(stab.==2), [norm(u) for u in eachcol(uh)], label="Unstable")

for (i, (om,uhb)) in enumerate(zip(Omsb, uhbs))
    scatterlines!(ax, om, [norm(u) for u in eachcol(uhb)],
        label="Branch $(i)")
end

axislegend(ax, nbanks=3, position=:ct)
xlims!(Om0, Om1)
ylims!(0, 3.75)
if Makie.current_backend()==GLMakie
   display(scr, fig);
else
   display(fig)
end
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

