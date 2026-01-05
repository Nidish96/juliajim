```@meta
EditURL = "../../examples/f_febeam.jl"
```

```@meta
CurrentModule = juliajim
```

## [Example F: Finite Element Beam with Nonlinear Support](@id ex_f)

This example considers a finite element model of a Beam (Euler-Bernoulli) that has a nonlinear attachment in one end and is clamped in the other.
The same code can be used for simulating 3 different types of nonlinear attachments:
1. A simple cubic spring of the form $\beta u^3$.
2. A stiff string attachment with reaction force of the form $T \frac{w_T}{\sqrt{\ell_s^2+w_T^2}}$
3. A non-smooth frictional joint.

## Preamble: Load Packages

````@example f_febeam
using GLMakie
using LinearAlgebra
using Random
using NonlinearSolve
using SparseArrays
using Revise

using juliajim.HARMONIC
using juliajim.MDOFUTILS
using juliajim.CONTINUATION
````

## Setup System

````@example f_febeam
Ey = 70e9;
rho = 2700.;
w, h = 3*25e-3, 12.5e-3;
ell = 1.0;
ρA = rho*w*h;
EI = Ey*w*h^3/12;
nothing #hide
````

Cubic Spring

````@example f_febeam
β = 500.0;  # Nonlinearity
fnl = (t, u, ud) -> return β.*u.^3, 3β.*u.^2, zeros(size(u));
nothing #hide
````

Stiffened String

````@example f_febeam
Ts = 1e2;
ls = 0.10;
fnl = (t, u, ud) -> return Ts.*u./sqrt.(ls^2 .+u.^2),
    Ts*ls^2 ./sqrt.(ls^2 .+u.^2).^3, zeros(size(u));
typ = :Inst;
nothing #hide
````

Frictional Support

````@example f_febeam
kt = 500;
fs = 1e-2kt;
fnl = (t,u,up,fp)-> if all(abs.(fp+kt*(u-up)).<fs)
    return fp+kt*(u-up), kt*ones(size(u)), -kt*ones(size(u)), ones(size(u)); else
        return fs*sign.(fp+kt*(u-up)), zeros(size(u)), zeros(size(u)), zeros(size(u));
end
typ = :Hyst;

Ne = 10;  # Number of elements
Nn = Ne+1;  # Number of nodes
Xn = range(0, ell, Nn);

Me(Le) = ρA*Le/420*[156 22Le 54 -13Le;
	            22Le 4Le^2 13Le -3Le^2;
	            54 13Le 156 -22Le;
	            -13Le -3Le^2 -22Le 4Le^2];
Ke(Le) = EI/Le^3*[12 6Le -12 6Le;
		  6Le 4Le^2 -6Le 2Le^2;
		  -12 -6Le 12 -6Le;
		  6Le 2Le^2 -6Le 4Le^2];

M = zeros(2Nn, 2Nn);
K = zeros(2Nn, 2Nn);
for ei in 1:Ne
    is = 2(ei-1)+1;
    ie = 2(ei+1);

    M[is:ie, is:ie] += Me(diff(Xn[ei:ei+1])[1]);
    K[is:ie, is:ie] += Ke(diff(Xn[ei:ei+1])[1]);
end
F = [zeros(2Nn-2); 1.;0.];
nothing #hide
````

Apply Clamped Boundary Condition

````@example f_febeam
Lb = I(2Nn)[:, 3:end];
Mb = Lb'M*Lb;
Kb = Lb'K*Lb;
Fb = Lb'F;

Wn = sqrt.(eigvals(Kb, Mb));
Zts = [0.2e-2, 0.1e-2];
ab = [1 ./2Wn[1:2] Wn[1:2]/2]\Zts;
Cb = ab[1]*Mb + ab[2]*Kb;
nothing #hide
````

Setup Nonlinearity

````@example f_febeam
mdl = MDOFGEN(Mb, Cb, Kb);
mdl = ADDNL(mdl, typ, fnl, Float64.(Lb[end-1:end-1,:]));
nothing #hide
````

## Setup HB

````@example f_febeam
h = 0:5;
N = 256;

Nhc = NHC(h);
_, _, zinds, rinds, iinds = HINDS(mdl.Ndofs, h);

Fl = zeros(Nhc*mdl.Ndofs);
Fl[rinds[1:mdl.Ndofs]] = Fb;

Om0 = 10.0;
Om1 = 450.0;
dOm = 10.0;
nothing #hide
````

Linear Initial Guess

````@example f_febeam
E, _ = HARMONICSTIFFNESS(mdl.M, mdl.C, mdl.K, Om0, h);
U0_(Fa) = E\ (Fa*Fl);
nothing #hide
````

Setup nonlinear function

````@example f_febeam
fun(Fa) = NonlinearFunction((r,u,p)->HBRESFUN!([u;p], mdl, Fa*Fl, h, N; R=r),
    jac=(J,u,p)->HBRESFUN!([u;p], mdl, Fa*Fl, h, N; dRdU=J),
    paramjac=(Jp,u,p)->HBRESFUN!([u;p], mdl, Fa*Fl, h, N; dRdw=Jp));
nothing #hide
````

## Forced Response Continuation

````@example f_febeam
Famp = 8.0  # 0.1, 1.0, 2.0, 4.0, 8.0, 16.0

U0 = U0_(Famp);
prob = NonlinearProblem(fun(Famp), U0, Om0; abstol=1e-6, reltol=1e-6);
sol = solve(prob, show_trace=Val(true));
nothing #hide
````

Continuation

````@example f_febeam
cpars = (parm=:arclength, nmax=2000, Dsc=:auto);

Om0 = 40.0;
Om1 = 100.0;
dOm = 2.5;
sols, _, _, _, _ = CONTINUATE(U0, fun(Famp), [Om0, Om1], dOm; cpars...);

Uh = [Lb[end-1,:]'reshape(u, mdl.Ndofs,:) for u in sols.u];
nothing #hide
````

### Stability Certification

````@example f_febeam
E0, _ = HARMONICSTIFFNESS(zeros(mdl.Ndofs,mdl.Ndofs), -2mdl.M,
    zeros(mdl.Ndofs, mdl.Ndofs), [1.], h);
indsh1 = [rinds[1:mdl.Ndofs]; iinds[1:mdl.Ndofs]];
E0 = collect(E0);
J = zeros(mdl.Ndofs*Nhc, mdl.Ndofs*Nhc);
stab = zeros(length(sols));
for iw in 1:length(sols)
    fun(Famp).jac(J, sols.u[iw], sols.p[iw])
    eVs = eigvals(J[indsh1,indsh1], sols.p[iw]*E0[indsh1,indsh1]);
    stab[iw] = sum(real(eVs).>=0);
end
````

## Plot Forced Response

````@example f_febeam
set_theme!(theme_latexfonts())
fsz = 18;
fig = Figure(fontsize=fsz);

ax = Axis(fig[1, 1], xlabel="Excitation Frequency (rad/s)",
    ylabel="Response (m)", yscale=log10);
scatterlines!(ax, sols.p./(stab.==0), norm.(Uh)/Famp)
scatterlines!(ax, sols.p./(stab.!=0), norm.(Uh)/Famp)

xlims!(ax, Om0, Om1)
    fig
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

