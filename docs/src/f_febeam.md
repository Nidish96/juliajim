```@meta
EditURL = "../../examples/f_febeam.jl"
```

```@meta
CurrentModule = juliajim
```

## [Example F: Finite Element Beam]

This example considers a finite element model of a Beam (Euler-Bernoulli) that has a nonlinear attachment in one end and is clamped in the other.

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

β = 500.0;  # Nonlinearity
fnl = (t, u, ud) -> return β.*u.^3, 3β.*u.^2, zeros(size(u));

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
mdl = ADDNL(mdl, :Inst, fnl, Float64.(Lb[end-1:end-1,:]));
nothing #hide
````

* Setup HB

````@example f_febeam
h = 0:1;
N = 128;

Nhc = NHC(h);
_, _, zinds, rinds, iinds = HINDS(mdl.Ndofs, h);

Fl = zeros(Nhc*mdl.Ndofs);
Fl[rinds[1:mdl.Ndofs]] = Fb;

Famp = 1.0;

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

* Setup nonlinear function

````@example f_febeam
fun(Fa) = NonlinearFunction((r,u,p)->HBRESFUN!([u;p], mdl, Fa*Fl, h, N; R=r),
    jac=(J,u,p)->HBRESFUN!([u;p], mdl, Fa*Fl, h, N; dRdU=J),
    paramjac=(Jp,u,p)->HBRESFUN!([u;p], mdl, Fa*Fl, h, N; dRdw=Jp));

Famp = 20.0  # 100

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
dOm = 5.0;
sols, _, _, _, _ = CONTINUATE(U0, fun(Famp), [Om0, Om1], dOm; cpars...);

Uh = [Lb[end-1,:]'reshape(u, mdl.Ndofs,:) for u in sols.u];
nothing #hide
````

*

````@example f_febeam
set_theme!(theme_latexfonts())
fsz = 18;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr) && Makie.current_backend()==GLMakie
   scr = GLMakie.Screen();
end

ax = Axis(fig[1, 1], xlabel="Excitation Frequency (rad/s)",
    ylabel="Response (m)");
scatterlines!(ax, sols.p, norm.(Uh)/Famp)

if Makie.current_backend()==GLMakie
   display(scr, fig);
else
   display(fig)
end
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

