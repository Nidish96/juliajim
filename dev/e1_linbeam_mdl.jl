using GLMakie
using LinearAlgebra
using ProgressMeter
using NonlinearSolve
using StatsBase

using Revise
using juliajim.CONTINUATION
using juliajim.MDOFUTILS
using juliajim.HARMONIC

# * System Setup
function SSETUP(Ne)
    Ey = 70e9;
    rho = 2700.;
    w, h = 3*25e-3, 12.5e-3;
    ell = 1.0;
    ρA = rho*w*h;
    EI = Ey*w*h^3/12;

    Me(Le) = ρA*Le/420*[156 22Le 54 -13Le;
	                22Le 4Le^2 13Le -3Le^2;
	                54 13Le 156 -22Le;
	                -13Le -3Le^2 -22Le 4Le^2];
    Ke(Le) = EI/Le^3*[12 6Le -12 6Le;
		      6Le 4Le^2 -6Le 2Le^2;
		      -12 -6Le 12 -6Le;
		      6Le 2Le^2 -6Le 4Le^2];

    Nn = Ne+1;  # Number of nodes
    Xn = range(0, ell, Nn);

    M = zeros(2Nn, 2Nn);
    K = zeros(2Nn, 2Nn);
    for ei in 1:Ne
        is = 2(ei-1)+1;
        ie = 2(ei+1);

        M[is:ie, is:ie] += Me(diff(Xn[ei:ei+1])[1]);
        K[is:ie, is:ie] += Ke(diff(Xn[ei:ei+1])[1]);
    end
    F = [zeros(2Nn-2); 1.;0.];

    # Apply Clamped Boundary Condition
    Lb = I(2Nn)[:, 3:end];
    Mb = Lb'M*Lb;
    Kb = Lb'K*Lb;
    Fb = Lb'F;

    # Damping
    W0, V0 = eigen(Kb, Mb);
    W0 = sqrt.(W0);
    Zts = [0.2e-2, 0.1e-2];
    ab = [1 ./2W0[1:2] W0[1:2]/2]\Zts;
    Cb = ab[1]*Mb + ab[2]*Kb;

    return Xn, Lb, Mb, Cb, Kb, Fb
end

# * Build System
Ne = 10;
Xn, Lb, Mb, Cb, Kb, Fb = SSETUP(Ne);
Hfun(Om) = (Kb+1im*Om*Cb-Om^2*Mb)\Lb[end-1,:];
dHfun(Om) = -((Kb+1im*Om*Cb-Om^2*Mb)\ (1im*Cb-2Om*Mb))*Hfun(Om);

Om0 = 40.;
Om1 = 90.;
Nom = 1000;
Oms = range(Om0, Om1, Nom);

Uh = hcat(Hfun.(Oms)...);
dUh = hcat(dHfun.(Oms)...);

U = [real(Uh);-imag(Uh)];
dU = [real(dUh);-imag(dUh)];

tgt = hcat(normalize.(eachcol([dU;ones(1,Nom)]))...);

mdl = MDOFGEN(Mb,Cb,Kb);

# * HB for response
h = (0:5);
N = 128;
Nhc = NHC(h);
_,_, zinds,rinds,iinds = HINDS(mdl.Ndofs, h);

Fl = zeros(mdl.Ndofs*Nhc);
Fl[rinds[1:mdl.Ndofs]] = Fb;

E0, _ = HARMONICSTIFFNESS(mdl.M, mdl.C, mdl.K, Om0, h);
U0_(Fa) = E0\ (Fl*Fa);

fun(Fa) = NonlinearFunction((r,u,p)->HBRESFUN!([u;p], mdl, Fa*Fl, h, N; R=r),
    jac=(J,u,p)->HBRESFUN!([u;p], mdl, Fa*Fl, h, N; dRdU=J),
    paramjac=(Jp,u,p)->HBRESFUN!([u;p], mdl, Fa*Fl, h, N; dRdw=Jp));

Famp = 1.0;

# Continuation
cpars = (nmax=200, Dsc=:auto, angopt=deg2rad(20));
dOm = 2.5;
sols, its, dss, xis, Dsc = CONTINUATE(U0_(Famp), fun(Famp), [Om0, Om1], dOm; cpars...);

set_theme!(theme_black())
fsz = 18;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr) && Makie.current_backend()==GLMakie
   scr = GLMakie.Screen();
end

ax = Axis(fig[1, 1], xlabel="Excitation Frequency (rad/s)", ylabel="Response (m)");
lines!(ax, Oms, abs.(kron([1 1im], Lb[end,:]')*U)[:])
scatterlines!(ax, sols.p,
    abs.(kron([1 1im], Lb[end,:]')*
         hcat(sols.u...)[[rinds[1:mdl.Ndofs]; iinds[1:mdl.Ndofs]],:])[:], color=:red)

if Makie.current_backend()==GLMakie
   display(scr, fig);
else
   display(fig)
end
   
