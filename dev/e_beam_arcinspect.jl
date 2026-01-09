using GLMakie
using LinearAlgebra
using ProgressMeter

using Revise
using juliajim.CONTINUATION

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

# * Compute in Steps
Emat(Om) = [Kb-Om^2*Mb Om*Cb;-Om*Cb Kb-Om^2*Mb];
dEmat(Om) = [-2Om*Mb Cb;-Cb -2Om*Mb];
Fv = kron([1;0], Fb);

atgt = 1e-1;
nxi = 0.95;

dOm0 = 0.5;
dOm = copy(dOm0);
Nsteps = 10000;
Us = zeros(4Ne+1, Nsteps);
dUs = zeros(4Ne+1, Nsteps);
Vs = zeros(4Ne+1, Nsteps);
Us[:, 1] = [Emat(Om0)\Fv; Om0];
dUs[:, 1] = normalize([-Emat(Om0)\dEmat(Om0)*Us[1:end-1,is]; 1]);

~, ~, V0 = svd([Emat(Om0) dEmat(Om0)*Us[1:end-1,is]]);
Vs[:, 1] = V0[:, end];

Dsc0 = norm.(eachcol([Emat(Om0) dEmat(Om0)\Fv]));
Dsc = copy(Dsc0);
ndxi = 0.5;

angs = [];
ns = Nsteps;
for is in 2:Nsteps
    Om = Us[end, is-1]+dOm;
    
    Us[:, is] = [Emat(Om)\Fv; Om];    
    dUs[:, is] = normalize([-Emat(Om)\dEmat(Om)*Us[1:end-1,is]; 1]);

    _, _, V = svd([Emat(Om) dEmat(Om)*Us[1:end-1,is]]);
    Vs[:, is] = V[:, end]*sign(Vs[:, is-1]'V[:, end]);

    Dsc_ = norm.(eachcol([Emat(Om) dEmat(Om)\Fv]));
    xiDsc = clamp.((Dsc_./Dsc).^ndxi, 0.5, 2.0);
    global Dsc = xiDsc.*Dsc;
    
    # ang = acos(normalize(Dsc.*dUs[:,is-1])'normalize(Dsc.*(Us[:,is]-Us[:,is-1])))/dOm;
    # ang = acos(Vs[:, is]'Vs[:, is-1])/dOm;

    ninds = findall(dUs[:, is-1].!=0);
    ang = acos(normalize(ones(length(ninds)))'normalize((Us[ninds,is]-Us[ninds,is-1])./dUs[ninds,is-1]))
    
    push!(angs, ang);    
    if is>2
        xi = clamp((atgt/ang)^(nxi), 0.5, 2.0);
        println("$is) Om=$Om, dOm=$dOm, ang=$ang, atgt=$atgt");
        global dOm = clamp(dOm*xi, dOm0/5, 5dOm0);
        
        if Us[end, is]>Om1
            global ns = is;
            break;
        end
    end

    # global atgt = sqrt(atgt*ang);
end
Us = Us[:, 1:ns];
dUs = dUs[:, 1:ns];
Vs = Vs[:, 1:ns];


# Plot
set_theme!(theme_latexfonts())
fsz = 18;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr) && Makie.current_backend()==GLMakie
   scr = GLMakie.Screen();
end

ax = Axis(fig[1, 1], xlabel="Excitation Frequency (rad/s)",
    ylabel="Response (m/N)", yscale=log10);
lines!(ax, Oms, abs.(Lb[end-1,:]'*Uh)[:]);
scatterlines!(ax, Us[end,:], abs.(kron([1,1im], Lb[end-1,:])'*Us[1:end-1,:])[:],
    color=:red);

ax = Axis(fig[2, 1], xlabel="Excitation Frequency (rad/s)",
    ylabel="Derivative (m/N/rad/s)", yscale=log10);
lines!(ax, Oms[2:end], ((x,y)->acos(x'y)).(eachcol(tgt[:,2:end]),
    eachcol(tgt[:,1:end-1])));
scatterlines!(ax, Us[end,2:end], angs);

if Makie.current_backend()==GLMakie
   display(scr, fig);
else
   display(fig)
end
   

