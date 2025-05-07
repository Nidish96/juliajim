using LinearAlgebra
using NLsolve
using UnPack
using Plots
include("../ROUTINES/HARMONIC.jl")
include("../ROUTINES/CONTINUATION.jl");

function EPMCRESFUN!(R, dRdUxw, dRda, Uxwa, pars, h,Nt,L)
    @unpack w0, mu, alpha = pars
    
    Nhc = sum(all(h.==0, dims=2) + 2*any(h.!=0, dims=2));
    a = Uxwa[end];
    A = 10^a;
    dAda = A*log(10);

    Asc = ones(Nhc,1);
    Asm = zeros(Nhc,1);
    if h[1]==0
        Asc[1]=0;
        Asm[1]=1;

        M = L[2:3,:]'*L[2:3,:];
    else
        M = L[1:2,:]'*L[1:2,:];
    end
    Usc = (L*Uxwa[1:end-3]).*(Asm+Asc*A);

    w = Uxwa[end-1];
    xi = Uxwa[end-2];

    # Linear Portion
    E, dEdw = HARMONICSTIFFNESS(1, 2*mu-xi, w0^2, [w], h);
    dEdxi, _ = HARMONICSTIFFNESS(0, -1, 0, [w], h);
    dEdw = dEdw[1];
    
    # AFT for Nonlniear Force
    ut = AFT(Usc, h, Nt, :f2t);
    ft = alpha*ut.^3;
    Fh = AFT(ft, h, Nt, :t2f);

    # Construct Residue
    if !(R === nothing)
        R[:] = [E*Usc + Fh;
                Uxwa[1:end-3]'*M*Uxwa[1:end-3]-1];
    end
    if !(dRdUxw === nothing)
        cst = AFT(I(Nhc), h, Nt, :f2t);
        dfdu = (3*alpha*ut.^2) .* cst;
        Jh = AFT(dfdu, h, Nt, :t2f);

        dRdUxw[:, :] = vcat(hcat((E+Jh)*((Asm+Asc*A).*L), dEdxi*Usc, dEdw*Usc),
                            hcat(2*Uxwa[1:end-3]'*M, 0, 0));
    end
    if !(dRda === nothing)
        dRda[:] = [(E+Jh)*((Asc*dAda).*L)*Uxwa[1:end-3];0];
    end
end

###
savfig = false;
pars = (w0=2, mu=0.01, alpha=0.1);

h = HSEL(9, 1);
Nt = 512;
Nhc = sum(all(h.==0, dims=2) + 2*any(h.!=0, dims=2));

L = I(Nhc);
L = L[:, 1:end .!= 3];

###
U0 = zeros(Nhc-1,1);
U0[2] = 1;

Uxw0 = vcat(U0, 2*pars.mu, pars.w0);

Ast = 0;
Aen = 1;
ΔA = 0.01;

UxwC = SOLVECONT((F,J,Jw,X)->EPMCRESFUN!(F, J, Jw, X, pars, h,Nt,L),
                 Uxw0, Ast, Aen, ΔA);

###
As = 10 .^ UxwC[end,:];
wres = UxwC[end-1,:];
X3 = (L[6,:] .- 1im .* L[7,:])'*UxwC[1:end-3,:];
X5 = (L[10,:] .- 1im .* L[11,:])'*UxwC[1:end-3,:];

wres_O1 = pars.w0 .+ 3*pars.alpha/(8*pars.w0) .* (As .^ 2);
wres_O2 = wres_O1 .- (pars.mu^2/(2*pars.w0) .+ 15*pars.alpha^2/(256*pars.w0^3) .* (As .^ 4) );

X3_O1 = pars.alpha/(32*pars.w0^2) .* (As .^ 3);
X3_O2 = X3_O1 .- (3im*pars.alpha*pars.mu/(64*pars.w0^3) .* (As .^ 3) .+ 21*pars.alpha^2/(1024*pars.w0^4) .* (As .^ 5));

X5_O1 = zeros(size(X3_O1));
X5_O2 = pars.alpha^2/(1024*pars.w0^4) .* (As .^ 5);

###
p1 = plot(As, [wres_O1 wres_O2 wres], label=["O(ϵ)" "O(ϵ²)" "EPMC"],
          lw=3, ls=[:solid :solid :dashdot], lc=["red" "blue" "black"],
          grid=true, xscale=:log10, 
          ylabel="ω̃")
p2 = plot(As, abs.([X3_O1 X3_O2 X3[:]]), label=false,
          lw=3, ls=[:solid :solid :dashdot], lc=["red" "blue" "black"],
          grid=true, xscale=:log10, yscale=:log10, 
          ylabel="H3")
p3 = plot(As, abs.([X5_O2 X5[:]]), legend=false,
          lw=3, ls=[:solid :dashdot], lc=["blue" "black"],
          grid=true, xscale=:log10, yscale=:log10,
          xlabel="Modal Amplitude", ylabel="H5")
plot(p1, p2, p3, layout=(3,1), size=(600,600))
if savfig
    png("./FIGS/D_duffcomp.png")
end
