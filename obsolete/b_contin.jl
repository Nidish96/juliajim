# * Preamble
# using PyPlot
using GLMakie
using LinearAlgebra
using NLsolve
using UnPack
using DSP

using Revise
includet("../src/HARMONIC.jl")
includet("../src/CONTINUATION.jl")

# * Residue Function
function RESFUN!(R, dRdU, dRdws, Uws, Fl, pars, h, Nt)
    @unpack m, c, k, bt = pars
    C = size(h,2);
    ws = Uws[end-C+1:end];
    Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));

    # Linear Portion
    E, dEdws = HARMONICSTIFFNESS(m, c, k, ws, h);

    # AFT For nonlinear force
    ut  = AFT(Uws[1:end-C], h, Nt, :f2t);
    ft = bt*ut.^3;
    Fnl    = AFT(ft, h, Nt, :t2f);

    # Construct Residue
    if !(R === nothing)        
        R[:] = E*Uws[1:end-C] + Fnl - Fl;
    end
    if !(dRdU === nothing)
        cst = AFT(I(Nhc), h, Nt, :f2t);    
        dfdat = (3*bt*ut.^2) .* cst;
        Jnl    = AFT(dfdat, h, Nt, :t2f);
        
        dRdU[:, :] = E + Jnl;
    end
    if !(dRdws === nothing)
        dRdws[:, :] = hvncat(2, map(dE->dE*Uws[1:end-C], dEdws)...);
    end
    return nothing;
end

# * Investigation
# pars = (m = 1.0, c = 0.01, k = 4.0, bt = 0.01);
pars = (m = 1.0, c = 0.015, k = 4.0, bt = 0.01);
famp = 1.25; # 1.25;

C = 1;
Nhmax = 3;
Nt = 128;

# h = HSEL(Nhmax, 1:C);
h = HSEL(Nhmax, 1.);
Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));

Fl = zeros(Nhc, 1);
Fl[2] = 1.0;

Wst = 0.1;
Wen = 4.0;

E,dEdws = HARMONICSTIFFNESS(pars.m, pars.c, pars.k, [Wst], h);
U0 = E\Fl * famp;

ds = 0.3;
Dscale = max.(ones(Nhc+1)*1e-4, abs.([U0; Wst]))[:,1];
# Dscale = ones(Nhc+1)*1e-1;
xΩ = SOLVECONT((F,J,Jw,Uw)->RESFUN!(F, J, Jw, Uw, Fl * famp, pars, h, Nt),
               U0, Wst, Wen, ds, Dscale=Dscale, DynDscale=true, itopt=3,
               stepmax=500);

# ** Stability
dRdU = zeros(Nhc, Nhc);
λ1s = zeros(ComplexF64, Nhc, size(xΩ,2));
λ2s = zeros(ComplexF64, 2*Nhc, size(xΩ,2));
for i=eachindex(xΩ[1,:])
    RESFUN!(nothing, dRdU, nothing, xΩ[:, i], Fl * famp, pars, h, Nt);
    E1,_ = HARMONICSTIFFNESS(0, 2*pars.m, pars.c, [xΩ[end, i]], h);
    E2,_ = HARMONICSTIFFNESS(0, 0, pars.m, [xΩ[end, i]], h);

    peph = PEP([dRdU, Matrix(E1), Matrix(E2)]);
    pevs,_ = polyeig(peph);
    is = sortperm(abs.(imag(pevs)));
    λ2s[:, i] = pevs[is];

    λ1s[:, i] = eigvals(-dRdU, Matrix(E1))
end
stab = real(λ2s[1,:]).<0;
istab = conv(.~stab, [1;1;1]);
istab = istab[2:end-1].!=0;

# * Plotting
mrk = :none;
lw = 3;

fsz = 18;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr)
   scr = GLMakie.Screen();
end

ax = Axis(fig[1, 1], xlabel="Excitation Frequency", ylabel="H1 Response (m)");
lines(xΩ[end,:]./stab, abs.(xΩ[2,:]-im*xΩ[3,:]));

display(scr, fig);


# p1 = plot(xΩ[end,:]./stab, abs.(xΩ[2,:]-im*xΩ[3,:]),
#           ls=:solid, lw=lw, markershape=mrk, ms=2,
#           label="stable", grid=true,
#           ylabel="H1 Response (m)")
# plot!(xΩ[end,:]./istab, abs.(xΩ[2,:]-im*xΩ[3,:]),
#       ls=:dash, lw=lw, markershape=mrk, ms=2,
#       label="unstable", grid=true)

# p2 = plot(xΩ[end,:]./stab, rad2deg.(angle.(xΩ[2,:]-im*xΩ[3,:])),
#           ls=:solid, lw=lw, markershape=mrk, ms=2,
#           label=false, grid=true, yticks=-180:45:0,
#           ylabel="H1 Phase (degs)")
# plot!(xΩ[end,:]./istab, rad2deg.(angle.(xΩ[2,:]-im*xΩ[3,:])),
#       ls=:dash, lw=lw, markershape=mrk, ms=2,
#       label=false, grid=true) 

# plot(p1, p2, layout=(2,1), size=(600,600))

# # TODO: Make the continuation better!!
