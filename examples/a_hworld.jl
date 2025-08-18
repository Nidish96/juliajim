# * Preamble
using GLMakie
using LinearAlgebra
using SparseArrays
using Arpack

using Revise
using juliajim.HARMONIC

# * AFT - 1 Freq Check
N = 128;
t = (0:N-1)*2π/N;
# yt = cos.(t) + 3sin.(2t) .+ 4;
yt = cos.(t) .+ 2;
yt = [yt yt];

h = collect(0:3)[:,:];
Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));

yf = AFT(yt, h, N, :t2f);
YT = AFT(yf, h, N, :f2t);

# ** AFT - 2 Freq Check
ts = Iterators.product(t, t);
yt = [cos(t1)+3sin(t2)+4 for (t1,t2) in ts];
h = [0 0;1 0;0 1;1 1];

yf = AFT([yt[:] yt[:]], h, N, :t2f);
YT = AFT(yf, h, N, :f2t);

# * HSEL
C = 2;
Nhmax = 4;

h = HSEL(Nhmax, 1:C);
Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));

# * HARMONIC STIFFNESS
M = I(2);
K = [2 -1;-1 2];
D = 0.001.*K + 0.01.*M;

C = 2;
Nhmax = 3;

h = HSEL(Nhmax, 1:C)
ws = [1, π];
ws = ws[1:C];

E, dEdw = HARMONICSTIFFNESS(M, D, K, ws, h);

# * Forced Response
C = 1;
Nhmax = 3;

h = HSEL(Nhmax, 1:C);
Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));
Fl = zeros(Nhc, 1);
Fl[2] = 1.0;
Fl = kron(Fl, [1;0]);

Omr, V = eigen(K, collect(M));
Omr = sqrt.(Omr);

Ωs = LinRange(0.5, 2, 500);
As = zeros(ComplexF64, size(Ωs));

for i = eachindex(Ωs)
    local E, dEdw
    E, dEdw = HARMONICSTIFFNESS(M, D, K, [Ωs[i]], h);
    U = E\Fl;
    As[i] = U[3]+im*U[5];
end

fsz = 18;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr)
   scr = GLMakie.Screen();
end

ax = Axis(fig[1, 1],
          xlabel="Excitation Frequency (rad/s)",
          ylabel="Response Amplitude (m)",
          xscale=Makie.pseudolog10,
          yscale=Makie.pseudolog10);
lines!(ax, Ωs, abs.(As))

display(scr, fig);
