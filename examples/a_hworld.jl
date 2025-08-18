# # Example a: Introduction to the AFT Routines
# This example is intended to showcase the Alternating Frequency Time (AFT) routines in [`juliajim.HARMONIC`](@ref).
# ## PreambleLoad the necessary packages
using GLMakie
using LinearAlgebra
using SparseArrays
using Arpack

using Revise
using juliajim.HARMONIC

# ## Alternating Frequency Time Transformation

# ### AFT - Single Frequency Check (Periodic signal)

# Here we will use the AFT routine to do a Fourier transform of a time-domain signal, `yt`, to compute its Fourier coefficients, `yf`.
# The general form assumed for a Fourier series is:

# ```math
# y = a_0 + \sum_{n=1}^H a_n \cos n\tau + \sin n\tau
# ```

# where $\tau$ is the scaled time such that the signal $y(\tau)$ is $2\pi$-periodic. For instance, if the excitation frequency is $\Omega$ and physical time is $t$ (such that we have $2\pi/\Omega$ as the time period), the scaled time coordinate $\tau$ is defined as $\tau=\Omega t$.

# The AFT routine provides a convenience utility for transforming from a discrete time array to an array of the Fourier coefficients which is written as

# ```math
# \begin{bmatrix} a_0 & a_1 & b_1 & a_2 & b_2 & \dots \end{bmatrix}.
# ```

# \\[
# m \ddot{x} + c \dot{x} = 0 
# \\]

# The same routine may also be used to do the opposite transformation (frequency coefficients to time array).

# We first discretize time by dividing the domain $[0, 2\pi]$ into `128` parts.
N = 128;
t = (0:N-1)*2π/N;
## yt = cos.(t) + 3sin.(2t) .+ 4;
yt = cos.(t) .+ 2;
yt = [yt yt];

h = collect(0:3)[:,:];
Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));

yf = AFT(yt, h, N, :t2f);
YT = AFT(yf, h, N, :f2t);

# ## AFT - 2 Freq Check
ts = Iterators.product(t, t);
yt = [cos(t1)+3sin(t2)+4 for (t1,t2) in ts];
h = [0 0;1 0;0 1;1 1];

yf = AFT([yt[:] yt[:]], h, N, :t2f);
YT = AFT(yf, h, N, :f2t);

# ## HSEL
C = 2;
Nhmax = 4;

h = HSEL(Nhmax, 1:C);
Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));

# ## HARMONIC STIFFNESS
M = I(2);
K = [2 -1;-1 2];
D = 0.001.*K + 0.01.*M;

C = 2;
Nhmax = 3;

h = HSEL(Nhmax, 1:C)
ws = [1, π];
ws = ws[1:C];

E, dEdw = HARMONICSTIFFNESS(M, D, K, ws, h);

# ## Forced Response
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

set_theme!(theme_latexfonts())
fsz = 18;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr) && Makie.current_backend()==GLMakie # src
   scr = GLMakie.Screen(); # src
end # src

ax = Axis(fig[1, 1],
          xlabel="Excitation Frequency (rad/s)",
          ylabel="Response Amplitude (m)",
          xscale=Makie.pseudolog10,
          yscale=Makie.pseudolog10);
lines!(ax, Ωs, abs.(As))

if Makie.current_backend()==GLMakie #src
   display(scr, fig); #src
else #src
   display(fig) #src
end #src

fig #md
