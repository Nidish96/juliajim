# ```@meta
# CurrentModule = juliajim.HARMONIC
# ```

# # [Example A: Introduction to the AFT Routines](@id ex_a)
# This example is intended to showcase the Alternating Frequency Time (AFT) routines in [`juliajim.HARMONIC`](@ref).

# At its core, it is just a bunch of utility routines that allow time-to-frequency and frequency-to-time transformations. These can also be used for Chebyshev polynomial expansions, as will be seen below.

# ## Preamble: Load the necessary packages
using GLMakie
using LinearAlgebra 
using SparseArrays
using Arpack

using Revise #src
using juliajim.HARMONIC

# ## Alternating Frequency Time Transformation

# ### AFT - Single Frequency (Periodic signal)

# Here we will use the AFT routine to do a Fourier transform of a time-domain signal, `yt`, to compute its Fourier coefficients, `yf`.
# The general form assumed for a Fourier series is:

# ```math
# y = a_0 + \sum_{n=1}^H a_n \cos n\tau + \sin n\tau
# ```

# where \(\tau\) is the scaled time such that the signal $y(\tau)$ is $2\pi$-periodic. For instance, if the excitation frequency is $\Omega$ and physical time is $t$ (such that we have $2\pi/\Omega$ as the time period), the scaled time coordinate $\tau$ is defined as $\tau=\Omega t$.

# The AFT routine provides a convenience utility for transforming from a discrete time array to an array of the Fourier coefficients which is written as

# ```math
# \begin{bmatrix} a_0 & a_1 & b_1 & a_2 & b_2 & \dots \end{bmatrix}.
# ```

# The same routine may also be used to do the opposite transformation (frequency coefficients to time array).

# We first discretize time by dividing the domain $[0, 2\pi]$ into `128` parts.
N = 128;
t = (0:N-1)*2π/N;
## yt = cos.(t) + 3sin.(2t) .+ 4;
yt = cos.(t) .+ 2;
yt = [yt yt];
# We have put two copies of the vector `yt` in columns to expose the fact that the [`AFT`](@ref) routine is vectorized: Having time series in `M` columns will return corresponding Fourier coefficients in `M` columns. 

# Now we set the list of harmonics of interest. This can also be an unordered list, but if you're including the zeroth harmonic, that should always be the first. 
h = 0:3
Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));
# The variable `Nhc` is the number of harmonic coefficients. For the zeroth harmonics we have only one coefficient (\(a_0\) above) and for every non-zero harmonic we have two, correponding to the sine and cosine (\(a_i,b_i\) above). The formula above looks more complicated than it must because this works even for the multi-frequency case. 

# Now we call the actual [`AFT`](@ref) routine. The first argument is the time series vector (or matrix of column vectors), the second argument is the harmonic indices of interest, third is the number of AFT samples (when inputting the time vector `N` has to be the same as the number of rows of `yt`), and the last argument is a `Symbol`. The last argument specifies if we're interested in a time-to-frequency (`:t2f`) or frequency-to-time (`:f2t`) transformation. 
yf = AFT(yt, h, N, :t2f);  # Time to Frequency
YT = AFT(yf, h, N, :f2t);  # Frequency to Time

# `yf` is the list of harmonics and `YT` is the time vector reconstructed from the harmonics. You should be able to verify that `yt` and `YT` are numerically the same.

# ### AFT - 2 Frequency Case (Quasi-periodic signal)

# Now we show how the same [`AFT`](@ref) routine can be used for multi-frequency Fourier representations. Our Fourier series representation is written as

# ```math
# y(t) = a_0 + \sum_{n=1}^H a_n cos( h1_n \tau_1 + h2_n \tau_2) + b_n sin( h1_n \tau_1 + h2_n\tau_2 ),\,\text{with } \tau_i = \Omega_i t.
# ```
# Here, \(h1, h2\) are index arrays storing the harmonic coefficient corresponding to the two frequencies present. The scaled time coordinates (aka Torus coordinates) \(\tau_1,\tau_2\) are defined as above such that the signal `y(t)` can be written as the torus function `Y(\tau_1,\tau_2)`, which is periodic on the 2D domain. Note that $y(t)$ and $Y(\tau_1, \tau_2)$ are not the same, although we can reconstruct \(y(t)\) from \(Y(\tau_1,\tau_2)\). 

# The usage of the routine is identical to before.
ts = Iterators.product(t, t);
yt = [cos(t1)+3sin(t2)+4 for (t1,t2) in ts];
h = [0 0;1 0;0 1;1 1];

yf = AFT([yt[:] yt[:]], h, N, :t2f);
YT = AFT(yf, h, N, :f2t);
yf

# #### Convenience Routine for Harmonic Selection: [`HSEL`](@ref)

# Since harmonic selection in the 2D (and general N-D) case is not as trivial as writing `0:3` for the 1D case (different combinations of the indices must be considered on the tensor-grid while accounting for redundancies. The convenience routine [`HSEL`](@ref) helps to do this - one can specify the maximum harmonic order (`Nhmax` below) and the number of components (`C` below) and obtain a \(H\times C\) matrix of indices.
C = 2;
Nhmax = 4;

h = HSEL(Nhmax, 1:C);  # Second argument is the list of frequencies
Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));

# #### Convenient Routine for Harmonic Indices: [`HINDS`](@ref)

# It is often necessary to also know the indices of the cosine and sine harmonics individually. The [`HINDS`](@ref) routine helps with this.
h = (0:4);
zinds, hinds, rinds0, rinds, iinds = HINDS(1, h)

# Here, `zinds` and `hinds` store the indices that have to be used for a complex representation (\(\begin{bmatrix} a_0&a_1-ib_1&a_2-ib_2&\dots \end{bmatrix}\)) and `rinds0`, `rinds,` and `iinds` store the indices \hat{t} have to be used for the real Fourier coefficients representation. Specifically, `rinds0` are the indices for the zeroth harmonic, `rinds` are the indices for the cosine harmonics and `iinds` are the indices for the sine harmonics.

# The routine can also be called when there are multiple Degrees-of-Freedom. See below for the case of 3 DoFs. This should also clarify the ordering convention used in this package.
h = (0:4);
zinds, hinds, rinds0, rinds, iinds = HINDS(3, h)

# Note that the routine also supports multi-frequency cases (`h` with multiple columns). 

# ## The "Harmonic" Stiffness

# Harmonic Balance is a numerical technique where the solution for a dynamical system (expressed in second order form) is expressed in terms of its Fourier series. The resulting algebraic system is solved numerically. For a linear system, however, the algebraic system can be simplified significantly and written in a succinct \(\mx{E} \vc{u} = \vc{f}\) form.

# Consider a linear dynamical system with \(n\) Degrees-of-Freedom of the form:

# ```math
# \mx{M} \ddot{\vc{x}} + \mx{C} \dot{\vc{x}} + \mx{K} \vc{x} = \vc{F}(t),\quad \vc{x}\in\mathbb{R}^n.
# ```

# We expand out the vector of unknowns \(\vc{x}\) as

# ```math
# \vc{x} = \vc{a}_0 + \sum_{n=1}^H \vc{a}_n \cos n\tau + \vc{b}_n \sin n\tau,\,\text{with } \vc{a}_i,\vc{b}_i\in \mathbb{R}^n,
# ```

# and write the vector of harmonic coefficients \(\vc{u}\) as

# ```math
# \vc{u} = \begin{bmatrix} \vc{a}_0\\ \vc{a}_1\\ \vc{b}_1\\ \vc{a}_2\\ \vc{b}_2\\ \vdots \end{bmatrix}.
# ```

# Substituting the harmonic ansatz into the governing equations and projecting the resulting residue onto the Fourier basis functions leads to:

# ```math
# \begin{bmatrix}
# \mx{K} & \\
# &\mx{K}-\Omega^2 \mx{M} & \Omega \mx{C} & \\
# &-\Omega \mx{C} & \mx{K}-\Omega^2 \mx{M} \\
# &&&\mx{K}-(2\Omega)^2 \mx{M} & (2\Omega) \mx{C}\\
# &&& -(2\Omega) \mx{C} & \mx{K}-(2\Omega)^2 \mx{M} \\
# &&&&&\ddots&&&&
# \end{bmatrix}
# \begin{bmatrix} \vc{a}_0\\ \vc{a}_1\\ \vc{b}_1\\ \vc{a}_2\\ \vc{b}_2\\ \vdots \end{bmatrix} =
# \begin{bmatrix} \vc{f}_{a0}\\ \vc{f}_{a1}\\ \vc{f}_{b1}\\ \vc{f}_{a2}\\ \vc{f}_{b2}\\ \vdots \end{bmatrix}.
# ```

# The [`HARMONICSTIFFNESS`](@ref) function provides a convenience routine for just this. It works for arbitrary number of frequency components also. The usage is as follows. Since the frequency jacobian of this matrix is often necessary, this is also returned by the routine. 
# There is also an in-place version of this in [`HARMONICSTIFFNESS!`](@ref).
M = I(2);
K = [2 -1;-1 2];
D = 0.001.*K + 0.01.*M;

C = 2;
Nhmax = 3;

h = HSEL(Nhmax, 1:C)
ws = [1, π];
ws = ws[1:C];

E, dEdw = HARMONICSTIFFNESS(M, D, K, ws, h);
E

# ## Forced Response Computation with [`HARMONIC`](@ref)
# We shall now use the above utilites to compute the forced response of a linear system. The code below should be self explanatory.
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
if !isdefined(Main, :scr) && Makie.current_backend()==GLMakie #src
    scr = GLMakie.Screen(); #src
end #src

ax = Axis(fig[1, 1],
    xlabel="Excitation Frequency (rad/s)",
    ylabel="Response Amplitude (m)",
    xscale=Makie.pseudolog10,
    yscale=Makie.pseudolog10);
lines!(ax, Ωs, abs.(As))

if Makie.current_backend()==GLMakie #src
    display(scr, fig); #src
else #src
    fig
end #src
