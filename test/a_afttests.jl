# * Preamble
using Revise
using Test

includet("../src/HARMONIC.jl");

# * Single Frequency Tests
# ** AFT
N = 128;
t = (0:N-1)*2π/N;

h = collect(0:3)[:,:];
Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));

yfsp = randn(Nhc, 2);

yt = yfsp[1,:]' .+
    sum([yfsp[2+2(i-1),:]'.*cos.(i*t)+yfsp[3+2(i-1),:]'.*sin.(i*t) for i in 1:h[end]]);
yf = AFT(yt, h, N, :t2f);

YT = AFT(yf, h, N, :f2t);

# ** Fourier Series Evaluation
yte = FSEVAL(h, t, yf)[1];

@testset "AFT-1" begin
    @test yf ≈ yfsp
    @test YT ≈ yt
    @test yte ≈ yt
end

# ** Fourier Differentiation Matrix
yf = AFT(cos.(t)+4sin.(2t), h,N, :t2f);
ydf = AFT(-sin.(t)+8cos.(2t), h,N, :t2f);

@testset "DFOUR-1" begin
    @test DFOUR(h)[1]*yf ≈ ydf
end

# ** Fourier Products
h = collect(0:7)[:,:];
Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));
    
Ni = round(Int, Nhc/3);
No = Nhc-Ni;
yf1 = [randn(Ni, 1); zeros(No,1)];
yf2 = [randn(Ni, 1); zeros(No,1)];

ypf = AFT(prod(AFT([yf1 yf2], h,N, :f2t), dims=2), h,N, :t2f);

@testset "Fourier Product" begin
    @test PRODMAT_FOUR(yf1, h)*yf2 ≈ ypf
end

# * Multi Frequency Tests

# ** HSEL
C = 2;
Nhmax = 3;

h = HSEL(Nhmax, 1:C)
Nhc = sum(all(h.==0, dims=2) + 2any(h.!=0, dims=2));

href = [0  1  0   2  2  3  1  2   1  0  1  0;
        0  0  1  -1  0  0  1  1  -2  2  2  3]';

@testset "HSEL" begin
    @test h == href
end

# ** AFT - 2 Frequency Case
yf = randn(Nhc, 1);

yt = AFT(yf, h,N, :f2t);
YF = AFT(yt, h,N, :t2f);

tt = Iterators.product(t,t);
yte = FSEVAL(h, [[t1 for (t1,t2) in tt][:] [t2 for (t1,t2) in tt][:]], yf)[1];

@testset "AFT-2" begin
    @test YF ≈ yf
    @test yte ≈ yt
end

# ** AFT - 3 Frequency Case
C = 3;
h = HSEL(Nhmax, 1:C);
Nhc = sum(all(h.==0, dims=2) + 2any(h.!=0, dims=2));

yf = randn(Nhc, 1);

yt = AFT(yf, h,N, :f2t);
YF = AFT(yt, h,N, :t2f);

tt = Iterators.product(t,t,t);
yte = FSEVAL(h, [[t1 for (t1,t2,t3) in tt][:] [t2 for (t1,t2,t3) in tt][:] [t3 for (t1,t2,t3) in tt][:]], yf)[1];

@testset "AFT-3" begin
    @test YF ≈ yf
    @test yte ≈ yt
end
