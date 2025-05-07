# * Preamble
using Revise
using Test

includet("../src/HARMONIC.jl");

# * Try
N = 128;

t = 2(0:N-1)/(N-1).-1;

h = collect(0:8)[:,:];

yt = t;
yf = ACT(yt, h,N, :t2f);

@testset "ACT-1" begin
    @test ACT(yf, h,N, :f2t) ≈ collect(yt)
end

