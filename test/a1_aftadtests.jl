using Base: Forward
using Revise
using Test
using ForwardDiff
using GLMakie
using FFTW
using LinearAlgebra
using BlockDiagonals

using juliajim.HARMONIC

# * Single Frequency Case
N = 128;
t = (0:N-1)*2π/N;

h = collect(0:3)[:,:];
Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));

yfsp = randn(Nhc, 2);

Jad = ForwardDiff.jacobian(U->AFT(AFT(U, h,N, :f2t).^3, h,N, :t2f), yfsp);

ut = AFT(yfsp, h,N, :f2t);
dfdut = 3ut.^2;
cst = AFT(1.0I(Nhc), h,N, :f2t);
Jan = BlockDiagonal([AFT(dfdut[:, i].*cst, h,N, :t2f) for i in 1:2]);

@testset "1D AFT Fwd AutoDiff" begin
    @test Jad ≈ Jan
end

# * Multi Frequency Tests
C = 2;
Nhmax = 3;

h = HSEL(Nhmax, 1:C)
Nhc = sum(all(h.==0, dims=2) + 2any(h.!=0, dims=2));

yfsp = randn(Nhc, 1);

Jad = ForwardDiff.jacobian(U->AFT(AFT(U, h,N, :f2t).^3, h,N, :t2f), yfsp);

ut = AFT(yfsp, h,N, :f2t);
dfdut = 3ut.^2;
cst = AFT(1.0I(Nhc), h,N, :f2t);
Jan = AFT(dfdut.*cst, h,N, :t2f);

@testset "2D AFT Fwd AutoDiff" begin
    @test Jad ≈ Jan
end
