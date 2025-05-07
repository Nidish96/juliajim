#TODO: Validate for higher dimensions
# * Preamble
using Revise
using Test

includet("../src/HARMONIC.jl");

# * Single Dimensional Tests
# ** ACT - 1
h = collect((0:8));
N = 64;
x = cos.(range(0, π, N));

U = randn(length(h), 1);
u = ACT(U, h,N, :f2t);
Ue = ACT(u, h,N, :t2f);
ue = ACT(Ue, h,N, :f2t);

@testset "ACT-1" begin
    @test Ue ≈ U
    @test ue ≈ u
end

# ** Differentiation
D1 = DCHEB(h);

u = 1 .+ 2x .+ 3cos.(2acos.(x)) .+ 4cos.(3acos.(x));
duan = 2 .+ 3*2sin.(2acos.(x))./sqrt.(1 .-x.^2) .+ 4*3sin.(3acos.(x))./sqrt.(1 .-x.^2);

dun = ACT(D1*ACT(u, h,N, :t2f), h,N, :f2t);

@testset "DCHEB - 1" begin
    @test dun[2:end-1] ≈ duan[2:end-1]
end

# ** Product Matrix
Nh = 17;
Ni = round(Int, Nh/3);
No = Nh-Ni;
h = collect(0:Nh-1);

U = [randn(Ni, 2); zeros(No, 2)];

pUr = ACT(prod(ACT(U, h,N, :f2t), dims=2), h,N, :t2f);

pU = PRODMAT_CHEB(U[:,1], h)*U[:,2];

@testset "PROD_CHEB - 1" begin
    @test pU ≈ pUr
end

# * Multi-Dimensional Tests
Nhmax = 3;
h = HSEL(Nhmax, [1.,2.]);

N = 64;
x = cos.(range(0, π, N));

U = randn(size(h,1), 1);
u = ACT(U, h,N, :f2t);
Ue = ACT(u, h,N, :t2f)

@testset "ACT - 2" begin
    @test Ue ≈ U
end
