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

# U = randn(size(h,1));
U = zeros(size(h,1));
U[4] = 1.0;
u = ACT(U, h,N, :f2t);
Ue = ACT(u, h,N, :t2f);
ue = ACT(Ue, h,N, :f2t);

@testset "ACT - 2" begin
    @test Ue ≈ U
end

# ** Analytical Check for f2t
x = cos.(range(0, π, N));

xx = Iterators.product(x, x);

hi = 4;
U = randn(size(h,1));
# U = zeros(size(h,1));
# U[hi] = 1.0;

uan = sum([[U[hi].*cos.(h[hi,1].*acos.(x1)+h[hi,2].*acos.(x2)) for (x1,x2) in xx]
           for hi in eachindex(h[:,1])]);

u = ACT(U, h,N, :f2t);

@testset "ACT - 2 - Analytical" begin
    @test u ≈ uan[:]
end
# ** Analytical
N = 8;

h = HSEL(Nhmax, [1.,1.]);
h = vcat(unique(eachrow(abs.(h)))'...);

x = cos.(range(0, π, N));

xx = Iterators.product(x, x);

hi = 4;
# U = randn(size(h,1));
U = zeros(size(h,1));
U[hi] = 1.0;

uan = sum([[U[hi].*cos.(h[hi,1].*acos.(x1)+h[hi,2].*acos.(x2)) for (x1,x2) in xx]
           for hi in eachindex(h[:,1])]);
Ue = ACT(uan[:], h,N, :t2f)
[U Ue]

# @test U ≈ Ue

yin = uan[:];
tmp0 = reshape(yin, rNs);
tmp1 = tmp0[[1:N; 2:N-1], :];
tmp2 = [tmp1 -tmp1[:, 2:N-1]];

Nt = 2N-2;
L = I(sum(all(h.==0, dims=2) + 2*any(h.!=0, dims=2)));
L = L[[1; 2:2:end], :];
exa = reshape(AFT(L'U, h,Nt, :f2t), Nt,Nt);

@test tmp2 ≈ exa
