# * Preamble
using Revise
using Test

includet("../src/HARMONIC.jl");

# * AFT - Single Frequency Case
N = 128;
t = (0:N-1)*2π/N;

h = collect(0:3)[:,:];
Nhc = sum(all(h.==0, dims=2) + 2*any(h .!= 0, dims=2));

yfsp = randn(Nhc, 2);

yt = yfsp[1,:]' .+
    sum([yfsp[2+2(i-1),:]'.*cos.(i*t)+yfsp[3+2(i-1),:]'.*sin.(i*t) for i in 1:h[end]]);
yf = AFT(yt, h, N, :t2f);

YT = AFT(yf, h, N, :f2t);

@testset "AFT-1" begin
    @test isapprox(yf, yfsp)
    @test isapprox(YT, yt)
end

# * HSEL
C = 2;
Nhmax = 3;

h = HSEL(Nhmax, 1:C)
Nhc = sum(all(h.==0, dims=2) + 2any(h.!=0, dims=2));

href = [0  1  0   2  2  3  1  2   1  0  1  0;
        0  0  1  -1  0  0  1  1  -2  2  2  3]';

@testset "HSEL" begin
    @test h == href
end


# * AFT - Multiple Frequency Case
yf = randn(Nhc, 1);

yt = AFT(yf, h,N, :f2t);
YF = AFT(yt, h,N, :t2f);


@testset "AFT-2" begin
    @test isapprox(YF, yf)
end

# * AFT - 2D
N = 8;
Nhmax = 3;
C = 2;

h = HSEL(Nhmax, 1:C);
Nhc = sum(all(h.==0, dims=2) + 2any(h.!=0, dims=2));

# yf = randn(Nhc, 1);
yf = zeros(Nhc, 1);
yf[rinds[findall(r->r in ([2,-1], [1,-2]), eachrow(h)).-1]] = [1,2];

yt = AFT(yf, h,N, :f2t);
YF = AFT(yt, h,N, :t2f);

@test isapprox(yf, YF)

zinds, hinds, rinds0, rinds, iinds = HINDS(1, h);

yff = [yf[rinds0]; yf[rinds].+1im*yf[iinds]];
YFF = [YF[rinds0]; YF[rinds].+1im*YF[iinds]];
