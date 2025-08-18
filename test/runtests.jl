using Test

@testset "Alltests" begin
    @testset "a_afttests" begin
        include("./a_afttests.jl");
    end

    @testset "b_acttests" begin
        include("./b_acttests.jl");
    end;

    @testset "a1_aftadtests" begin
        include("./a1_aftadtests.jl");
    end

    @testset "c_contintest" begin
        include("./c_contintest.jl");
    end
end
