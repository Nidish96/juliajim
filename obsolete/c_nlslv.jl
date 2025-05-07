using NLsolve
using UnPack
using Plots
include("../ROUTINES/CONTINUATION.jl");

function RESFUN!(R, dRdA, dRdΩ, AΩ, pars)
    @unpack F, c, ω0, α = pars
    if R !== nothing
        R[:] .= c^2/4 - (F/(2ω0*AΩ[1]))^2 + (ω0-abs(AΩ[2])+3α*AΩ[1]^2/8ω0)^2;
    end
    if dRdA !== nothing
        dRdA[:] .= F^2/(2ω0^2*AΩ[1]^3) + 3α*abs(AΩ[2])*(ω0-abs(AΩ[2])+3α*AΩ[1]^2/8ω0)/2ω0;
    end
    if dRdΩ !== nothing
        dRdΩ[:] .= -2(ω0-abs(AΩ[2])+3α*AΩ[1]^2/8ω0)*sign(AΩ[2]);
    end
end

pars = (F = 0.1, c = 0.01, ω0 = 2.0, α = 0.1);

# Continuation Settings
Ω₀ = 1.0;
Ω₁ = 3.0;
Δs = 0.2;  # Step Size
Δsmax = 0.2;
Δsmin = 0.01;
maxPts = 130;
itopt = 5;

Up = [0.01; Ω₀];  # First point Predictor
tgt = [0.0; sign(Ω₁-Ω₀)];

# 1. Correct First Point
U, R, J, info =
    NSOLVE((R,J,aΩ) -> RESFUN!(R, view(J, :,1:1), view(J, :,2:2), aΩ, pars), Up; N=1);
# Setup Lists
Us = [U];
UPs = [Up];
TGTs = [tgt];

# 2. Estimate Tangent
jQR = qr(J');
tgt = jQR.Q[:,end];
Δ = (-1)^(length(U)-1)*prod(diag(jQR.R));
tgt .*= sign(Δ);

# 3. Tangent Predictor
Up = U + tgt*Δs;

@time begin
    for is ∈ 1:maxPts
        global U, R, J, info, jQR, tgt, Δ, Δs, Up
        # 4. Correct Predicted Point
        U, R, J, info =
            NSOLVE((R,J,aΩ) -> RESFUN!(R, view(J, :,1:1), view(J, :,2:2), aΩ, pars), Up; N=1);
        # Append to Lists
        push!(Us, U);
        push!(UPs, Up);
        push!(TGTs, tgt);

        if U[2]*sign(Ω₁-Ω₀) ≥ Ω₁*sign(Ω₁-Ω₀)
            break;
        end
        
        # 5. Estimate Tangent
        jQR = qr(J');
        tgt = jQR.Q[:,end];
        Δ = (-1)^(length(U)-1)*prod(diag(jQR.R));
        tgt .*= sign(Δ);
        #    tgt .*= sign(tgt'*TGTs[end]);

        println("$is $(info.itns) $(tgt[2])")    
        # 7. Step Length Adaptation
        Δs = max(min(Δs*itopt/info.itns, Δsmax), Δsmin);
        
        # 6. Tangent Predictor
        Up = U + tgt*Δs;
    end
end

Us = hcat(Us...);
UPs = hcat(UPs...);
TGTs = hcat(TGTs...);

# gr()
# pyplot()
# gaston()
plot(Us[2,:], Us[1,:], marker=:circle, ms=2, label="sol", overwrite_figure=true)
plot!(UPs[2,:], UPs[1,:], marker=:circle, ms=2, label="pred")
