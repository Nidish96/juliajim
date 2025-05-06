using LinearAlgebra
using NLsolve
using SparseArrays
using Arpack
using Printf

function SOLVECONT(func!, x0, lam0, lam1, ds; itopt=5,
                   dsmin=nothing, dsmax=nothing,
                   stepmax=100, Dscale=nothing,
                   DynDscale=false)
    Nx = length(x0);
    if dsmin === nothing
        dsmin = ds/5;
    end
    if dsmax === nothing
        dsmax = ds*5;
    end
    if DynDscale && Dscale === nothing
        Dscale = ones(Nx+1,1);
    end
    
    # Initial Point (solving this without scaling)
    res = nlsolve(only_fj!((F,J,x)->func!(F,J,nothing,[x;lam0])),
                  x0, show_trace=true);
    x0 = res.zero;
    J0 = zeros(Nx, Nx+1);
    func!(nothing, view(J0, :, 1:Nx), view(J0, :, Nx+1), [x0;lam0]);
    if Dscale !== nothing
        J0 .*= Dscale';
    end
    
    # Initialize Tangent (in scaled space)
    z = -J0[:, 1:Nx]\J0[:, Nx+1];
    dxn = sign(lam1-lam0);
    α = dxn/sqrt(1.0 + z'*z)
    pars = Dict("tgt"=> [z;1] .* α, "xl0"=>[x0;lam0]);
    if Dscale !== nothing
        pars["xl0"] ./= Dscale;
    end
    
    # Initiate Continuation
    xl = [[x0;lam0]];  # Storing initial solution (physical space)
    J = zeros(Nx, Nx+1);

    println("Starting Continuation");

    println("===========================================");
    if Dscale === nothing
        println("λ, ds, Itns");
    else
        println("λˢ (λ), ds, Itns");
    end
    while (xl[end][end]-lam1)*(lam1-lam0) < 0
        if DynDscale
            pars["xl0"] .*= Dscale;
            pars["tgt"] .*= Dscale;
            Dscale = max.(abs.(pars["xl0"]), Dscale);
            pars["xl0"] ./= Dscale;
            pars["tgt"] ./= Dscale;
        end
        # First order initial guess (in scaled space)
        xlg = pars["xl0"] + pars["tgt"]*ds;
        
        res = nlsolve(only_fj!((F,J,xl)->EXTRES!(F,J,xl, func!,ds,pars,Dscale=Dscale)),
                      xlg, show_trace=false);
        
        push!(xl, res.zero);
        if Dscale !== nothing
            xl[end] = xl[end] .* Dscale;
        end
        
        # Compute Tangent (scaled space)
        α0 = α;
        z0 = z;
        func!(nothing, view(J, :, 1:Nx), view(J, :, Nx+1), xl[end]);
        if Dscale !== nothing
            J .*= Dscale';
        end
        z = -J[:, 1:Nx]\J[:, Nx+1];
        dxn = sign(α0*(z'*z0+1.0));
        α = dxn/sqrt(1+z'*z);
        pars["tgt"] = [z;1]*α;
        pars["xl0"] = copy(xl[end]);
        if Dscale !== nothing
            pars["xl0"] ./= Dscale;
        end
        
        # Print Line to Screen
        if Dscale === nothing
            @printf("%.2f, %.2f, %d, %d\n", res.zero[end], ds, res.iterations, dxn);
        else
            @printf("%.2f (%.2f), %.2f, %d, %d\n", res.zero[end], xl[end][end], ds,
                    res.iterations, dxn);
        end

        # Adapt Step Size & Set Initial guess
        ds = max(min(ds*(itopt/res.iterations)^2, dsmax), dsmin);
        
        if length(xl) >= stepmax
            break
        end
    end
    
    return (hcat(xl...));
end

function EXTRES!(R, dRdxl, xl, func!, ds, pars; Dscale=nothing)
    Nx = length(xl)-1;
    if Dscale !== nothing
        xl .*= Dscale;
    end
    if R !== nothing && dRdxl !== nothing
        func!(view(R, 1:Nx), view(dRdxl, 1:Nx,1:Nx), view(dRdxl, 1:Nx, Nx+1), xl);
    elseif dRdxl === nothing
        func!(view(R, 1:Nx), nothing, nothing, xl);
    elseif R === nothing 
        func!(nothing, view(dRdxl, 1:Nx,1:Nx), view(dRdxl, 1:Nx, Nx+1), xl);
    end
    if Dscale !== nothing
        xl ./= Dscale;
        if dRdxl !== nothing
            dRdxl .*= Dscale';
        end
    end

    # Orthonormal Constraint
    if R !== nothing
        R[end] = (pars["tgt"]'*(xl-pars["xl0"]))[] - ds;
    end
    if dRdxl !== nothing
        dRdxl[end, :] = pars["tgt"]';
    end
end

function NSOLVE(func!, U0; N=nothing, maxIter=20, tol=1e-10, show_trace::Bool=false)
    if N === nothing
        N = length(U0);
    end
    n = length(U0);
    
    R = zeros(N);
    J = zeros(N,n);
    func!(R, J, U0);
    Delta_U = -J\R;
    errs = Dict(:u => norm(Delta_U), :r => norm(R));
    its = 0;
    if show_trace
        println("Iter\t norm(Delta_U)\t norm(R)");
        println("----\t --------\t -------");
        println("$its\t $(errs[:u])\t $(errs[:r])");
    end
    U = U0;
    while its == 0 || (errs[:r]>tol && its<maxIter)
        U += Delta_U;
        func!(R, J, U);
        Delta_U[:] = -J\R;
        its += 1;
        errs[:u] = norm(Delta_U);
        errs[:r] = norm(R);
        if show_trace
            println("$its\t $(errs[:u])\t $(errs[:r])");
        end
    end
    info = (itns = its, );
    return (U, R, J, info);
end
