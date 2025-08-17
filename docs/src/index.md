```@meta
CurrentModule = juliajim
```

# Welcome to the juliajim Documentation!

Welcome to the documentation page. 

!!! note "Why should you be interested?"
    Here's my vision for the package: **I would like to use a singly specified dynamical system to subject it to different types of analyses. This includes transient simulations, steady-state periodic/quasi-periodic, etc.**
	
	If this is interesting to you, you may find value in the project. 

Note that this is still under heavy development. It is functional but the documentation needs a lot of work. Please reach out to me in case you're interested in contributing, share the vision, and have started loving Julia! :-)

# Introduction

This is a set of Julia routines to help with computational research in nonlinear structural dynamics. Most of the content is derived from my earlier efforts (functional but badly maintained) with MATLAB/Octave at [https://github.com/Nidish96/octave-jim](https://github.com/Nidish96/octave-jim).

## Why yet another package?

Each programming language/environment/specification you choose (MATLAB, Python, Julia, C/C++, etc.) has its own host of libraries that can do a lot of what juliajim sets out to do. Notable examples include [https://trilinos.github.io/](Trilinos) (C++) and the [https://sciml.ai/](SciML ecosystem) (Julia), which can do what juliajim sets out and (FAR!) more.

My vision for this tool is this: **I'd like to be able to specify a system once and conduct different types of analysis on it.** This mainly includes nonlinear static solves, transient time-marching (also shooting for periodic solutions), harmonic balance (for periodic and quasi periodic solutions), and variants therein (potentially more).

Although a lot of the existing tools come very close to acheiving this, they don't span nearly all the kinds of nonlinearities. For instance, it is presently quite non-trivial to use any existing code as is for simulating systems with something as simple as an elastic dry-friction element (a hysteretically saturated linear spring).

Furthermore, very few tools provide a lot of convenience for harmonic balance simulations. A notable exception to this is the MATLAB tool [https://github.com/maltekrack/NLvib](NLvib) (a very well designed set of minimal routines that work really well). The main reason that juliajim is not just a fork of NLvib is Julia - I'm convinced of the several advantages of Julia and want to ensure a useful package exists in Julia that I can throw myself behind. Furthermore, since I started using Julia MATLAB has started feeling really annoying (read: https://mateusaraujo.info/2024/04/03/matlab-is-dead-long-live-julia/ ).

### Other Projects to Look at
+ [https://github.com/maltekrack/NLvib](NLvib): A Matlab toolbox for harmonic balance, shooting, and continuation.
+ [https://github.com/tmd-lab/tmdsimpy](tmdsimpy): A Python package with very similar goals as juliajim.
+ [https://github.com/Nidish96/octave-jim](octave-jim): My old code - its a mess that works.
+ [https://github.com/QuantumEngineeredSystems/HarmonicBalance.jl](HarmonicBalance.jl): A Julia toolbox for harmonic balance. 

## Design Philosophy

I am taking care to ensure that the project is as small as possible in order to ensure that the routines developed here can be part of something bigger. Furthermore, the different parts of the package (Harmonic, Continuation) are written in order to be fully functional independently. 
[`MDOFGEN`](@ref) and its suite of routines (including frequency domain residues, time domain marchers, etc.) are my way of using the two main parts of the package to provide a functional interface for nonlinear dynamics research.

### Conventions
  * All scripts are tagged by an alphabet. 
  * All routines are upper case. 
  
# Getting Started
If you're just here for the...
  * Fourier Utilities, start by looking at [`AFT`](@ref).
  * Continuation Utilities, start by looking at [`CONTINUATE`](@ref).
  * Integrated nonlinear dynamics, start by looking at [`MDOFGEN`](@ref).
      * Harmonic Balance Utilities (start by looking at [`HBRESFUN!`](@ref))
      * Transient Simulation Utilities (start by looking at [NEWMARKMARCH](@ref))
      * All of the above.
# Examples with Documentation
Few examples of increasing complexity. Find them in *examples*.
1. [../../examples/a_hworld.jl](a_hworld.jl)
2. [../../examples/b_duffhb.jl](b_duffhbl.jl)
3. [../../examples/c_jenkhb.jl](c_jenkhb.jl)
4. [../../examples/d_mdofgen_instnl.jl](d_mdofgen_instnl.jl)
5. [../../examples/d1_mdofgen_hystnl.jl](d1_mdofgen_hystnl.jl)

And more!

## Example code snippet
```@example
using juliajim

N = 8
t = (0:N-1)*2pi/N;
y = cos.(t);
h = 0:3;
Yh = AFT(y, h,N, :t2f);
```

