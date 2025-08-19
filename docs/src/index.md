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

## Installation

The development is not yet complete and the package is not registered either. But you can try out the package by installing it manually from Github as follows:
```julia
] add https://github.com/Nidish96/juliajim
```
or, with `Pkg`
```julia
using Pkg
Pkg.add(url="https://github.com/Nidish96/juliajim")
```
After this you can just do `using juliajim` and get cracking.

## Why yet another package?

Each programming language/environment/specification you choose (MATLAB, Python, Julia, C/C++, etc.) has its own host of libraries that can do a lot of what juliajim sets out to do. Notable examples include [Trilinos](https://trilinos.github.io/) (C++) and the [SciML ecosystem](https://sciml.ai/) (Julia), which can do what juliajim sets out and (FAR!) more.

My vision for this tool is this: **I'd like to be able to specify a system once and conduct different types of analysis on it.** This mainly includes nonlinear static solves, transient time-marching (also shooting for periodic solutions), harmonic balance (for periodic and quasi periodic solutions), and variants therein (potentially more).

Although a lot of the existing tools come very close to acheiving this, they don't span nearly all the kinds of nonlinearities. For instance, it is presently quite non-trivial to use any existing code as is for simulating systems with something as simple as an elastic dry-friction element (a hysteretically saturated linear spring).

Furthermore, very few tools provide a lot of convenience for harmonic balance simulations. A notable exception to this is the MATLAB tool [NLvib](https://github.com/maltekrack/NLvib) (a very well designed set of minimal routines that work really well). The main reason that juliajim is not just a fork of NLvib is Julia - I'm convinced of the several advantages of Julia and want to ensure a useful package exists in Julia that I can throw myself behind. Furthermore, since I started using Julia MATLAB has started feeling really annoying (read: https://mateusaraujo.info/2024/04/03/matlab-is-dead-long-live-julia/ ).

### Other Projects to Look at
+ [NLvib](https://github.com/maltekrack/NLvib): A Matlab toolbox for harmonic balance, shooting, and continuation.
+ [tmdsimpy](https://github.com/tmd-lab/tmdsimpy): A Python package with very similar goals as juliajim.
+ [octave-jim](https://github.com/Nidish96/octave-jim): My old code - its a mess that works.
+ [HarmonicBalance.jl](https://github.com/QuantumEngineeredSystems/HarmonicBalance.jl): A Julia toolbox for harmonic balance. 
+ [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl): A great set of continuation and bifurcation routines.

## Design Philosophy

I am taking care to ensure that the project is as small as possible in order to ensure that the routines developed here can be part of something bigger. Furthermore, the different parts of the package ([`HARMONIC`](@ref), [`CONTINUATION`](@ref)) are written in order to be fully functional independently. For instance, the HB routines in the former portion of the toolbox can be used to write residue that can be used with [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl) for doing numerical continuation if the provided routines in [`CONTINUATION`](@ref) are deemed insufficient. 

[`MDOFUTILS`](@ref) and its suite of routines (including frequency domain residues, time domain marchers, etc.) are my way of using the two main parts of the package to provide a functional interface for nonlinear dynamics research.

### Conventions
  * All scripts are tagged by an alphabet. 
  * All routines are upper case. 
  
# Getting Started
If you're just here for the...
  * Fourier Utilities, start by looking at [`AFT`](@ref juliajim.HARMONIC.AFT).
  * Continuation Utilities, start by looking at [`CONTINUATE`](@ref juliajim.CONTINUATION.CONTINUATE).
  * Integrated nonlinear dynamics, start by looking at [`MDOFUTILS`](@ref).
      * Harmonic Balance Utilities, start by looking at [`HBRESFUN!`](@ref juliajim.MDOFUTILS.HBRESFUN!)
      * Transient Simulation Utilities, start by looking at [`NEWMARKMARCH`](@ref juliajim.MDOFUTILS.NEWMARKMARCH)
      * All of the above.
