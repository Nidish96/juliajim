```@meta
CurrentModule = juliajim
```

# Welcome to the juliajim Documentation!

Welcome to the documentation page. 

!!! note "Quick-Glimpse Tutorial"
    This tutorial just offered a quick glimpse on Julia's built-in documentation system, make sure to read the docs for more.

# Introduction
Hello world.
## Conventions
  * All scripts are tagged by an alphabet. 
  * All routines are upper case. 
# Getting Started
If you're just here for the...
  * Fourier Utilities (start by looking at [`AFT`](@ref))
  * Continuation Utilities (start by looking at [`CONTINUATE`](@ref))
  * Integrated nonlinear dynamics (start by looking at [`MDOFGEN`](@ref))
      * Harmonic Balance Utilities (start by looking at [`HBRESFUN!`](@ref))
      * Transient Simulation Utilities (start by looking at [NEWMARKMARCH](@ref))
      * All of the above.
# Examples with Documentation
Few examples of increasing complexity. Find them in *examples*.

```@example
using juliajim

N = 8
t = (0:N-1)*2pi/N;
y = cos.(t);
h = 0:3;
Yh = AFT(y, h,N, :t2f);
```
