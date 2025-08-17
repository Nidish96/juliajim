# This is juliajim Documentation

Welcome to the documentation page. 

!!! note "Quick-Glimpse Tutorial"
    This tutorial just offered a quick glimpse on Julia's built-in documentation system, make sure to read the docs for more.

```@docs
juliajim
jimgreet()
```
# Introduction
Welcome welcome!
## Conventions
  * All scripts are tagged by an alphabet. 
  * All routines are upper case. 
# Getting Started
If you're just here for the...
  * Fourier Utilities
  * Continuation Utilities
  * Harmonic Balance Utilities
  * Transient Simulation Utilities
  * All of the above
# Examples with Documentation
Few examples of increasing complexity. Find them in /examples/.
# Documentation
## Harmonic/Galerkin Utilities
```@docs
AFT
HSEL
FSEVAL
HARMONICSTIFFNESS
HARMONICSTIFFNESS!
DFOUR
DFOUR!
PRODMAT_FOUR
ACT
DCHEB
PRODMAT_CHEB
```
## Continuation Utilities
### Routines
```@docs
EXTRESFUN!
EXTRESFUN_scaled!
CONTINUATE
```
## MDOF Problem Utilities
### Structs
```@docs
NONLINEARITY
MDOFGEN
```
### Routines
```@docs
ADDNL
NLEVAL!
HBRESFUN!
HBRESFUN_A!
EPMCRESFUN!
NLFORCE
NEWMARKMARCH
```
### Struct(s)
```@docs
myNLSoln
```
