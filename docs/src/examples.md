# Overview

The examples are organized in an order of increasing complexity.

## Example code snippet
```@example
using juliajim.HARMONIC

N = 8
t = (0:N-1)*2pi/N;
y = cos.(t);
h = 0:3;
Yh = AFT(y, h,N, :t2f);
```
