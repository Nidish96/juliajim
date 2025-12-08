# ```@meta
# CurrentModule = juliajim
# ```
#
# ## [Example D1: Deflation to Locate Multiple roots](@id ex_d1)
#
# It is sometimes the case that we are trying to solve for the solution of a system of equations that has multiple zeros where one is (or a few are) already known. Case in point in dynamical systems: computing bifurcated branches emanating from a "main" branch of solutions.
#
# We use the method of deflation to discourage the solver from converging to the previously known solution(s) by penalizing the residue by the distance to that solution. Let us suppose that $u_0$ is the previously known solution point and $R(u)$ is the residue function that needs to be solved ($u\in \mathbb{R}^n$, $R: \mathbb{R}\to \mathbb{R}$). We write the deflated residue
# ```math
# R_d(u) = \frac{1}{||u-u_0||_2^2} R(u),
# ```
# where $||x||_2$ denotes the 2-norm. The limit $u\to u_0$ may or may not exist for $R_d(u)$ based on the speed at which $R(u)$ goes to zero at $u_0$, and this may affect the behavior of $R_d(u)$. Notwithstanding this, we find using the 2-norm as above to be quite effective for most problems.
#
# $R_d(u)$ discourages the solver from approaching the known solution $u_0$, allowing us to converge to another solution if it exists and if our initial guess is close enough.
#
# [`juliajim`](@ref juliajim) implements [`DEFLATEDRES!`](@ref juliajim.MDOFUTILS.DEFLATEDRES!) which is a deflation wrapper routine that can be used on top of any existing residue function to achieve deflation.
#
# This example demonstrates this through a simple 2-state system of polynomial equations with two known roots:
# ```math
# R(u; p) = \begin{bmatrix} u_2\\-p u_1-u_2+\frac{u_1^3}{4} \end{bmatrix},
# ```
# which is solved by $(u_1,u_2)=(-2\sqrt{p},0), (0,0), (2\sqrt{p},0)$.
#
# ## Load Packages
using NonlinearSolve
using LinearAlgebra

using juliajim.MDOFUTILS

# ## Setup Residue function
function fres!(up, r=nothing, j=nothing, jp=nothing)
    if !(r === nothing)
        r[:] = [up[2];-up[3]*up[1]-up[2]+up[1]^3/4];
    end
    if !(j === nothing)
        j[:, :] = [0 1;-up[3]+3up[1]^2/4 -1];
    end
    if !(jp === nothing)
        jp[:] = [0;-1];
    end
end


# ## Obtain the first solution
par = 3.0;
u0 = [3.0,0];

fun1 = NonlinearFunction((r,u,p)->fres!([u;p], r),
    jac=(J,u,p)->fres!([u;p], nothing,J),
    paramjac=(Jp,u,p)->fres!([u;p], nothing,nothing,Jp));
prob1 = NonlinearProblem(fun1, u0, par);
sol1 = solve(prob1, show_trace=Val(true));

# This returns the solution
sol1
# which is exactly $(2\sqrt{p},0)$ as expected.

# ## Deflated Solution
uC = [sol1.u; par];

fun2 = NonlinearFunction((r,u,p)-> DEFLATEDRES!([u;p], uC, fres!; R=r),
    jac=(J,u,p)-> DEFLATEDRES!([u;p], uC, fres!; J=J),
    paramjac=(Jp,u,p)-> DEFLATEDRES!([u;p], uC, fres!; Jp=Jp));
prob2 = NonlinearProblem(fun2, u0, par);
sol2 = solve(prob2, show_trace=Val(true));

# This returns the trivial solution $(0,0)$ for the above initial guess:
sol2

# ## Adding multiple deflators
# [`DEFLATEDRES!`](@ref juliajim.MDOFUTILS.DEFLATEDRES!) allows specifying multiple solutions (in a list) for deflation. Let us specify both the solutions found above.

uCs = [[sol1.u;par], [sol2.u;par]];

fun3 = NonlinearFunction((r,u,p)-> DEFLATEDRES!([u;p], uCs, fres!, R=r),
    jac=(J,u,p)-> DEFLATEDRES!([u;p], uCs, fres!, J=J),
    paramjac=(Jp,u,p)-> DEFLATEDRES!([u;p], uCs, fres!, J=J));
prob3 = NonlinearProblem(fun3, u0, par);
sol3 = solve(prob3, show_trace=Val(true));

# Now we recover the solution $(-2\sqrt{p},0)$!
sol3

# Now isn't that cool? We've used exactly the same initial guess and have progressively recovered all the solutions of the system.

# ## What if we kept going?
# In the practical setting we are sometimes not destined to know a priori how many solutions to expect. So let us keep going to see what happens.

uCs = [[s.u;par] for s in [sol1,sol2,sol3]];
fun4 = NonlinearFunction((r,u,p)-> DEFLATEDRES!([u;p], uCs, fres!, R=r),
    jac=(J,u,p)-> DEFLATEDRES!([u;p], uCs, fres!, J=J),
    paramjac=(Jp,u,p)-> DEFLATEDRES!([u;p], uCs, fres!, J=J));
prob4 = NonlinearProblem(fun4, u0, par);
sol4 = solve(prob4, show_trace=Val(true))

# NonlinearSolve will try out its list of methods and then finally get stalled. While this implies the absence of another solution in this case (because we know this), in the general case this merely implies the absence of another solution in the basin of attraction of the provided initial guess.
# The general problem of numerically eliminating possible solutions is quite complex!
