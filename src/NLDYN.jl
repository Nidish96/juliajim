# * Preamble
using Random
using LinearAlgebra
using TensorOperations

# * Exports
export NORMALFORMFIT

# * Normal Form Fitting
"""
# Description
NORMALFORMFIT fits a quadratic normal form to the given residue
function (in R2) numerically. It starts with a cloud of randomly
distributed points with given amplitude standard deviation and
fits a general cubic model such that the dynamics may be
approximated as:
       xdot = A x + B x x + C x x x .
  Then the quadratic terms are eliminated by posing a near
identity transformation from x to z as
	x = z + H z z .
  The simplified normal form equations in terms of z are
	zdot = A z + C3 z z z .

  A good reference for this is sec. 19.10 in Wiggins (1990).   
# Arguments
- xyfun       : (xy) Taking in a 2-vector and returning a 2-residue vector.
		Must have equilibrium at origin.
- d_amp       : If scalar, standard deviation for amplitude.
                If vector, [d_min, d_max] for the amplitude.
- Npts        : Number of points to be used for fit.
- constrained : (default true)
# Outputs
- A           : (2,2) The linear Jacobian matrix.
- B2n         : (2,3) The quadratic coefficient matrix such that
                      the nonlinear influence is written as
                         B2*[z1^2; z1*z2; z2^2]
- C3n         : (2,4) The cubic coefficient matrix such that the
                      nonlinear influence is written as
                         C3*[z1^3; z1^2*z2; z1*z2^2; z2^3].
- Hn          : (2,3) The cubic near identity transformation
                      matrix such that the near identity
                      transformation can be written as
                         x = z + Hn [z1^2; z1*z2; z2^2]
- C3          : (2,2,2,2) The cubic coefficient matrix in tensor
                          form.
- H           : (2,2,2) The quadratic near identity
                        transformation coefficients in tensor form.
"""
function NORMALFORMFIT(xyfun, d_amp, Npts; constrained=true)
    if Npts<10
        @warn "For the linear regression to work as expected" *
              "(to have a full rank regressor matrix), a minimum" *
              "of 10 points is recommended."
    end
    if length(d_amp)==1
        xys = randn(2, Npts)*d_amp;
    else
        ds = rand(Npts).*diff(d_amp[1:2]) .+ d_amp[1];
        ths = rand(Npts)*2π;
        xys = (ds.*[cos.(ths) sin.(ths)])';  # 2, Npts
    end
    Rps = hcat(xyfun.(eachcol(xys))...)';  # Npts, 2 residue

    XX = hcat(ones(Npts), xys[1,:], xys[2,:],
        xys[1,:].^2, xys[1,:].*xys[2,:], xys[2,:].^2,
        xys[1,:].^3, xys[1,:].^2.0.*xys[2,:], xys[1,:].*xys[2,:].^2, xys[2,:].^3);
    xycofs = XX\Rps;

    A0 = xycofs[1,:]';
    A1 = xycofs[2:3,:]';
    A2 = xycofs[4:6,:]';
    A3 = xycofs[7:10,:]';
    As = [A0, A1, A2, A3];

    # * Coefficients in Tensor Notation
    A = A1;

    B = zeros(2,2,2);
    B[:,1,1] = A2[:,1];
    B[:,2,2] = A2[:,3];
    B[:,1,2] = A2[:,2]/2;
    B[:,2,1] = A2[:,2]/2;

    C = zeros(2,2,2,2);
    C[:,1,1,1] = A3[:, 1];
    C[:,1,1,2] = A3[:, 2]/3;
    C[:,1,2,1] = A3[:, 2]/3;
    C[:,2,1,1] = A3[:, 2]/3;
    C[:,1,2,2] = A3[:, 3]/3;
    C[:,2,1,2] = A3[:, 3]/3;
    C[:,2,2,1] = A3[:, 3]/3;
    C[:,2,2,2] = A3[:, 4];

    # * Estimate Near Identity Transformation H
    δ = I(2);
    O = zeros(2,2,2,2);
    @tensor begin
        O[i, j, p, q] = 2δ[p, i]*A[q, j] - δ[q, j]*A[i, p];
    end
    Omx = reshape(O, 4,4);
    Bmx = reshape(B, 4,2);
    if constrained
        L = I(6);
        L = L[[1:4;3:6], :];
        Hmx = reshape(L*((kron(I(2), Omx)*L)\Bmx[:]), 4,2);
    else
        Hmx = Omx\Bmx;
    end
    H = reshape(Hmx, 2,2,2);

    # * Obtain Transformed Coefficient Tensor
    B2 = reshape(Bmx-Omx*Hmx, 2,2,2);

    C3 = zeros(2,2,2,2);
    @tensor begin
        C3[i,j,k,l] = C[i,j,k,l] + 2B[i,j,p]*H[p,k,l]-
                      2H[i,p,j]*(A[p,q]*H[q,k,l]+B[p,k,l]);
    end

    # * Convert to a Palatable form for implementation
    Hn = zeros(2, 3);  # Coefs for x^2, xy, y^2
    Hn[:, 1] = H[:, 1,1];
    Hn[:, 2] = H[:, 1,2]+H[:, 2,1];
    Hn[:, 3] = H[:, 2,2];

    B2n = zeros(2, 3);  # Coefs for x^2, xy, y^2
    B2n[:, 1] = B2[:, 1,1];
    B2n[:, 2] = B2[:, 1,2]+B2[:, 2,1];
    B2n[:, 3] = B2[:, 2,2];

    C3n = zeros(2, 4);  # Coefs for x^3, x^2y, xy^2, y^3
    C3n[:, 1] = C3[:, 1,1,1];
    C3n[:, 2] = C3[:, 1,1,2] + C3[:, 1,2,1] + C3[:, 2,1,1];
    C3n[:, 3] = C3[:, 1,2,2] + C3[:, 2,1,2] + C3[:, 2,2,1];
    C3n[:, 4] = C3[:, 2,2,2];

    return A, B2n, C3n, Hn;
end
