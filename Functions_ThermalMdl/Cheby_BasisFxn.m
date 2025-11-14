function [x,Phi_x] = Cheby_BasisFxn(N_Grid,am,bm,ap,bp,N_basis)

% ========================================================================
% Cheby_BasisFxn.m
% ------------------------------------------------------------------------

% Purpose:
% - This function constructs Chebyshev basis functions that satisfy the 
%   specified (homogeneous Robin - in this case) boundary conditions (BC) 
%   at x = -1 and x = 1, and evaluates them at Chebyshev–Gauss–Lobatto nodes

% - Homogenous Robin BCs:
%   See (R1 - Eq. 4.10, page 143)
%   at x = -1 (left):  am*u(-1) + bm*u'(-1) = 0   
%   at x =  1 (right): ap*u(1)  + bp*u'(1)  = 0  

% Inputs:
% N_Grid  - no. of spatial grid points (NB: just for the numerical quadrature)
% am,bm   - a_ = -h/k, b_ = 1 --> BC coefficients at x = -1
% ap,bp   - a+ =  h/k, b+ = 1 --> BC coefficients at x =  1
% N_basis - no. of basis functions in a direction

% Outputs:
% x     - Chebyshev–Gauss–Lobatto nodes in [-1,1]
% Phi_x - Chebyshev basis function

% References: 
% R1 - Spectral methods algorithms, analysis and applications.
% R2 - Spectral methods in Matlab, Trefethen, 2000
% R3 - Simulation of ODE/PDE Models with MATLAB®, OCTAVE and SCILAB 

% Comments expanded by: Godwin K. Peprah, 2025   

% ========================================================================

%% Chebyshev Gauss-Lobatto (extreme) points in [-1,1]
% (R2 - Eq. 5.2, p.42, R3 - p.257)
% xj = cos((j*pi)/N) with j = 0, 1, ..., N
theta = pi*((0:N_Grid)/N_Grid);      %angles
x     = cos(theta);   

%% Compute Chebyshev basis (interpolation) fxn
% Phi_k(x)   = Tk(x)+ak*Tk+1(x)+bk*Tk+2(x)     (R1 - Eq. 4.26 p.148)
% with Tk(x) = cos(k*theta)  (R1 - Eq. 3.210, p.106, R3 - p.259)

% Initialize output matrix
Phi_x = zeros(N_basis+1,N_Grid+1);           
% Check if any BCs are specified
if (am^2 + bm^2 + ap^2 + bp^2 == 0)
    % Case 1: No BCs specified (am = bm = ap = bp = 0)
    % Default to the standard Chebyshev polynomials Tk(x)
    % (These already satisfy homogeneous Neumann BCs u_{\pm} 1) = 0$)
    for k = 0:1:N_basis
        Phi_x(k+1,:) = cos(k*theta);
    end
else
    % Case 2: General Robin Boundary Conditions
    % Calculate ak and bk for each k to satisfy the BCs
    for k = 0:1:N_basis
        % Determinant (R1 - p.148)
        DETk = 2*ap*am + ((k+1)^2+(k+2)^2)*(am*bp-ap*bm) - ...
               2*bm*bp*(k+1)^2*(k+2)^2;   

        % Constants {ak,bk} (R1 - Lemma 4.3, Eq. 4,28, p.148)
        ak = 4*(k+1)*(ap*bm + am*bp)/DETk;
        bk = (-2*am*ap + (k^2 + (k+1)^2)*(ap*bm - am*bp) + ...
             2*bm*bp*k^2*(k+1)^2)/DETk;

        % Construct the basis function \Phi_k from T_k, T_{k+1}, T_{k+2}
        Phi_x(k+1,:) = cos(k*theta) + ak*cos((k+1)*theta) + bk*cos((k+2)*theta);
    end
end

end
