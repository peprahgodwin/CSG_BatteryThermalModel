function [x, DM] = chebdif(N_diff, N_deri)

% ========================================================================
% chebdif.m
% ------------------------------------------------------------------------

% Purpose:
% - This function computes the Chebyshev-Gauss-Lobatto nodes (grid points) 
%   and spectral differentiation matrices D^(1), D^(2), ..., D^(N_deri) 
%   on those nodes

% Notes:
% - The code implements two strategies for enhanced accuracy suggested by
%   W. Don and S. Solomonoff (SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268, 1994):
%   (a) Use of trigonometric identities to avoid computing differences (x(k)-x(j))
%   (b) Use of the "flipping trick" for numerical stability
%
% - Note (May 2003, from original authors):
%   It may be slightly better not to implement (a) and (b). See
%   "Spectral Differencing with a Twist" by R. Baltensperger and
%   M.R. Trummer (SIAM J. Sci. Comp.)
%   This code implements strategies (a) and (b)

% Inputs:
% N_diff - no. of nodes or grid points (size of the N x N differentiation matrix)
%          This is N in the N-th order polynomial
% N_deri - no. of derivatives required (e.g., 2 for D^1 and D^2).
%          Must be 0 < N_deri <= N_diff-1

% Outputs:
% x  - Chebyshev Guass-Lobatto nodes on [-1, 1]
% DM - 3D array (N_diff x N_diff x N_deri) where:
%      DM(:,:,1) = 1st derivative matrix (D1)
%      DM(:,:,2) = 2nd derivative matrix (D2)
%      ...
%      DM(:,:,ell) = ell-th derivative matrix


% Original code by: J.A.C. Weideman, S.C. Reddy 1998
% Comments expanded by: Godwin K. Peprah, 2025  

% ========================================================================

%% Setup grid and indices
I = eye(N_diff);         %identity matrix    
L = logical(I);          %logical mask for diagonal elements (fast indexing)

% Indices for "flipping trick" (strategy b)
n1 = floor(N_diff/2); 
n2 = ceil(N_diff/2);     

% Index vector k = [0, 1, ..., N_diff-1]'
k  = (0:N_diff-1)';      
% Angles theta for Chebyshev-Gauss-Lobatto nodes
th = k*pi/(N_diff-1);

% Compute Chebyshev-Gauss-Lobatto nodes (x_k = cos(th_k))
% Note: The sin() formulation is used, which is equivalent to cos(th)
% x_k = sin(pi/2 - th_k) = cos(th_k)
% This is part of the Don & Solomonoff numerical stability strategy
x = sin(pi * (N_diff-1:-2:1-N_diff)' / (2*(N_diff-1)));

%% Compute DX = [x_j - x_k] using trig identity (strategy a)

% Create a matrix T where T(j,k) = th(j)/2
T = repmat(th/2,1,N_diff);   

% Compute matrix of differences DX(j,k) = x_j - x_k
% This uses the trig identity:
%   x_j - x_k = cos(th_j) - cos(th_k)
%             = -2*sin((th_j+th_k)/2)*sin((th_j-th_k)/2)
% Here: T'+ T = (th_k + th_j)/2 and T'- T = (th_k - th_j)/2
% The sign is flipped, so DX(j,k) = x_j - x_k.
DX = 2 * sin(T'+T) .* sin(T'-T);   % DX(j,k) = x_j - x_k

% Apply the "flipping trick" 
% This enforces the anti-symmetry DX(j,k) = -DX(k,j) to maintain numerical precision
DX = [DX(1:n1,:); 
     -flipud(fliplr(DX(1:n2,:)))];

% Set diagonal elements to 1 to avoid division by zero. These will be corrected to 0 later
DX(L) = ones(N_diff,1);

%% Compute differentiation matrices (D^1, D^2, ...)

% C is the matrix C(j,k) = (-1)^(j+k)
% (NB: toeplitz((-1).^k) is a fast way to compute this)
C = toeplitz((-1).^k);
% Scale C to C(j,k) = (c_j/c_k)*(-1)^(j+k)
% where c_j = 2 for j=0 or N_diff-1, and 1 otherwise
C(1,:)      = C(1,:)*2; 
C(N_diff,:) = C(N_diff,:)*2;     
C(:,1)      = C(:,1)/2;
C(:,N_diff) = C(:,N_diff)/2;

% Z is the matrix of 1/(x_j - x_k) for j ~= k
Z    = 1./DX;
Z(L) = zeros(N_diff,1);    %diagonals are 0

% Initialize D. It will be iteratively built
% D starts as D^(0) = I
D = eye(N_diff);       

% Pre-allocate the output 3D matrix
DM = zeros(N_diff, N_diff, N_deri);

% Loop to compute D^(1), D^(2), ..., D^(N_deri)
    for ell = 1:N_deri
        % Compute D^(ell) from D^(ell-1)
        % This is the core recurrence relation for off-diagonal elements
        % D = D^(ell-1) at the start of this line
        D = ell*Z.*(C.*repmat(diag(D),1,N_diff) - D);   

        % Correct the main diagonal elements
        % Property: sum(D, 1) = 0 (column sums are zero)
        % This line sets d_jj = -sum(d_kj for k ~= j)
        % NB: sum(D') is a row vector of the column sums of D
        % (with the incorrect diagonal from the line above).
        % This line replaces the diagonal with the correct values
        D(L) = -sum(D');                                   
        % Store the computed D^(ell) matrix
        DM(:,:,ell) = D;                                

end

end


