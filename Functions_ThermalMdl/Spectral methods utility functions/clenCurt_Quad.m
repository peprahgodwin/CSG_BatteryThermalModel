function [x,w] = clenCurt_Quad(N_Grid)

% ========================================================================
% clenCurt_Quad.m
% ------------------------------------------------------------------------

% Purpose:
% - This function performs numerical quadrature (integration) using Clenshaw-Curtis.
%   This method finds the integral of a function on [-1,1] by
%   approximating it as the sum of f(x_j)*w_j
% - It is based on a fixed set of Chebyshev-Gauss-Lobatto nodes (x_j) and 
%   their corresponding weights (w_j)
% - This code implements the direct summation formula (e.g., Trefethen,
%   Eq. 16.4) to find the weights
% - This numerical method can be classed as the formula of optimal order 
%   based on the fixed set of Chebyshev nodes (xj) as opposed to the Gauss 
%   formula of optimal order based on optimally chosen nodes [R2 - p.126]. 
% - NB: Code taken directly from (R2 - p.128)

% Inputs:
% N_Grid - no. of spatial grid points for numerical integration

% Outputs:
% x - Chebyshev Guass-Lobatto points
% w - quadrature weights

% References: 
% R1 - Spectral methods algorithms, analysis and applications.
% R2 - Spectral methods in Matlab, Trefethen, 2000
% R3 - Simulation of ODE/PDE Models with MATLAB®, OCTAVE and SCILAB 

% Comments expanded by: Godwin K. Peprah, 2025   

% ========================================================================

%% Compute Chebyshev-Gauss-Lobatto (extreme) nodes (points) in [-1,1]
% xj = cos((j*pi)/N) with j = 0,1,...,N
  theta = pi*(0:N_Grid)'/N_Grid;   %angles
  x     = cos(theta);

  %% Compute quadrature weights
  w  = zeros(1,N_Grid+1); 
  ii = 2:N_Grid;                   %indices for interior points
  v  = ones(N_Grid-1,1);

  if mod(N_Grid,2) == 0  %Case N is even:
      % Set endpoint weights
    w(1)        = 1/(N_Grid^2-1); 
    w(N_Grid+1) = w(1);
    % Summation for interior weights
    for k = 1:N_Grid/2-1 
        v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); 
    end
    % Final term in the sum
    v = v - cos(N_Grid*theta(ii))/(N_Grid^2-1);
  else    %Case N is odd:
    % Set endpoint weights
    w(1)        = 1/N_Grid^2; 
    w(N_Grid+1) = w(1);
    % Summation for interior weights
    for k = 1:(N_Grid-1)/2 
        v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); 
    end
  end

  % Assign computed interior weights (scaled)
  % w_j = 2*v_j / N for j = 1, ..., N-1
  w(ii) = 2*v/N_Grid;

end



