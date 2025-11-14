
function [A,B,C,G,r_tilde,z_tilde,N_Grid,N_Grid_centre,r_kron,z_kron,Psi_mn,Tp_t,Tp_b,Tp_c,Tp_s] ...
         = SpectralGalerkin_ThermalMdl(p,N_basis,h_coeffs)
 
% ========================================================================
% SpectralGalerkin_ThermalMdl.m
% ------------------------------------------------------------------------

% Purpose:
% - This function implements a control-oriented (cylindrical) cell thermal 
%   model using the Chebyshev Spectral–Galerkin (CSG) method
% - The model features independent boundary controls for the top, bottom, 
%   core/jelly, and lateral surface of the cell, enabling direction-specific 
%   cooling actuation
% - The governing heat PDE (Eq. 1) is solved using the 'method of
%   boundary lifting', expressed as T = T_h + T_p, where T_h satisfies the
%   homogeneous boundary conditions (BCs) and T_p enforces the inhomogeneous
%   (time-varying) boundary effects
% - This transformation reduces the original PDE to a reduced-order system
%   of ODEs in continuous descriptor form:  G·ẋ = A·x + B·u + F·w (Eq. 38),
%   where: G — mass matrix, A — stiffness matrix, B — input (boundary) matrix
%          F — volumetric heat-generation (disturbance) vector

% Inputs:
% p        - struct with geometric, thermophysical and input parameters
% N_basis  - no. of basis functions per direction (r,z)
% h_coeffs - vector [h_t, h_b, h_c, h_s] of convection coefficients

% Ouputs:
% A, B, C, G      - descriptor state-space matrices (continuous time)
% r_tilde/z_tilde - radial/axial spatial chebyshev grid nodes (coordinates) in domain [-1:1]
% N_Grid/N_Grid_centre - no. of spatial grid points for numerical integration/grid point of centre node
% r_kron/z_kron        - radial/axial physical-space coordinates in kronecker form
% Psi_mn               - 2D matrix containing basis functions (each row is one function)
% Tp_t/Tp_b/Tp_c/Tp_s  - particular solutions for each cell boundary 

% Reference: 
% Thermal Modelling of Battery Cells for Optimal Tab and Surface Cooling Control

% ========================================================================

%% Define scaled parameters & BCs
% Scaling factors to map [R_in, R_out] -> [-1, 1] and [0, L] -> [-1, 1]
alpha = 2/p.R;        %r-direction (radial) scaling
beta  = 2/p.H;        %z-direction (axial) scaling

% Convection coefficients
% h_t  = h_coeffs(1) // h_b = h_coeffs(2)
% h_c  = h_coeffs(3) // h_s = h_coeffs(4)
% Define Robin BC parameters 
% r-direction 
ap =  h_coeffs(4);   bp = alpha*p.kr;  %r = R_out, lateral surface
am = -h_coeffs(3);   bm = alpha*p.kr;  %r = R_in,  core
% z-direction 
cp =  h_coeffs(1);   dp = beta*p.kz;   %z = L,     top
cm = -h_coeffs(2);   dm = beta*p.kz;   %z = 0,     btm    

%% Basis functions and quadrature setup
% Grid points for numerical integration (Clenshaw-Curtis)
N_Grid = 41;     %no. of grid points per axis (N = 40, so 41 pts --> 1681 2D pts)
N_Grid_centre = round(N_Grid/2);    %centre node's index 

% Get 1D radial/axial basis functions (Phi_m(r), Phi_n(z)) evaluated at grid points
% These functions already satisfy the homogeneous Robin BCs
[r_tilde, Phi_r] = Cheby_BasisFxn(N_Grid-1,am,bm,ap,bp,N_basis-1);  %Eq. 15a
[z_tilde, Phi_z] = Cheby_BasisFxn(N_Grid-1,cm,dm,cp,dp,N_basis-1);  %Eq. 15b

% Build 2D tensor-product basis functions Psi_mn, using Kronecker product
% Psi_mn(r,z) = Phi_m(r).Phi_n(z)
N_states = N_basis^2;                                 %no. of model states
Psi_mn   = zeros(N_states,(N_Grid)*(N_Grid));         %pre-allocate 
k = 0;
for m = 1:N_basis             %index for r basis fxn, Phi_m(r)
    for n = 1:N_basis         %index for z basis fxn, Phi_n(z)
        Psi_mn(k+1,:) = kron(Phi_r(m,:),Phi_z(n,:));  
        k = k+1;
    end
end

%% Parameters for numerical quadrature (integration) (Clenshaw-Curtis)
[~,w] = clenCurt_Quad(N_Grid-1);        %quadrature weights
wr    = repelem(w,1,N_Grid)';           %radial weights 
wz    = repmat(w,1,N_Grid)';            %axial weights

%% Chebyshev differentiation matrices (Dm) (Operators) 
[~,Dm] = chebdif(N_Grid,1);             %Dm 1st order in [-1,1]
D2m    = Dm^2;                          %Dm 2nd order in [-1,1]
% Scale operators to physical domain
D_r    = alpha*Dm;                      %∂/∂r̃,   radial (1st order) 
D2_r   = alpha^2*D2m;                   %∂²/∂r̃², radial (2nd order) 
D2_z   = beta^2*D2m;                    %∂²/∂z̃², axial (2nd order) 
I      = eye(size(D2m));                

%% Conversion to original coordinates (physical space)
% Map grid points from [-1, 1] back to [R_in, R_out] and [0, L]
r = (1+r_tilde+alpha*p.r1)/alpha;       %Eq. 5a, radial coordinates
z = (z_tilde+1)/beta;                   %Eq. 5b, axial coordinates                  
% Flattened (Kronecker) vectors of physical space
r_kron = repelem(r,1,N_Grid)';          %radial coordinate 
z_kron = repmat(z,1,N_Grid)';           %axial coordinate

%% Construct state matrices (Galerkin projection)
% This loop performs the Galerkin projection onto each test function psi_m
% (m-loop) against each basis function psi_n (n-loop)

% Initialise matrices
G = zeros(N_states,N_states);
A = zeros(N_states,N_states); 
B = zeros(N_states,5);          %B matrix including disturbance vector(heat generation)
ny = 4;                         %no. of outputs
C = zeros(ny,N_states);         %output matrix

% Assign matrix values
for m = 0:N_states-1                   %for test fxns
    psi_eta = Psi_mn(m+1,:)';          %get m-th 2D basis function η(r,z)

    % Apply strong-form spatial operator L(Phi_eta)
    psi_rr = (kron(D2_r,I)*psi_eta);   %1st term in 40b, d^2(phi_r*phi_z)/dr_tilde^2 
    psi_r  = (kron(D_r,I)*psi_eta);    %2nd term in 40b, d(phi_r*phi_z)/dr_tilde
    psi_zz = (kron(I,D2_z)*psi_eta);   %3rd term in 40b, d^2(phi_r*phi_z)/dz_tilde^2
    L_psi_eta = (p.kr*r_kron.*psi_rr + p.kr*psi_r + p.kz*r_kron.*psi_zz);

    for n = 0:N_states-1               %for basis fxns
        psi_mn = Psi_mn(n+1,:)';       %get n-th 2D basis function Phi_m(r)*Phi_n(z)
        G(m+1,n+1) = p.rho * p.cp * ...
                     sum(r_kron.* psi_eta.* psi_mn.* wr.* wz) ...
                     * 1/alpha * 1/beta;  %Eq.40a
        A(m+1,n+1) = sum(L_psi_eta.*psi_mn.*wr.*wz) ...
                     * 1/alpha * 1/beta;  %Eq.40b, strong-form
    end

    % Disturbance input, F 
    % Projection of heat source Q onto test function eta
    B(m+1,1) = sum(r_kron.*psi_eta.*wr.*wz) * 1/alpha * 1/beta;   %Eq.40h

    % C (output) matrix (midpoints), Eq.40f-g
    % y_i = C_i x = sum_mn x_mn psi_mn(r_i,z_i) --> C_(i,m) = psi_mn(r_i,z_i)
    % This 'samples' the m-th basis function at sensor 'i'
    psisq    = reshape(psi_eta, N_Grid, N_Grid);  
    C(1,m+1) = psisq(end, N_Grid_centre);      %mid z = 0,     (z_tilde = -1, r_tilde = 0),  btm
    C(2,m+1) = psisq(1,  N_Grid_centre);       %mid z = L,     (z_tilde = 1 , r_tilde = 0),  top
    C(3,m+1) = psisq(N_Grid_centre, end);      %mid r = R_in,  (z_tilde = 0 , r_tilde = -1), core
    C(4,m+1) = psisq(N_Grid_centre, 1);        %mid r = R_out, (z_tilde = 0 , r_tilde = 1),  lateral surf
   
end

%% Determine expansion coefficients for particular solutions (Tp)
% Section III-D-2

% Matrices Big_Phi Φ, Eq, 28
P_t = zeros(N_basis,N_basis);
P_b = zeros(N_basis,N_basis);
P_l = zeros(N_basis,N_basis);
P_r = zeros(N_basis,N_basis);
for k = 1:N_basis
    phi_rk = Phi_r(k,:);
    phi_zk = Phi_z(k,:);
    for i = 1:N_basis
        phi_ri = Phi_r(i,:);
        phi_zi = Phi_z(i,:);
        P_t(k,i) = sum(r.* phi_rk.* phi_ri.* w) * 1/alpha;    %Big_Phi_h
        P_b(k,i) = sum(r.* phi_rk.* phi_ri.* w) * 1/alpha;    %Big_Phi_h
        P_l(k,i) = sum(p.r1* phi_zk.* phi_zi.* w) * 1/beta;   %Big_Phi_v
        P_r(k,i) = sum(p.r2* phi_zk.* phi_zi.* w) * 1/beta;   %Big_Phi_v
    end
end

% Source term, S vectors, Eq. 28
S_t = zeros(N_basis,1);
S_b = zeros(N_basis,1);
S_l = zeros(N_basis,1);
S_r = zeros(N_basis,1);
for i = 1:N_basis
    phi_ri = Phi_r(i,:);
    phi_zi = Phi_z(i,:);
    S_t(i,1) = sum(r.* phi_ri.* w) * 1/alpha;       %S_h
    S_b(i,1) = sum(r.* phi_ri.* w) * 1/alpha;       %s_h
    S_l(i,1) = sum(p.r1* phi_zi.* w) * 1/beta;      %S_v
    S_r(i,1) = sum(p.r2* phi_zi.* w) * 1/beta;      %S_v
end

%s_1 - b_2
k1 = cp + dp;    
k2 = cp + 2*dp;  
k3 = dm - cm;    
k4 = cm - 2*dm;     
j1 = ap + bp;    
j2 = ap + 2*bp;  
j3 = bm - am;    
j4 = am - 2*bm;  
scal_k1 = (k1/(k1*k4 - k2*k3));  %t_1 
scal_k4 = (k4/(k1*k4 - k2*k3));  %b_2
scal_j1 = (j1/(j1*j4 - j2*j3));  %s_1
scal_j4 = (j4/(j1*j4 - j2*j3));  %c_2

% D vectors, Eq. 29-31  
% D1    = scal_j4*(P_l\(u_s*S_l) - (j2/j4)*(P_r\(u_c*S_r))); 
% d^s_1 = scal_j4*((P_l\S_l));             *u_s                 
% d^c_1 = scal_j4*-1*((j2/j4)*(P_r\S_r));  *u_c    

% D2    = scal_j1*(P_r\(u_c*S_r) - (j3/j1)*(P_l\(u_s*S_l)));
% d^s_2 = scal_j1*-1*((j3/j1)*(P_l\S_l));  *u_s
% d^c_2 = scal_j1*((P_r\S_r));             *u_c
 
% D3    = scal_k4*(P_b\(u_t*S_b) - (k2/k4)*(P_t\(u_b*S_t))); 
% d^t_1 = scal_k4*((P_b\s_b));             *u_t 
% d^b_1 = scal_k4*-1*((k2/k4)*(P_t\s_t));  *u_b

% D4    = scal_k1*(P_t\(u_b*S_t) - (k3/k1)*(P_b\(u_t*S_b))); 
% d^t_2 = scal_k1*-1*((k3/k1)*(P_b\S_b));  *u_t
% d^b_2 = scal_k1*((P_t\S_t));             *u_b 

%% Model decomposition for control: Section III-D-3 
% Eq. 31, d-vectors
% ds
ds_1 = scal_j4 * ((P_l\S_l));                 
ds_2 = scal_j1 * -1 * ((j3/j1) * (P_l\S_l));      
% dc
dc_1 = scal_j1 * ((P_r\S_r));                 
dc_2 = scal_j4 * -1 * ((j2/j4) * (P_r\S_r));      
% dt
dt_1 = scal_k4 * ((P_b\S_b));                     
dt_2 = scal_k1 * -1 * ((k3/k1) * (P_b\S_b));      
% db
db_1 = scal_k1 * ((P_t\S_t));                 
db_2 = scal_k4 * -1 * ((k2/k4) * (P_t\S_t));      

% Construct full 2D T_p(r,z) on the grid
% Eq. 34
Tp_s = (kron(r_tilde',Phi_z'*ds_1) + kron((r_tilde').^2,Phi_z'*ds_2)); 
Tp_c = (kron(r_tilde',Phi_z'*dc_2) + kron((r_tilde').^2,Phi_z'*dc_1)); 
Tp_t = (kron(Phi_r'*dt_1,z_tilde') + kron(Phi_r'*dt_2,(z_tilde').^2));    
Tp_b = (kron(Phi_r'*db_2,z_tilde') + kron(Phi_r'*db_1,(z_tilde').^2));  

%% Input matrix - B 
% Apply strong-form operator L to each T_p and project unto test fxns

for m = 0:N_states-1
    psi_eta   = Psi_mn(m+1,:)';            %test fxns, η(r,z)

    % L(T_p) operators
    L_Tp_t = p.kr*r_kron.*(kron(D2_r,I)*Tp_t) + p.kr.*(kron(D_r,I)*Tp_t) + ...
             p.kz*r_kron.*(kron(I,D2_z)*Tp_t);
    
    L_Tp_b = p.kr*r_kron.*(kron(D2_r,I)*Tp_b) + p.kr.*(kron(D_r,I)*Tp_b) + ...
             p.kz*r_kron.*(kron(I,D2_z)*Tp_b);
             
    L_Tp_c = p.kr*r_kron.*(kron(D2_r,I)*Tp_c) + p.kr.*(kron(D_r,I)*Tp_c) + ...
             p.kz*r_kron.*(kron(I,D2_z)*Tp_c);
             
    L_Tp_s = p.kr*r_kron.*(kron(D2_r,I)*Tp_s) + p.kr.*(kron(D_r,I)*Tp_s) + ...
             p.kz*r_kron.*(kron(I,D2_z)*Tp_s);

    % Projections unto psi_eta
    B(m+1,2) = sum(L_Tp_t.* psi_eta.* wr.* wz) * 1/alpha * 1/beta;  %Eq. 40d
    B(m+1,3) = sum(L_Tp_b.* psi_eta.* wr.* wz) * 1/alpha * 1/beta;  %Eq. 40e
    B(m+1,4) = sum(L_Tp_c.* psi_eta.* wr.* wz) * 1/alpha * 1/beta;  
    B(m+1,5) = sum(L_Tp_s.* psi_eta.* wr.* wz) * 1/alpha * 1/beta;  %Eq. 40c
end


end




