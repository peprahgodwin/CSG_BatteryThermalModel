function [sg] = MdlSetup_ThermalMdl(p,sg)

% ========================================================================
% MdlSetup_ThermalMdl.m
% ------------------------------------------------------------------------

% Purpose:
% - This function creates the discrete-time state-space thermal model by 
%   calling the 'SpectralGalerkin_ThermalMdl' function for the
%   continuous-time matrices (G, A, B, C) and stores them in struct sg
% - Discretised system: x(k+1) = Ad*x(k) + Bd_con*u(k) + Fd*w(k)
%                        y(k)  = C*x(k) + Du
% - It contains the convection coefficients which can be varied depending
%   on the type of cooling
% - The function also calculates the initial conditions and the output offset Du
%
% Inputs:
% p  - struct containing parameters
% sg - struct containing model setup (N_basis, t_step)
%
% Outputs:
% sg - updated struct containing state-space matrices
%
% Reference: 
% Thermal Modelling of Battery Cells for Optimal Tab and Surface Cooling Control

% ========================================================================

%% Convection coefficients (h):

% - Convection coefficients h, are varied to show the type of convection. 
% - Convective forced liquid cooling via cooling plate with h = 400[W/(m^2K)] 
%   assumed for areas exposed to tab/surface cooling each and mild forced 
%   convection to ambient air with h = 30[W/(m^2K)], to areas not exposed 
%   to tab/surface cooling
h_t = 400;           %z = L, top           
h_b = 400;           %z = 0, btm/base       
h_c = 0;             %R = Rin, jelly roll, negligible cooling
h_s = 400;           %R = R_out, lateral surface 
h_coeffs = [h_t; h_b; h_c; h_s];

%% Build continuous and discrete State-Space thermal models  
sg.N_states = sg.N_basis^2;   %no. of states
sg.N_inputs = 5;              %no. of inputs [disturbance, u_t, u_b, u_c, u_s]   

% Get continuous-time descriptor system matrices (G*dx/dt = A*x + B*u + F*w)
[A,B,C,G,~,~,N_Grid,N_Grid_centre,~,~,Psi_mn,Tp_t,Tp_b,Tp_c,Tp_s] = ...
    SpectralGalerkin_ThermalMdl(p,sg.N_basis,h_coeffs);
 
% Discretise (ZOH) system to x(k+1) = Ad*x(k) + Bd*u(k) + Fw(k)
sg.Ad = expm((G\A) * sg.t_step);
sg.Bd = (G\A) \ ((sg.Ad - eye(size(sg.Ad))) * (G\B));
sg.C  = C;
% Split Bd into disturbance (heat gen) and controls (boundaries)
sg.Fd     = sg.Bd(:,1);
sg.Bd_con = sg.Bd(:,2:sg.N_inputs);

%% Calculate initial conditions
% Initial inputs are based on initial external temps (T_inf) at each edge
T_inf0_stb = [p.T_0; p.T_0; p.T_0; p.T_0];    %initial ambient temps
u_t     =  h_t * T_inf0_stb(1);               %u = h*T_inf
u_b     = -h_b * T_inf0_stb(2);
u_c     = -h_c * T_inf0_stb(3);
u_s     =  h_s * T_inf0_stb(4);
sg.u_th = [u_t; u_b; u_c; u_s];

% Initial particular solution Tp and feedthrough term Du
T_p    = [Tp_t, Tp_b, Tp_c, Tp_s] * sg.u_th; 
% 'Sample' T_p at same sensor locations as 'C'
T_esq  = reshape(T_p, N_Grid, N_Grid);  
sg.Du    = zeros(4,1);
sg.Du(1) = T_esq(end, N_Grid_centre);   %mid z = 0,     (z_tilde = -1, r_tilde = 0),  btm
sg.Du(2) = T_esq(1,   N_Grid_centre);   %mid z = L,     (z_tilde = 1 , r_tilde = 0),  top
sg.Du(3) = T_esq(N_Grid_centre, end);   %mid r = R_in,  (z_tilde = 0 , r_tilde = -1), core
sg.Du(4) = T_esq(N_Grid_centre, 1);     %mid r = R_out, (z_tilde = 0 , r_tilde = 1),  lateral surf 

% Initial state vector and output vector y0
% Eq. 9: T = T_h + T_p = p.T = sg.xt * Psi_mn + Tp, with xt = cmn 
% Ax = b --> x = A\b, sg.Psi_mn'sg.xt_0 = p.T_0 - T_p0 -->
sg.xt_0 = Psi_mn'\(p.T_0 - T_p);     %Eq. 36
sg.yt_0 = C*sg.xt_0 + sg.Du;         %Eq. 38b
        

end
    












