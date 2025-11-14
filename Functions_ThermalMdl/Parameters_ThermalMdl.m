function [p] = Parameters_ThermalMdl(p)

% ========================================================================
% Parameters_ThermalMdl.m
% ------------------------------------------------------------------------

% Purpose:
% - This function defines cell geometric parameters and thermo-physical 
%   properties for the Spectral-Galerkin thermal model, and stores them in 
%   struct p
% - The parameters are based on the large format LFP cylindrical cell 
%   defined in Section IV-B of R1 

% Inputs:
% p - struct for thermal model parameters

% Outputs:
% p - updated struct with parameters

% Reference: 
% Thermal Modelling of Battery Cells for Optimal Tab and Surface Cooling Control

% ========================================================================

%% Cylindrical cell geometric parameters
p.r1 = 0.004;                       %[m]   inner radius (jelly roll core)
p.r2 = 0.032;                       %[m]   outer radius (cell casing)   
p.z1 = 0;                           %[m]   z = 0, bottom 
p.z2 = 0.198;                       %[m]   z = L, top 
p.H  = p.z2 - p.z1;                 %[m]   cell height (L)
p.R  = p.r2 - p.r1;                 %[m]   thickness of jelly role (R_out-R_in) 
p.Vb = pi*(p.r2^2 - p.r1^2)*p.z2;   %[m^3] volume of cell, hollow cylinder (pi*radius^2*height) 

%% Thermo-physical properties
p.rho = 2118;                       %[kg m^-3]      density 
p.cp  = 795;                        %[J kg^-1 K^-1] specific heat capacity 
p.kr  = 0.666;                      %[W m^-1 K^-1]  radial thermal conductivity 
p.kz  = 66.6;                       %[W m^-1 K^-1]  axial thermal conductivity 


end
    


