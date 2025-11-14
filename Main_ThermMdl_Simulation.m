
% =========================================================================
% Main_ThermMdl_Simulation.m
% -------------------------------------------------------------------------

% Purpose:
% - This function is the main executable script for the battery thermal model
% - It simulates the 2D thermal dynamics of a cylindrical battery cell using 
%   a pre-built Spectral-Galerkin state-space model
%
% Script workflow:
% - Initialises workspace
% - Loads all physical parameters (p) and simulation inputs (p.t, p.Q)
% - Calls MdlSetup_ThermalMdl to build the discrete-time matrices and initial states
% - Runs a time-step simulation to get the temperature response
% - Plots the resulting temperatures at key sensor locations
%
% Reference: 
% Thermal Modelling of Battery Cells for Optimal Tab and Surface Cooling Control

% =========================================================================

%% Initialise workspace

clear; clc; close all;
addpath(genpath('./Functions_ThermalMdl'));    %add all subfolders

%% Load model parameters and inputs
p = struct;                    %create parameters struct
p = Parameters_ThermalMdl(p);  %append parameters to struct p
p = Inputs_ThermalMdl(p);      %append inputs to struct p

% Define model settings
sg.t_step  = 1;                %[s] time step for simulation
sg.N_basis = 5;                %no. of basis functions per direction (r,z)

%% Build State-Space thermal Model
sg = MdlSetup_ThermalMdl(p,sg);

%% Initialise and run simulation 
% Pre-allocate simulation vectors
sg.xt      = zeros(sg.N_states, length(p.t));       %state (coefficients)
sg.u       = zeros(size(sg.Bd_con,2),length(p.t));  %control inputs    
sg.w       = zeros(1, length(p.t));                 %disturbance 
sg.yt      = zeros(size(sg.C, 1), length(p.t));     %output (temperatures) 
% Initial xt and yt
sg.xt(:,1) = sg.xt_0;                               
sg.yt(:,1) = sg.yt_0;
% Run simulation loop 
for k = 1:length(p.t)-1
    sg.u(:,k)    = [sg.u_th(1); sg.u_th(2); sg.u_th(3); sg.u_th(4)];   %constant inputs
    sg.w(k)      = p.Q(k)/p.Vb;                    %volumetric heat gen
    sg.xt(:,k+1) = sg.Ad*sg.xt(:,k) + sg.Bd_con*sg.u(:,k) + sg.Fd*sg.w(k);  %x(k+1) = Adx(k) + Bdu(k) + Fdw(k)
    sg.yt(:,k+1) = sg.C*sg.xt(:,k+1) + sg.Du;      %y(k) = Cx(k) + Du(k) = T_h + T_p
end

%% Plot results
figure;
hold on;
plot(p.t, sg.yt(2,:),'r', 'LineWidth', 2);  %Top 
plot(p.t, sg.yt(1,:),'k--','LineWidth', 2); %Btm
plot(p.t, sg.yt(3,:), 'k','LineWidth', 2);  %Core
plot(p.t, sg.yt(4,:), 'b','LineWidth', 2);  %Surf
xlim([min(p.t) max(p.t)]);
xlabel('Time [s]', 'fontsize', 16, 'fontname', 'Times New Roman', 'Interpreter', 'latex');
ylabel('Temperature [$^\circ$C]', 'fontsize', 16, 'fontname', 'Times New Roman', 'Interpreter', 'latex');
title('Battery (Wall) temperatures', 'fontsize', 16, 'fontname', 'Times New Roman', 'Interpreter', 'latex');
ax = gca;  
set(ax, 'FontSize', 15);
ax.Box = 'on';
legend('$T_{t,\mathrm{mid}}$', '$T_{b,\mathrm{mid}}$', '$T_{c,\mathrm{mid}}$','$T_{s,\mathrm{mid}}$','FontSize', 14, 'Location', 'best', 'fontname', 'Times New Roman', 'Interpreter', 'latex');
hold off
grid on;

