function [p] = Inputs_ThermalMdl(p)

% ========================================================================
% Inputs_ThermalMdl.m
% ------------------------------------------------------------------------

% Purpose:
% - This function defines the simulation inputs: 
%   time vector p.t, the volumetric, heat generatio profile p.Q, 
%   and the initial cell temperature p.T_0, and stores them in struct p
% - It is structured to allow easy switching between different
%   input profiles (e.g., WLTP cycle, pulse tests) by changing
%   the 'simulation_case' 

% Inputs:
% p - struct containing thermal model parameters

% Outputs:
% p - updated struct with inputs:
%   .T_0 - initial uniform cell temperature 
%   .t   - time vector 
%   .Q   - volumetric heat generation vector 

% Reference: 
% Thermal Modelling of Battery Cells for Optimal Tab and Surface Cooling Control

% ========================================================================

%% Iniitial uniform temperatures at each edge, T(t = 0,r,z)
p.T_0 = 15;         %[deg C]  initial cell temp 

%% Select simulation case
% Valid options: 'wltp', 'pulse1','pulse2'

simulation_case = 'wltp';   %choose one: 'wltp' or 'pulse1' or 'pulse2'

switch simulation_case
    case 'wltp'
        %% WLTP drive cycle
        % Loads pre-processed time and heat generation data from a .mat file
        load('wltp_time_n_q.mat', 't_wltp', 'Q_wltp');
        p.t = t_wltp;
        p.Q = Q_wltp;

    case 'pulse1'
        %% Pulse power profile 1 (step changes)
        p.tfin = 3500;           %[s] end time
        p.t    = (0:1:p.tfin)';  %[s] total time span

        % Define stages
        p.q1 = 40;               %[W/m^3] heat for stage 1
        p.t1 = 400;              %[s] end time  of stage 1
        p.q2 = 80;               %[W/m^3] heat for stage 2
        p.t2 = 850;              %[s] end time  of stage 2
        p.q3 = 0;                %[W/m^3] heat for stage 3

        % Vectorized heat generation (logical indexing), much faster than 
        % a 'for' loop
        p.Q = zeros(size(p.t));                         %initialise
        p.Q(p.t <= p.t1) = p.q1;                        %stage 1
        p.Q(p.t > p.t1 & p.t <= (p.t1 + p.t2)) = p.q2;  %stage 2
        p.Q(p.t > (p.t1 + p.t2)) = p.q3;                %stage 3
        % or using for loop
        % for i = 1:length(p.t)
        %   p.Q(i) = p.q1*heaviside(p.t1-p.t(i)) ...         %stage 1
        %            + p.q2*heaviside(p.t(i)-p.t1)*...       %stage 2
        %            heaviside((p.t1+p.t2)-p.t(i)) ...       %stage 2 (cont.)
        %            + p.q3*heaviside(p.t(i)-(p.t1+p.t2));   %stage 3
        % end

    case 'pulse2'
        %% Pulse power profile 2 (high power pulse)
        p.tfin = 2000;           %[s] end time
        p.t    = (0:1:p.tfin)';  %[s] total time span

        % Define stages
        p.q1 = 50;              %[W/m^3] heat for stage 1
        p.t1 = 150;             %[s] end time  of stage 1
        p.q2 = 1000;            %[W/m^3] heat for stage 2 (pulse)
        p.t2 = 50;              %[s] duration  of stage 2 (pulse)
        p.q3 = 0;               %[W/m^3] heat for stage 3

        p.Q = zeros(size(p.t));                         %initialise
        p.Q(p.t <= p.t1) = p.q1;                        %stage 1
        p.Q(p.t > p.t1 & p.t <= (p.t1 + p.t2)) = p.q2;  %stage 2
        p.Q(p.t > (p.t1 + p.t2)) = p.q3;                %stage 3
        %or
        % for i = 1:length(p.t)
        %   p.Q(i) = p.q1*heaviside(p.t1-p.t(i)) ...         %stage 1
        %            + p.q2*heaviside(p.t(i)-p.t1)*...       %stage 2
        %            heaviside((p.t1+p.t2)-p.t(i)) ...       %stage 2 (cont.)
        %            + p.q3*heaviside(p.t(i)-(p.t1+p.t2));   %stage 3
        % end
    otherwise
        error("Unknown 'simulation_case': %s. \nCheck Inputs_ThermalMdl.m", ...
              simulation_case);
end


end

