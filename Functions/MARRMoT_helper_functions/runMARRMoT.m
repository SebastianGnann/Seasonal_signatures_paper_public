function [Qmod,AET] = runMARRMoT(...
    P,T,PET,Q,... % forcing
    model_name,model_parameters,initial_storages,... % model specifications
    doPlot,showStats) % plotting
% RUNMARRMOT run a MARRMoT model with a single parameter set (Knoben et
% al., 2019)
%
% INPUT
% P: precipitation [mm] (time,variable)
% T: temperature [Degrees C]
% PET: potential evaporation [mm]
% Q: streamflow [mm]
% model_name: name of model according to MARRMoT handbook
% model_parameters: parameters
% intitial_storages: initial values of storages
% doPlot: plot results?
% showStats: show water balance etc.
% 
% OUTPUT
% Qmod: modelled flow [mm]
% AET: actual evaporation [mm]
%
% References
% Knoben, W.J., Freer, J.E., Fowler, K.J., Peel, M.C. and Woods, R.A., 
% 2019. Modular Assessment of Rainfall–Runoff Models Toolbox (MARRMoT)
% v1.2: an open-source, extendable framework providing implementations of 
% 46 conceptual hydrologic models as continuous state-space formulations. 
% Geoscientific Model Development, 12(6), pp.2463-2480.
%
% ---
% 
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk
%
% Acknowledgements: Wouter Knoben

% specify input data
input_climatology.precip   = P(:,2);      % Daily data: P rate  [mm/d]
input_climatology.temp     = T(:,2);      % Daily data: mean T  [degree C]
input_climatology.pet      = PET(:,2);    % Daily data: Ep rate [mm/d]
input_climatology.delta_t  = 1; % default

% specify model
model = model_name; 

% parameter values
input_theta = model_parameters;

% initial storage values
input_s0 = initial_storages; % [mm]

% define the solver settings
input_solver.name              = 'createOdeApprox_IE';     % Use Implicit Euler to approximate ODE's
input_solver.resnorm_tolerance = 0.1;                      % Root-finding convergence tolerance
input_solver.resnorm_maxiter   = 6;                        % Maximum number of re-runs

if showStats
    % run the model and extract all outputs
    tic
    [output_ex,~,~,~] = ...                                          % Fluxes leaving the model: simulated flow (Q) and evaporation (Ea)
        feval(model,...                                        % Model function name
        input_climatology,...                                  % Time series of climatic fluxes in simulation period
        input_s0,...                                           % Initial storages
        input_theta,...                                        % Parameter values
        input_solver);                                         % Details of numerical time-stepping scheme
    toc
    %     output_in,...                                          % Internal model fluxes
    %     output_ss,....                                         % Internal storages
    %     output_waterbalance...                                 % Water balance check
    
else
    [output_ex] = ...                                          % Fluxes leaving the model: simulated flow (Q) and evaporation (Ea)
        feval(model,...                                        % Model function name
        input_climatology,...                                  % Time series of climatic fluxes in simulation period
        input_s0,...                                           % Initial storages
        input_theta,...                                        % Parameter values
        input_solver);                                         % Details of numerical time-stepping scheme
end

% analyze the outputs
t = Q(:,1);
Qmod  = output_ex.Q';
AET = output_ex.Ea';

% plotting
if doPlot
    figure;
    plot(t,Q(:,2))
    hold on
    plot(t,Qmod)
    xlabel('Year'); ylabel('Flow [mm]')
    datetick
    legend('Observed','Modelled')
end

end