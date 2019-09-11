function [equilibrium_storages] = calculateInitialStorages(...
    P,T,PET,Q,....
    model_name,model_parameters,initial_storages)
% CALCULATEINITIALSTORAGES calculate initial storage by looping over the
% first year and using the storage at the end of the year as a new input
% until that storage stabilises.
%
% INPUT
% P: precipitation [mm] (time,variable)
% T: temperature [Degrees C]
% PET: potential evaporation [mm]
% Q: streamflow [mm]
% model_name: name of model according to MARRMoT handbook
% model_parameters: parameters
% intitial_storages: initial values of storages
%
% OUTPUT
% equilibrium_storages: initial equilibrium storages
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
input_climatology_warmup.precip   = P(1:365,2);      % Daily data: P rate  [mm/d]
input_climatology_warmup.temp     = T(1:365,2);      % Daily data: mean T  [degree C]
input_climatology_warmup.pet      = PET(1:365,2);    % Daily data: Ep rate [mm/d]
input_climatology_warmup.delta_t  = 1; % default

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

% run the model and extract all outputs
[~, ~, output_ss] = ...
    feval(model,...                                        % Model function name
    input_climatology_warmup,...                           % Time series of climatic fluxes in simulation period
    input_s0,...                                           % Initial storages
    input_theta,...                                        % Parameter values
    input_solver);                                         % Details of numerical time-stepping scheme

store_names = fieldnames(output_ss);
input_ss_old = input_s0; % initial old storages
input_ss_new = input_s0; % updated new storages
for i = 1:length(store_names)
    input_ss_new(i) = output_ss.(store_names{i})(end); % make old storage new initial storage
end

max_iter = 20; % maximum amount of years run for warm up
iter = 1;

% loop over first year until we reach an approx. equilibrium or max. iter.
while iter<=max_iter && any(abs(input_ss_new - input_ss_old)./input_ss_new > 0.01)
    iter = iter + 1;  
    input_ss_old = input_ss_new; 
    
    [~, ~, output_ss] = ...
        feval(model,...                                        
        input_climatology_warmup,...                           
        input_ss_old,...                                      
        input_theta,...                                        
        input_solver);  
    
    for i = 1:length(store_names)
        input_ss_new(i) = output_ss.(store_names{i})(end); % make old storage new initial storage
    end
end

equilibrium_storages = input_ss_new;

end