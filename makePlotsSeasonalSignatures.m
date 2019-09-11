%% Plot results from seasonal signature analysis
%
%   - shows how to create the plots shown in the paper, with a few
%   examples.
%   - to reproduce all the plots shown in the paper, you will need some 
%   metadata which is not included here (see paper for data sources).
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

% clc
% clear all
% close all

%% add directories for functions to path

if exist('./BrewerMap') == 7
    addpath(genpath('./BrewerMap'));
else
    error('BrewerMap toolbox needed. Can be download from https://github.com/DrosteEffect/BrewerMap and should be in a folder named BrewerMap in the same directory.')
end

%% load data and results

% load results
results_UK = load('seasonal_signatures_UK'); % seasonal_signatures_UK_2
seasonal_signatures_UK = results_UK.seasonal_signatures_UK; 
results_CAMELS = load('seasonal_signatures_CAMELS'); % seasonal_signatures_CAMELS_2
seasonal_signatures_CAMELS = results_CAMELS.seasonal_signatures_CAMELS;

% store a few variables which are used very frequently
ID_UK = seasonal_signatures_UK.ID;
ID_CAMELS = seasonal_signatures_CAMELS.ID;

% phase shift, amplitude ratio and offset (mean)
A_vec_UK = seasonal_signatures_UK.A_vec; % change in amplitude
phi_vec_UK = seasonal_signatures_UK.phi_vec; % phase shift
ERR_vec_UK = seasonal_signatures_UK.ERR_vec; % proxy for mean(Q)/mean(P-PET)
A_vec_CAMELS = seasonal_signatures_CAMELS.A_vec; % change in amplitude 2nd method
phi_vec_CAMELS = seasonal_signatures_CAMELS.phi_vec; % phase shift 2nd method
ERR_vec_CAMELS = seasonal_signatures_CAMELS.ERR_vec; % proxy for mean(Q)/mean(P-PET)

% confidence intervals / uncertainty (assumed negligible when using method 1)
A_vec_confInt_UK = seasonal_signatures_UK.A_vec_confInt; % change in amplitude
phi_vec_confInt_UK = seasonal_signatures_UK.phi_vec_confInt; % phase shift
ERR_confInt_UK = seasonal_signatures_UK.ERR_vec_confInt; % proxy for mean(Q)/mean(P-PET)
A_vec_confInt_CAMELS = seasonal_signatures_CAMELS.A_vec_confInt; % change in amplitude
phi_vec_confInt_CAMELS = seasonal_signatures_CAMELS.phi_vec_confInt; % phase shift
ERR_confInt_CAMELS = seasonal_signatures_CAMELS.ERR_vec_confInt; % proxy for mean(Q)/mean(P-PET)

% annual period
period_1y = 365;
w_1y = 2*pi/period_1y;

%% plot results - UK

% scatter plots
% shows benchmark catchments and example catchments
plotSeasonalSignaturesDataBenchmark(A_vec_UK,phi_vec_UK,w_1y,...
    'attribute',seasonal_signatures_UK.isBenchmark,'attribute_name','isBenchmark',...
    'ID',ID_UK,... 
    'A_confInt',A_vec_confInt_UK,'phi_confInt',phi_vec_confInt_UK,...
    'x_limits',[0 1.2],'y_limits',[0 140],...
    'colour_scheme','Reds','flip_colour_scheme',false,...
    'c_limits',[-1 1],'c_lower_limit_open',false,'c_upper_limit_open',false,...
    'figure_title','(a)','figure_name','UK_isBenchmark')

% climate and catchment attributes
plotSeasonalSignaturesData(A_vec_UK,phi_vec_UK,w_1y,...
    'attribute',seasonal_signatures_UK.I_m,'attribute_name','I_m',...
    'ID',ID_UK,... 
    'A_confInt',A_vec_confInt_UK,'phi_confInt',phi_vec_confInt_UK,...
    'x_limits',[0 1.2],'y_limits',[0 140],...
    'colour_scheme','Reds','flip_colour_scheme',true,...
    'c_limits',[-1 1],'c_lower_limit_open',false,'c_upper_limit_open',true,...
    'figure_title','(a)','figure_name','UK_I_m')

%% plot results - CAMELS

% climate and catchment attributes
plotSeasonalSignaturesData(A_vec_CAMELS,phi_vec_CAMELS,w_1y,...
    'attribute',seasonal_signatures_CAMELS.I_m,'attribute_name','I_m',...
    'ID',ID_CAMELS,... 
    'A_confInt',A_vec_confInt_CAMELS,'phi_confInt',phi_vec_confInt_CAMELS,...
    'x_limits',[0 1.2],'y_limits',[0 250],...
    'colour_scheme','Reds','flip_colour_scheme',true,...
    'c_limits',[-1 1],'c_lower_limit_open',false,'c_upper_limit_open',true,...
    'figure_title','(a)','figure_name','CAMELS_I_m')
