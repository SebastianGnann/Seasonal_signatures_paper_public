%% Plot results from seasonal signature analysis
%
% 	- uses the example time series to show how to create the plots shown in 
%   the paper.
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

if exist('data_Example_struc') % check if data is already loaded
else
    Example_data_struc = load('./Seasonal_signatures_paper_public/Data_and_results/Example_data_new.mat'); 
end

% store a few variables which are used very frequently 
Example_data = Example_data_struc.Example_data_new;
ID_Example = Example_data.ID;
Lat_Example = Example_data.Latitude;
Lon_Example = Example_data.Longitude;
% P_Example = Example_data.P_daily; % [Matlab_date data [mm/d]]
% PET_Example = Example_data.PET_daily; % [Matlab_date data [mm/d]]
% Q_Example = Example_data.Q_daily; % [Matlab_date data [mm/d]]

% load results
results_Example = load('seasonal_signatures_Example'); %seasonal_signatures_Example_2
seasonal_signatures_Example = results_Example.seasonal_signatures_Example;

% phase shift, amplitude ratio and offset (mean)
A_vec_Example = seasonal_signatures_Example.A_vec; % change in amplitude
phi_vec_Example = seasonal_signatures_Example.phi_vec; % phase shift
ERR_vec_Example = seasonal_signatures_Example.ERR_vec; % proxy for mean(Q)/mean(P-PET)

% confidence intervals / uncertainty (assumed negligible when using method 1)
A_vec_confInt_Example = seasonal_signatures_Example.A_vec_confInt; % change in amplitude
phi_vec_confInt_example = seasonal_signatures_Example.phi_vec_confInt; % phase shift
ERR_confInt_Example = seasonal_signatures_Example.ERR_vec_confInt; % proxy for mean(Q)/mean(P-PET)

% metadata (climate indices)
I_m_Example = seasonal_signatures_Example.I_m; % moisture index
I_mr_Example = seasonal_signatures_Example.I_mr; % moisture index seasonality
f_s_Example = seasonal_signatures_Example.f_s; % snow fraction

% annual period
period_1y = 365;
w_1y = 2*pi/period_1y;

%% plot results 

% maps (imaginary coordinates are in the UK)
plotMapSingleValueUK(Lat_Example,Lon_Example,A_vec_Example,...
    'attribute_name','Amplitude ratio','ID',ID_Example,...
    'colour_scheme','Spectral','flip_colour_scheme',false,...
    'c_limits',[0 1],'c_lower_limit_open',false,'c_upper_limit_open',true,...
    'figure_title','','figure_name','Example_amplitude_ratio')

% scatter plots (climate and catchment attributes)
plotSeasonalSignaturesData(A_vec_Example,phi_vec_Example,w_1y,...
    'attribute',I_m_Example,'attribute_name','I_m',...
    'ID',ID_Example,... 
    'A_confInt',A_vec_confInt_Example,'phi_confInt',phi_vec_confInt_example,...
    'x_limits',[0 1.2],'y_limits',[0 140],...
    'colour_scheme','Spectral','flip_colour_scheme',true,...
    'c_limits',[-1 1],'c_lower_limit_open',false,'c_upper_limit_open',true,...
    'figure_title','','figure_name','Example_I_m') 
