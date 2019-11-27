%% Calculate seasonal signatures - example file
%
%   - the raw data can be downloaded from the sources mentioned in the
%   paper. The processed data (Matlab struc files) cannot be shared here. 
%   If you're interested in code for processing the data, please email me
%   (email shown below).
%   - uses an example file (Example_data_new.mat) to show how the code can
%   be used. The data structure file should contain the time series 
%   (P,PET,Q) in the following format: [Matlab date, variable [mm/d]] and
%   optionally some catchment attributes (e.g. frac_snow).
%
% References
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

%% load catchments
if exist('data_Example_struc') % check if data is already loaded
else
    Example_data_struc = load('./Seasonal_signatures_paper_public/Data_and_results/Example_data_new.mat'); 
end

Example_data = Example_data_struc.Example_data_new;
ID_Example = Example_data.ID;
Lat_Example = Example_data.Latitude;
Lon_Example = Example_data.Longitude;
P_Example = Example_data.P_daily; % [Matlab_date data [mm/d]]
PET_Example = Example_data.PET_daily; % [Matlab_date data [mm/d]]
Q_Example = Example_data.Q_daily; % [Matlab_date data [mm/d]]
frac_snow = Example_data.frac_snow_annual;

%% trim time series (from Oct 1989 to Sep 2009)
date_start = [1989,10,1]; % change to 1999 for example catchments
date_start = datenum(date_start);
date_end = [2009,09,30];
date_end = datenum(date_end);

P_Example = getSubPeriod(date_start,date_end,P_Example,ID_Example);
PET_Example = getSubPeriod(date_start,date_end,PET_Example,ID_Example);
Q_Example = getSubPeriod(date_start,date_end,Q_Example,ID_Example);

%% fit annual cycle to data, then calculate phase shift and amplitude ratio

% initialise arrays
% phase shift, amplitude ratio and offset (mean)
A_vec = NaN(length(ID_Example),1); % change in amplitude
A_vec2 = NaN(length(ID_Example),1); % change in amplitude 2nd method
phi_vec = NaN(length(ID_Example),1); % phase shift
phi_vec2 = NaN(length(ID_Example),1); % phase shift 2nd method
ERR_vec = NaN(length(ID_Example),1); % proxy for mean(Q)/mean(P-PET) "effective runoff ratio"
ERR_vec2 = NaN(length(ID_Example),1); % proxy for mean(Q)/mean(P-PET)

% confidence intervals / uncertainty
A_vec_confInt = NaN(length(ID_Example),1); % change in amplitude
A_vec_confInt2 = NaN(length(ID_Example),1); % change in amplitude
phi_vec_confInt = NaN(length(ID_Example),1); % phase shift
phi_vec_confInt2 = NaN(length(ID_Example),1); % phase shift
ERR_vec_confInt = NaN(length(ID_Example),1); % proxy for mean(Q)/mean(P-PET)
ERR_vec_confInt2 = NaN(length(ID_Example),1); % proxy for mean(Q)/mean(P-PET)

% sine curve parameters
param_P = NaN(length(ID_Example),3); %
param_Q = NaN(length(ID_Example),3); %

% confidence intervals (zeros for now)
param_P_confInt = zeros(length(ID_Example),6); %
param_Q_confInt = zeros(length(ID_Example),6); %

% Fourier modes
fft_max_F = NaN(length(ID_Example),1); %
fft_max_Q = NaN(length(ID_Example),1); % 

% define period
period_1y = 365;
w_1y = 2*pi/period_1y;

% load('inSubset.mat'); % subset used in modelling study

for i = 1:length(ID_Example) % nr of catchment
    
    % extract data for catchment i
    ID = ID_Example(i);
    P = P_Example{i};
    PET = PET_Example{i};
    Q = Q_Example{i};
    
    % different data checks, e.g. if record is complete, has a significant
    % snow fraction, is a benchmark record, is part of the subset
    if ~any(isnan(Q(:,2))) && ~isempty(Q) && frac_snow(i)<0.001 % && ...
        
        % calculate seasonal signatures
        % (1) linear regression
        [A_vec(i),phi_vec(i),ERR_vec(i),~,~,~] = ...
            sig_seasonal(Q(:,2),P(:,2),PET(:,2),Q(:,1));
        
        % (2) cross-covariance method
        [A_vec2(i),phi_vec2(i),ERR_vec2(i),...
            A_vec_confInt2(i),phi_vec_confInt2(i),ERR_vec_confInt2(i)] = ...
            sig_seasonal2(Q(:,2),P(:,2),PET(:,2),Q(:,1));
        
        % calculate strongest Fourier mode
        [fft_max_F(i), fft_max_Q(i)] = ...
            calcFourierSpectrum([P(:,1) P(:,2)-PET(:,2)],Q,w_1y);
        
        % plot two example catchments
        title_str = '(a) 12345 - Example catchment';
        plotSeasonalityCatchment([P(:,1) P(:,2)-PET(:,2)],Q,w_1y,ID,title_str)        
        plotFourierSpectrum([P(:,1) P(:,2)-PET(:,2)],Q,w_1y,ID,title_str)        
        
    end
    
end

%% plot results
% compare between methods
plotBivariate(A_vec,A_vec2,...
    'x_name','Amplitude ratio cross-cov.','y_name','Amplitude ratio lin. regression',...
    'ID',ID_Example,'x_limit',[0 1.5],'y_limit',[0 1.5],...
    'show_corr',true,'show_fit',false,'show_hist',false,...
    'figure_title','(a)','figure_name','scatter_amplitude_ratio_Example')
plotBivariate(phi_vec,phi_vec2,...
    'x_name','Phase shift cross-cov.','y_name','Phase shift lin. regression',...
    'ID',ID_Example,'x_limit',[0 150],'y_limit',[0 150],...
    'show_corr',true,'show_fit',false,'show_hist',false,...
    'figure_title','(b)','figure_name','scatter_phase_shift_Example')
plotBivariate(ERR_vec,ERR_vec2,...
    'x_name','Mean(Q)/Mean(P-E_p) cross-cov.','y_name','Mean(Q)/Mean(P-E_p) lin. regression',...
    'ID',ID_Example,'x_limit',[-1.2 1.2],'y_limit',[-1.2 1.2],...
    'show_corr',true,'show_fit',false,'show_hist',false,...
    'figure_title','(c)','figure_name','scatter_effective_runoff_ratio_Example')

%% save results
% method 2
seasonal_signatures_Example.ID = ID_Example;
seasonal_signatures_Example.A_vec = A_vec;
seasonal_signatures_Example.phi_vec = phi_vec;
seasonal_signatures_Example.ERR_vec = ERR_vec;
seasonal_signatures_Example.A_vec_confInt = A_vec_confInt;
seasonal_signatures_Example.phi_vec_confInt = phi_vec_confInt;
seasonal_signatures_Example.ERR_vec_confInt = ERR_vec_confInt;
seasonal_signatures_Example.I_m = -1; 
seasonal_signatures_Example.I_mr = 0;
seasonal_signatures_Example.f_s = frac_snow;
seasonal_signatures_Example.BFI = 1;
seasonal_signatures_Example.fft_max_F = fft_max_F;
seasonal_signatures_Example.fft_max_Q = fft_max_Q;
save('./Seasonal_signatures_paper_public/Data_and_results/seasonal_signatures_Example.mat','seasonal_signatures_Example')

% method 2
seasonal_signatures_Example.ID = ID_Example;
seasonal_signatures_Example.A_vec = A_vec2;
seasonal_signatures_Example.phi_vec = phi_vec2;
seasonal_signatures_Example.ERR_vec = ERR_vec2;
seasonal_signatures_Example.A_vec_confInt = A_vec_confInt2;
seasonal_signatures_Example.phi_vec_confInt = phi_vec_confInt2;
seasonal_signatures_Example.ERR_vec_confInt = ERR_vec_confInt2;
seasonal_signatures_Example.I_m = -1;
seasonal_signatures_Example.I_mr = 0;
seasonal_signatures_Example.f_s = frac_snow;
seasonal_signatures_Example.BFI = 1;
seasonal_signatures_Example.fft_max_F = fft_max_F;
seasonal_signatures_Example.fft_max_Q = fft_max_Q;
save('./Seasonal_signatures_paper_public/Data_and_results/seasonal_signatures_Example_2.mat','seasonal_signatures_Example')
