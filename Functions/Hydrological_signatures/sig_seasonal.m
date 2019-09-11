function [A,phi,ERR,...
    A_confInt,phi_confInt,ERR_confInt] = ...
    sig_seasonal(Q,P,PET,t,period)
%SIG_SEASONAL Calculate seasonal signatures, i.e. amplitude ratio and phase
% shift between P-PET and Q. Uses linear equation system to fit sine curve.
%
% INPUT
% Q: streamflow
% P: precipitation
% PET: potential evapotranspiration
% t: date vector
% period: period of cycle to be extracted
%
% OUTPUT
% A: amplitude_ratio
% phi: phase_shift
% ERR: effective runoff ratio 
% A_confInt: confidence interval amplitude_ratio
% phi_confInt: confidence interval phase_shift
% ERR_confInt: confidence interval effective runoff ratio
%
% ---
% 
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

if nargin<5
    period = 365; % 1y (ignoring leap years)
end

F = P-PET; % proxy for incoming signal (forcing)
w = 2*pi/period; % angular frequency 

% linear regression fitting 
[A_F,phi_F,mean_F,~,~] = fitSineCurve([t F],w);
[A_Q,phi_Q,mean_Q,~,~] = fitSineCurve([t Q],w);
% param = [amplitude period phase mean]
param_F = [A_F w phi_F mean_F]; 
param_Q = [A_Q w phi_Q mean_Q];
% TODO: add confidence intervals *
% param_confInt = [amplitude_lower amplitude_upper w_lower w_upper ... 
%                   phase_lower phase_upper mean_lower mean_upper]
param_F_confInt = [param_F(1) param_F(1) param_F(2) param_F(2) ...
    param_F(3) param_F(3) param_F(4) param_F(4)]; 
param_Q_confInt = [param_Q(1) param_Q(1) param_Q(2) param_Q(2) ...
    param_Q(3) param_Q(3) param_Q(4) param_Q(4)]; 

% determine phase shift, amplitude ratio and effective runoff ratio
[A,phi,ERR,...
    A_confInt,phi_confInt,ERR_confInt] = ...
    getSeasonalSignatures(...
    param_F, param_Q,...
    param_F_confInt, param_Q_confInt);

end

% * Parameter uncertainty is negligible if the other method 
% (cross-covariance)is used, therefore we assume it to be zero for now.
% Note that we are not interested in how well the sine wave fits the time  
% series (e.g. the peaks), but in the uncertainty of the seasonal harmonic. 
