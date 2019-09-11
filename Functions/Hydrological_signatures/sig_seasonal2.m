function [A,phi,ERR,...
    A_confInt,phi_confInt,ERR_confInt] = ...
    sig_seasonal2(Q,P,PET,t,period)
%SIG_SEASONAL2 Calculate seasonal signatures, i.e. amplitude ratio and 
% phase shift between P-PET and Q. Uses cross-covariance method to fit sine
% curve.
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
% w = 2*pi/period; % angular frequency 
n = 3*period; % parameter for cross-covariance fitting

% cross-covariance fitting method
[param_F, ~, ~, ~, ~, param_F_confInt] = TSDecomposePaperSeb([t F], n);
[param_Q, ~, ~, ~, ~, param_Q_confInt] = TSDecomposePaperSeb([t Q], n);
% param = [amplitude period phase mean]
param_F(4) = mean(F); % approximate mean of sine wave by mean of Q
param_Q(4) = mean(Q); 
% param_confInt = [amplitude_lower amplitude_upper w_lower w_upper ... 
%                   phase_lower phase_upper mean_lower mean_upper]
param_F_confInt(7:8) = mean(F); % assume zero uncertainty in mean
param_Q_confInt(7:8) = mean(Q); 

% determine phase shift, amplitude ratio and effective runoff ratio
[A,phi,ERR,...
    A_confInt,phi_confInt,ERR_confInt] = ...
    getSeasonalSignatures(...
    param_F, param_Q,...
    param_F_confInt, param_Q_confInt);

end
