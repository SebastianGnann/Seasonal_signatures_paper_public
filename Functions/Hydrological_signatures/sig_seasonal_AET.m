function [A,phi,ERR,...
    A_confInt,phi_confInt,ERR_confInt] = ...
    sig_seasonal_AET(Q,P,PET,t,period)
%SIG_SEASONAL_AET Calculate seasonal signatures, i.e. amplitude ratio and 
% phase shift between P-AET and Q. Uses linear equation system to fit sine 
% curve. AET is estimated using the Budyko framework.
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

% F = P-PET; % proxy for incoming signal (forcing)
w = 2*pi/period; % angular frequency 

% linear regression fitting
[A_P,phi_P,mean_P,~,~] = fitSineCurve([t P],w);
[A_PET,phi_PET,mean_PET,~,~] = fitSineCurve([t PET],w);
% get estimate for mean AET using the Budyko framework
mean_AET = BudykoCurve(mean_PET./mean_P)*mean_P; 
frac_AET = mean_AET./mean_PET; % fraction of PET becoming AET
% get AET sine wave
delta_PET = A_PET/mean(PET); % dimensionless amplitude
delta_AET = 1 - ((1-delta_PET)/frac_AET);
A_AET = delta_AET*mean_AET;
phi_AET = phi_PET;
% get forcing sine wave (P-AET)
[A_F,phi_F] = getCombinedSineWave(A_P,-A_AET,phi_P,phi_AET);
mean_F = mean_P - mean_AET;
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

%{ 
% plot AET vs PET to check results
figure; hold on;
plot(mean_P + A_P.*sin(w.*t + phi_P))
plot(-(mean_PET + A_PET.*sin(w.*t + phi_PET)))
plot(-(mean_AET + A_AET.*sin(w.*t + phi_AET)))
plot(mean_F + A_F.*sin(w.*t + phi_F));
F_comb = (mean_P + A_P.*sin(w.*t + phi_P)) + (-(mean_AET + A_AET.*sin(w.*t + phi_AET)));
plot(F_comb,'--')
xlim([0 9*365])
legend('P','PET','AET','F','F_{combined}') 
%}

function [A_comb,phi_comb] = getCombinedSineWave(A_1,A_2,phi_1,phi_2)
% helper function to get combined sine wave parameters from two sine waves
% with the same angular frequency
A_comb = sqrt(...
    (A_1.*cos(phi_1) + A_2.*cos(phi_2)).^2 ...
    +(A_1.*sin(phi_1) + A_2.*sin(phi_2)).^2);
phi_comb = atan2((A_1.*sin(phi_1) + A_2.*sin(phi_2)), ...
    (A_1.*cos(phi_1) + A_2.*cos(phi_2)));
end

% alternatively fit the sine curve to the difference (same result)
% F_comb = P-PET; % proxy for incoming signal (forcing)
% [A_F_comb,phi_F_comb,mean_F_comb,~,~] = fitSineCurve([t F_comb],w,doPlot);