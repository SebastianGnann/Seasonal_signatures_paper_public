function [A,phi,ERR,A_confInt,phi_confInt,ERR_confInt] = ...
    getSeasonalSignatures(...
    param_F,param_Q,param_F_confInt,param_Q_confInt)
%GETSEASONALSIGNATURES Calculate phase shift, amplitude ratio and effective 
% runoff ratio from fitted sinusoids.
%
% INPUT
% param_F: sine curve parameters F (amplitude, angular frequency, phase shift, offset)
% param_Q: sine curve parameters Q (amplitude, angular frequency, phase shift, offset)
% param_F_confInt: confidence intervals sine curve parameters F
% param_Q_confInt: confidence intervals sine curve parameters Q
%
% OUTPUT
% phi: vector
% A: vector with amplitude ratio (A) with phase shifts (phi=acos(A))
% ERR: effective runoff ratio (mean(Q)/mean(P-PET))
% A_confInt: confidence interval amplitude ratio
% phi_confInt: confidence interval phase shifts 
% ERR_confInt: confidence interval effective runoff ratio
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

% calculate amplitude ratio
A = param_Q(1)./param_F(1);
% calculate phase shift ... assume it's always less than 360 degrees
phase_shift = param_F(3)-param_Q(3);
% if Q "leads" P
if phase_shift < 0
    phase_shift = 2*pi+phase_shift; 
end
phi = phase_shift./param_Q(2); 
% calculate "effective runoff ratio" ERR (mean(Q)/mean(P-PET)) to check 
% whether the input and output volumes are equal
ERR = param_Q(4)./param_F(4); 

% conf. intervals (different for addition and multiplication)
A_confInt_F = param_F(1) - param_F_confInt(1);
A_confInt_Q = param_Q(1) - param_Q_confInt(1);
A_confInt = abs(A)*sqrt((A_confInt_Q/param_Q(1))^2 + (A_confInt_F/param_F(1))^2);
phi_confInt_F = param_F(3) - param_F_confInt(5);
phi_confInt_Q = param_Q(3) - param_Q_confInt(5);
phi_confInt = sqrt(phi_confInt_F^2 + phi_confInt_Q^2)./param_Q(2); 
ERR_confInt_F = param_F(4) - param_F_confInt(7);
ERR_confInt_Q = param_Q(4) - param_Q_confInt(7);
ERR_confInt = abs(ERR)*sqrt((ERR_confInt_Q/param_Q(4))^2 + (ERR_confInt_F/param_F(4))^2);

end

