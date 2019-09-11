function [BFI_UKIH] = calc_BFI_UKIH(Q_tot)
%CALC_BFI_LYNEHOLLICK: Calculate baseflow and BFI using UKIH "smoothed
% minima" method (UK Institute of Hydrology, 1980)
%
% INPUT
% Q_tot: total streamflow 
%
% OUTPUT
% BFI calculated using UKIH method
%
% References
% Institute of Hydrology (Great Britain), 1980. Low Flow Studies Reports. 
% Institute of Hydrology.
%
% ---
% 
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

% remove NaN values
if any(isnan(Q_tot))
    disp('NaN values are assigned the median value for baseflow calculation')
end
Q_tot(isnan(Q_tot)) = nanmedian(Q_tot);

% calculate baseflow
[B_UKIH, t_UKIH] = UKIH_Method(Q_tot, 5); % 5 is the default parameter

% use minimum baseflow to fill in missing values at the beginning
B_temp = min(Q_tot)*ones(size(Q_tot));
B_temp(t_UKIH) = B_UKIH;
B_UKIH = B_temp;

% calculate BFI
BFI_UKIH = sig_BFI(B_UKIH,Q_tot);

end