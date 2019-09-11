function [rising_limb_density] = sig_rising_limb_density(Q)
%SIG_RISING_LIMB_DENSITY Calculate rising limb density, the ratio of the 
% number of rising limbs and the total amount of time the hydrograph is 
% rising.
%
% INPUT
% Q: streamflow
%
% OUTPUT
% rising_limb_density
%
% References
% Sawicz, K., Wagener, T., Sivapalan, M., Troch, P.A. and Carrillo, G.,
% 2011. Catchment classification: empirical analysis of hydrologic 
% similarity based on catchment function in the eastern USA. Hydrology and
% Earth System Sciences, 15(9), pp.2895-2911.
%
% ---
% 
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

nr_peaks = length(findpeaks(Q)); % number of peaks in time series (higher than two surrounding values), returns 0 if no peaks
differences = diff(Q); % differences between one value and the next
up = differences(differences > 0); % sequence of positive changes (hydrograph going up)

try % try & catch structure in case tmp_peak is undefined or tmp_up is 0
    rising_limb_density = nr_peaks / length(up); % number of peaks / number of time steps that flow goes up
catch
    rising_limb_density = NaN;
end

end

