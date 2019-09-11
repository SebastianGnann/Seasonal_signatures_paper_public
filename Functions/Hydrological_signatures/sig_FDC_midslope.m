function [mid_slope] = sig_FDC_midslope(Q)
%SIG_FDC_MIDSLOPE Calculate mid slope of flow duration curve (FDC) from 
% streamflow (e.g. Sawicz et al., 2011).
%
% INPUT
% Q: streamflow
%
% OUTPUT
% mid_slope: slope of the flow duration curve (between the log-transformed 
% 33rd and 66th streamflow percentiles)
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

% get ranks as a proxy for exceedance probabilities
Q_sorted = sort(Q);
Q_ranked = tiedrank(Q_sorted); 
FDC = 1 - Q_ranked./length(Q_ranked);

% slow of normalised FDC between 0.33 and 0.66 percentile
indices = 1:length(FDC);
bound_66 = max(indices(FDC >= 0.66));
bound_33 = max(indices(FDC >= 0.33));

% mid_slope = (Q_sorted(bound_67) - Q_sorted(bound_33))./mean(Q_sorted)./...
%     (FDC(bound_67) - FDC(bound_33)) ;

% using log Q values
mid_slope = (log(Q_sorted(bound_66)) - log(Q_sorted(bound_33)))./...
     (FDC(bound_66) - FDC(bound_33));

% optional plotting
% x = FDC(bound_66:bound_33);
% c = log(Q_sorted(bound_33)) - FDC(bound_33)*mid_slope;
% y = mid_slope.*x + c;
% % 
% figure
% % plot(FDC,Q_sorted./mean(Q_sorted))
% plot(FDC,(Q_sorted))
% hold on
% plot(x,exp(y))
% legend('FDC','Mid slope')

% if flow is very intermittent (e.g. 66th percentile is 0)
if isempty(mid_slope)
    mid_slope = NaN;
end

% calculate difference between fitted mid FDC and actual FDC
% residuals = mean(y - Q_sorted(bound_66:bound_33)./mean(Q_sorted));

end

