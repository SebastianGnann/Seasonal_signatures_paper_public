function [Q_x] = sig_x_percentile(Q,x)
%SIG_X_PERCENTILE Calculate x-th flow percentile from streamflow.
% Following Addor et al., 2018, Q95 is a high flow measure, i.e. the 95%
% NON-exceedance probability. 
% If we have many values (>1000) the specific implementation should be
% unimportant (see e.g. Westerberg et al., 2011).
%
% INPUT
% Q: streamflow
% x: x-th percentile(s) (e.g. 95 for Q95)
%
% OUTPUT
% Q_x: x-th flow percentile (flow that is not reached for x of the time)
%
% References
% Addor, N., Nearing, G., Prieto, C., Newman, A.J., Le Vine, N. and Clark, 
% M.P., 2018. A ranking of hydrological signatures based on their 
% predictability in space. Water Resources Research, 54(11), pp.8792-8812.
% Westerberg, I.K., Guerrero, J.L., Younger, P.M., Beven, K.J., Seibert, 
% J., Halldin, S., Freer, J.E. and Xu, C.Y., 2011. Calibration of 
% hydrological models using flow-duration curves. Hydrology and Earth 
% System Sciences, 15(7), pp.2205-2227.
%
% ---
% 
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

p = 1 - x/100; % transform to get exceedance probability

% get ranks as a proxy for exceedance probabilities
Q_sorted = sort(Q);
Q_ranked = tiedrank(Q_sorted);
FDC = 1 - Q_ranked./length(Q_ranked);

% find x-th flow percentile
indices = 1:length(FDC);
bound_x = NaN(size(p));
for i = 1:length(p)
     % if flow is highly ephemeral, Q5 f.ex. is not well defined
    if isempty(max(indices(FDC >= p(i))))
    else
        bound_x(i) = max(indices(FDC >= p(i)));
    end
end
Q_x = NaN(size(p));
Q_x(~isnan(bound_x)) = Q_sorted(bound_x(~isnan(bound_x)));

end

