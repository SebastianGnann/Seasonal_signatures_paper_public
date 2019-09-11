function [high_Q_dur] = sig_high_Q_dur(Q)
%SIG_HIGH_Q_DUR Calculate high flow duration, defined as days with 9
% times the median daily flow.
%
% INPUT
% Q: streamflow
%
% OUTPUT
% high_Q_dur
%
% References
% Clausen, B. and Biggs, B.J.F., 2000. Flow variables for ecological 
% studies in temperate streams: groupings based on covariance. Journal of 
% hydrology, 237(3-4), pp.184-197.
%
% ---
% 
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

Q_high = 9*median(Q);
high_Q = Q>Q_high;

% find consecutive days with high flows
start1 = strfind([0,high_Q'],[0 1]);
end1 = strfind([high_Q',0],[1 0]);
interval_lengths = end1 - start1 + 1;

% calculate high flow duration
high_Q_dur = mean(interval_lengths);

end