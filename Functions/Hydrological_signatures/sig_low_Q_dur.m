function [low_Q_dur] = sig_low_Q_dur(Q)
%SIG_LOW_Q_DUR Calculate low flow duration, defined as days with 0.2
% times the mean daily flow.
%
% INPUT
% Q: streamflow
%
% OUTPUT
% low_Q_dur
%
% References
% Olden, J.D. and Poff, N.L., 2003. Redundancy and the choice of hydrologic
% indices for characterizing streamflow regimes. River Research and 
% Applications, 19(2), pp.101-121.
%
% ---
% 
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

Q_low = 0.2*mean(Q);
low_Q = Q<Q_low;

% find consecutive days with high flows
start1 = strfind([0,low_Q'],[0 1]);
end1 = strfind([low_Q',0],[1 0]);
interval_lengths = end1 - start1 + 1;

% calculate high flow duration
low_Q_dur = mean(interval_lengths);

end