function [no_Q_dur] = sig_no_Q_dur(Q)
%SIG_NO_Q_DUR Calculate low flow duration, defined as days with 0 flow.
%
% INPUT
% Q: streamflow
%
% OUTPUT
% no_Q_dur
%
% ---
% 
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

no_Q = Q==0;

% find consecutive days with high flows
start1 = strfind([0,no_Q'],[0 1]);
end1 = strfind([no_Q',0],[1 0]);
interval_lengths = end1 - start1 + 1;

% calculate high flow duration
if isempty(interval_lengths)
    no_Q_dur = 0;
else
    no_Q_dur = mean(interval_lengths);
end

end