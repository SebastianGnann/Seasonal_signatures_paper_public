function [BFI] = sig_BFI(B,Q)
%CALCULATEBFI Calculate BFI from streamflow and baseflow.
%
% INPUT
% B: baseflow
% Q: streamflow
%
% OUTPUT
% BFI: base flow index (= V_B/V_Q)
%
% ---
% 
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

if length(Q) ~= length(B)
    error('Records must be of the same length.')
end

% calculate BFI
BFI = nansum(B)./nansum(Q);

end

