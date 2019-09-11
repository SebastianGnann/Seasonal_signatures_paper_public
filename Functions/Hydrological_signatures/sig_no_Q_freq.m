function [no_Q_freq] = sig_no_Q_freq(Q)
%SIG_NO_Q_FREQ Calculate no flow frequency, defined as days 0 flow.
%
% INPUT
% Q: streamflow
%
% OUTPUT
% no_Q_freq
%
% References
% Addor, N., Nearing, G., Prieto, C., Newman, A.J., Le Vine, N. and Clark, 
% M.P., 2018. A ranking of hydrological signatures based on their 
% predictability in space. Water Resources Research, 54(11), pp.8792-8812.
%
% ---
% 
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

no_Q_num = length(Q(Q==0));
no_Q_freq = no_Q_num/length(Q);

end

