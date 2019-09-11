function [low_Q_freq] = sig_low_Q_freq(Q)
%SIG_LOW_Q_FREQ Calculate low flow frequency, defined as days with 0.2 
% times the mean daily flow.
%
% INPUT
% Q: streamflow
%
% OUTPUT
% low_Q_freq
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
low_Q_num = length(Q(Q<Q_low));
low_Q_freq = low_Q_num/length(Q);

end

