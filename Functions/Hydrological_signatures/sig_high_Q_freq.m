function [high_Q_freq] = sig_high_Q_freq(Q)
%SIG_HIGH_Q_FREQ Calculate high flow frequency, defined as days with 9 
% times the median daily flow.
%
% INPUT
% Q: streamflow
% x: x-th percentile(s) (e.g. 0.95 for Q95)
%
% OUTPUT
% high_Q_freq
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
high_Q_num = length(Q(Q>Q_high));
high_Q_freq = high_Q_num/length(Q);

end

