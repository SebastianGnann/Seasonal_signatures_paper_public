function [B, t] = UKIH_Method(Q,n_days)
%UKIH_METHOD UKIH "smoothed minima" method (UK Institute of Hydrology, 
% 1980). Returns a time series that is shorter than the original time 
% series since a few days at the beginning are omitted. Needs manual
% adjustment.
%
% INPUT
% Q: streamflow
% n_days: length of data blocks
%
% OUTPUT
% B: baseflow
% t: time (start of record = 1)
%
% References
% Institute of Hydrology (Great Britain), 1980. Low Flow Studies Reports. 
% Institute of Hydrology.
%
% ---
% 
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

if nargin < 2
    n_days = 5; % default is 5 
end

n = length(Q);
% B = NaN(n,1);
% B(1) = 0; % initial condition

Q_min5 = NaN(round(n/n_days),1); % 5 day minima
min_i = NaN(round(n/n_days),1); % corresponding indices
ind5 = 1; % minima counter
TP = 0; % turning points
t_TP = 0; % corresponding time/index
indTP = 1; % TP counter

for i = 1:n_days:floor(n/n_days)*n_days % divide in non-overlapping n-day blocks
    [Q_min5(ind5), min_i(ind5)] = min(Q(i:i+(n_days-1))); % find minimum
    if ind5 <= 2 % need at least three minima
    elseif      Q_min5(ind5-1)*0.9 < Q_min5(ind5-2) ...
            &&  Q_min5(ind5-1)*0.9 < Q_min5(ind5) % check if baseflow ordinate
        TP(indTP) = Q_min5(ind5-1);
        t_TP(indTP) = i - n_days - 1 + min_i(ind5-1); % get corresponding index
        indTP = indTP + 1;
    end
    ind5 = ind5 + 1;
end

t = [t_TP(1):t_TP(end)]';
B = interp1q(t_TP',TP',t); % linear interpolation
Qt = Q(t);
B(B>Qt) = Qt(B>Qt); % constrain B, so that B is never larger than Q

end

