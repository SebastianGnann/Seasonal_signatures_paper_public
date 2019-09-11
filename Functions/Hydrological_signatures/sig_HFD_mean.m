function [HFD_mean] = sig_HFD_mean(Q,t)
%SIG_HFD_MEAN Calculate mean half flow date (date on which the cumulative 
% discharge since October first reaches half of the annual discharge).
%
% INPUT
% Q: streamflow
% t: date vector
%
% OUTPUT
% HFD_mean
%
% References
% Court, A., 1962. Measures of streamflow timing. Journal of Geophysical 
% Research, 67(11), pp.4335-4339.
%
% ---
% 
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

% get individual years
t_string = datetime(t,'ConvertFrom','datenum');
[years, months, days] = ymd(t_string);
year_start = min(years);
year_end = max(years);
year_list = [year_start:year_end]';

Q_temp = Q;
% Q_annual = NaN(year_end-year_start,1);
% Q_daily = NaN(365,year_end-year_start);
HFD = NaN(year_end-year_start,1);

% extract years
for y = 2:length(year_list) % since we use water years, we always start in the "2nd year"
    year = year_list(y);
    Q_water_year = ...
        [Q_temp(years==year-1 & months>=10); ...
        Q_temp(years==year & months<10)];
    Q_half_sum = 0.5*sum(Q_water_year);
    Q_cumsum = cumsum(Q_water_year);
    aux_index = 1:length(Q_water_year);
    HFD_aux = aux_index(Q_cumsum>=Q_half_sum);
    HFD(y-1) = HFD_aux(1);
end

% get mean half flow date
HFD_mean = mean(HFD);

end

