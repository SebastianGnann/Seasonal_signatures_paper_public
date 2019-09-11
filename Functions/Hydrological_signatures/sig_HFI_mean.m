function [HFI_mean] = sig_HFI_mean(Q,t)
%SIG_HFI_MEAN Calculate mean half flow interval (time span between the 
% date on which the cumulative discharge since October first reaches a 
% quarter of the annual discharge and the date on which the cumulative 
% discharge since October first reaches three quarters of the annual 
% discharge)
%
% INPUT
% Q: streamflow
% t: date vector
%
% OUTPUT
% HFI_mean
%
% References
% Court, A., 1962. Measures of streamflow timing. Journal of Geophysical 
% Research, 67(11), pp.4335-4339.
% get individual years
%
% ---
% 
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

t_string = datetime(t,'ConvertFrom','datenum');
[years, months, days] = ymd(t_string);
year_start = min(years);
year_end = max(years);
year_list = [year_start:year_end]';

Q_temp = Q;
% Q_annual = NaN(year_end-year_start,1);
% Q_daily = NaN(365,year_end-year_start);
HFI = NaN(year_end-year_start,1);

% extract years
for y = 2:length(year_list) % since we use water years, we always start in the "2nd year"
    year = year_list(y);
    Q_water_year = ...
        [Q_temp(years==year-1 & months>=10); ...
        Q_temp(years==year & months<10)];
    Q_25_sum = 0.25*sum(Q_water_year);
    Q_75_sum = 0.75*sum(Q_water_year);
    Q_cumsum = cumsum(Q_water_year);
    aux_index = 1:length(Q_water_year);
    HFI_aux_25 = aux_index(Q_cumsum>Q_25_sum);
    HFI_aux_75 = aux_index(Q_cumsum>Q_75_sum);
    if isempty(HFI_aux_25) || isempty(HFI_aux_75) % if there is no flow
        %disp('')
    else
        HFI_25 = HFI_aux_25(1);
        HFI_75 = HFI_aux_75(1);
        HFI(y-1) = HFI_75 - HFI_25;
    end
end

% get mean half flow date (ignoring NaNs)
HFI_mean = nanmean(HFI);

end

