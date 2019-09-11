function [streamflow_elasticity] = sig_streamflow_elasticity(Q,P,t)
%SIG_STREAMFLOW_ELASTICITY Calculate streamflow elasticity, the sensitivity
% of a catchment’s streamflow response to changes in precipitation at the 
% annual time scale
%
% INPUT
% Q: streamflow
% P: precipitation
% t: date vector
%
% OUTPUT
% streamflow_elasticity
%
% References
% Sawicz, K., Wagener, T., Sivapalan, M., Troch, P.A. and Carrillo, G.,
% 2011. Catchment classification: empirical analysis of hydrologic 
% similarity based on catchment function in the eastern USA. Hydrology and
% Earth System Sciences, 15(9), pp.2895-2911.
%
% ---
% 
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

% get annual values
t_string = datetime(t,'ConvertFrom','datenum');
[years, months, days] = ymd(t_string);
year_start = min(years);
year_end = max(years);
year_list = [year_start:year_end]';

Q_temp = Q;
P_temp = P;
Q_annual = NaN(year_end-year_start,1);
% Q_daily = NaN(365,year_end-year_start); 
P_annual = NaN(year_end-year_start,1);
% P_daily = NaN(365,year_end-year_start); 

% extract years
for y = 2:length(year_list) % since we use water years, we always start in the "2nd year"
    year = year_list(y);
    Q_water_year = ...
        [Q_temp(years==year-1 & months>=10); ...
        Q_temp(years==year & months<10)];
    Q_annual(y-1) = mean(Q_water_year);
%     Q_daily(:,y-1) = Q_water_year(1:365); % skip last day of each leap year
    P_water_year = ...
        [P_temp(years==year-1 & months>=10); ...
        P_temp(years==year & months<10)];
    P_annual(y-1) = mean(P_water_year);
%     P_daily(:,y-1) = P_water_year(1:365); % skip last day of each leap year
end

% calculate elasticity = dQ/dP * P/Q
dQ = Q_annual - mean(Q);
dP = P_annual - mean(P);
streamflow_elasticity = median((dQ./dP)*(mean(P)/mean(Q)));

end

