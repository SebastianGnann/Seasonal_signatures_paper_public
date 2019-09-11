%% Create example time-series using single linear reservoir
%
%   - this time series is used to show how the code can be used. It's based
%   on a single linear reservoir forced with artificial rainfall.
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

% clc
% close all
% clear all

%% single linear reservoir with synthetic recharge
% define parameters
t_start = datenum(1979,10,1); % 10y warm up period
t_end = datenum(2009,09,30);
dt = 1;
t = [t_start:dt:t_end]'; % time
len_t = length(t);
a_vec = flip(logspace(-2,-1,5)); % decay constant % 1/tau
na = length(a_vec);
S0 = 50.0; % initial storage
% P = ones(length(t),1); % precip
w = 2*pi/365;
P_0 = 1; % amplitude
P_sin = P_0.*(1 + sin(w.*(t-t_start))); % periodic rainfall
rng(33); % random seed
P_prime = 10.*abs(randn(len_t,1)); % random rainfall
PET = zeros(len_t,1); % pot. evapotranspiration - zero for convenience
P = P_sin + P_prime - PET; % total "effective" rainfall

% initialise arrays
phase_shift = NaN(na,1);
phase_shift_days = NaN(na,1);

% solve for storage and outflow
for j = 1:na
    
    a = a_vec(j);
    
    % rainfall response
    S_prime = NaN(len_t,1);
    for k = 1:len_t
        S_prime(k) = sum(P_prime([1:k]).*exp(-a.*(t([k:-1:1])-t_start))).*dt;
    end
    
    % exponential component (from initial storage)
    S_ana_1 = S0*exp(-a*(t-t_start));
    % exponential component from sinusoidal offset 
    S_ana_2 = P_0/a*(1-exp(-a*(t-t_start))) ;
    % sinusoidal component from sinusoidal input
    prefac = P_0/(w^2+a^2);
    S_ana_3 = prefac.*(w.*exp(-a.*(t-t_start)) ...
        + sqrt(w^2+a^2)*(sin(w.*(t-t_start) + atan(-w/a))));
    % random rainfall input
    S_ana_4 = S_prime;
    
    % sum up
    S_ana = S_ana_1 + S_ana_2 + S_ana_3 + S_ana_4;
    
    % numerical
    S_num = NaN(len_t,1);
    S_num(1) = S0 + (P(1)+P(2))/2*dt; % pick start value
    
    for i = 2:length(t)
        % use explicit scheme for now
        S_num(i) = S_num(i-1) - a*dt*S_num(i-1) + (P(i)+P(i-1))/2*dt;
    end
    
    % remove warm up period
    Q_Example = a.*S_ana(end-20*365-4:end); % account for leap years
    Q_num_Example = a.*S_num(end-20*365-4:end); % account for leap years
    P_Example = P(end-20*365-4:end);
    PET_Example = PET(end-20*365-4:end); 
    t_Example = t(end-20*365-4:end);
    
    % determine phase lag using seasonal signatures 
    [~,phase_shift(j),~,~,~,~] = ...
            sig_seasonal(Q_Example,P_Example,PET_Example,t_Example);
        
    phase_shift_analytical = atan(-w/a);
    phase_shift_days(j) = -phase_shift_analytical/w;
    
end

%% plot results
% plot peak lag
figure
grid on
hold on
plot(1./a_vec,phase_shift,'linewidth',2)
plot(1./a_vec,phase_shift_days,'--','linewidth',2)
plot(1./a_vec,1./a_vec,'k --','linewidth',1)
xlabel('\tau [d]')
ylabel('Peak Time [d]')
legend('phase shift',...
    'analytical phase shift','1:1 line','location','nw')
% set(gca,'Xdir','reverse')

% plot last time series as an example
figure
yyaxis right
bar(t_Example,P_Example,0.5)
set(gca,'Ydir','reverse')
axis([t_Example(1) t_Example(end) 0  3*(max(P_Example))])
datetick
ylabel('Precipitation [mm]')
hold on
grid on
yyaxis left
plot(t_Example,Q_Example,'b -','linewidth',2)
plot(t_Example,Q_num_Example,'r --','linewidth',2) % compare to numerical solution
axis([t_Example(1) t_Example(end) 0.5*min(Q_Example) 1.5*max(Q_Example)])
datetick
xlabel('Date'); ylabel('Flow [mm]')
legend('Analytical Solution','Numerical Solution','Precipitation','location','se')

%% save results
Example_data_new.ID = 12345;
Example_data_new.Latitude = 51.456;
Example_data_new.Longitude = -2.602;
Example_data_new.Q_daily = {[t_Example, Q_Example]};
Example_data_new.P_daily = {[t_Example, P_Example]};
Example_data_new.PET_daily = {[t_Example, PET_Example]};
Example_data_new.Q_b_daily = {[t_Example, NaN(size(t_Example))]};
Example_data_new.T_daily = {[t_Example, NaN(size(t_Example))]};
Example_data_new.frac_snow_annual = 0.0;
save('./Seasonal_signatures_paper_public/Data_and_results/Example_data_new.mat','Example_data_new');
