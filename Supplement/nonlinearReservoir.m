%% Numerical solution for non-linear reservoir with sinusoidal input
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

clc
% close all
% clear all

%% add directories for functions to path

if exist('./BrewerMap') == 7
    addpath(genpath('./BrewerMap'));
else
    error('BrewerMap toolbox needed. Can be download from https://github.com/DrosteEffect/BrewerMap and should be in a folder named BrewerMap in the same directory.')
end

%% specify input, parameters, etc.

% time and periodicity
period = 365; % 1y
t_end = period*300;
dt = 1;
t = [0:dt:t_end]'; % time vector
len_t = length(t);
w = 2*pi/period; % angular frequency

% rainfall (or recharge) input
P_0 = 0.5; % amplitude
P_sin = P_0.*(1 + sin(w.*t)); % periodic rainfall, no phase shift
% rng(66); % random seed
% P_prime = abs(randn(len_t,1)); % random rainfall 
P = P_sin;% + P_prime; % total rainfall
% ET = zeros(len_t,1); % evapotranspiration

% reservoir parameters
%{
reservoir eq. (Kirchner, 2009): Q = Q_ref*((S-S_ref)/m)^n
m = Qref^(2-beta)/((2-beta)*alpha)
n = 1/(2-beta)
%} 
S0 = 1.0; % initial storage
S_ref = 0; % reference storage
Q_ref = 1; % reference discharge
alpha_vec = flip(logspace(0,-4,25))'; 
beta_vec = [0.1:0.1:1.9]'; 
n_vec = 1./(2-beta_vec); 
% a_vec = ((2-beta_vec).*alpha_vec).^(1./(2-beta_vec)); 

% preallocate vectors to determine phase shift
n_alpha = length(alpha_vec);
n_beta = length(n_vec);
maxP = NaN(n_alpha,n_beta);
maxSnon = NaN(n_alpha,n_beta);
minP = NaN(n_alpha,n_beta);
minSnon = NaN(n_alpha,n_beta);
diff_non = NaN(n_alpha,n_beta);
diff_max_non = NaN(n_alpha,n_beta);
diff_min_non = NaN(n_alpha,n_beta);
amplitude_change = NaN(n_alpha,n_beta);
a_mat = NaN(n_alpha,n_beta);
n_mat = NaN(n_alpha,n_beta);
t_rec_mat = NaN(n_alpha,n_beta);
phase_lag_days = NaN(n_alpha,1);

%% get solutions

for j = 1:n_alpha
    
    for k = 1:n_beta
        
        a = ((2-beta_vec(k)).*alpha_vec(j)).^(1./(2-beta_vec(k))); % 2.5*10^-5;
        m = Q_ref.^(2-beta_vec(k))./((2-beta_vec(k)).*alpha_vec(j)); % 200;
        n = n_vec(k); % 2;
        
%         a_mat(j,k) = alpha_vec(j);
%         n_mat(j,k) = beta_vec(k);
        
        % non-linear reservoir
        S_non = NaN(len_t,1);
        S_non(1) = S0;
        for i = 2:length(t)
            S_non(i) = S_non(i-1) - dt*Q_ref*((S_non(i-1)-S_ref)/m)^n + (P(i)+P(i-1))/2*dt; %a*dt*S_non(i-1)^n
            % water balance control
            if S_non(i)<0 
                S_non(i) = 0;
            end
        end
        
        % determine phase shift
        [maxPValue, maxP(j,k)] = max(P(end-366:end));
        [minPValue, minP(j,k)] = min(P(end-366:end));
        [maxSnonValue, maxSnon(j,k)] = max(S_non(end-366:end));
        [minSnonValue, minSnon(j,k)] = min(S_non(end-366:end));
        if maxSnon(j,k) < maxP(j,k)
            diff_max_non(j,k) = - maxSnon(j,k) + maxP(j,k);
            diff_min_non(j,k) = - minSnon(j,k) + minP(j,k);
            diff_non(j,k) = (diff_max_non(j,k) + diff_min_non(j,k))/2;
        else
            diff_max_non(j,k) = maxSnon(j,k) - maxP(j,k);
            diff_min_non(j,k) = minSnon(j,k) - minP(j,k);
            diff_non(j,k) = (diff_max_non(j,k) + diff_min_non(j,k))/2;
        end
        amplitude_change(j,k) = (Q_ref.*((maxSnonValue-S_ref)./m).^n - Q_ref.*((minSnonValue-S_ref)./m).^n)./(maxPValue-minPValue);
        
        % define recession time: Q50 to Q90
        Q90 = quantile(Q_ref.*((S_non(end-366:end)-S_ref)./m).^n,0.1);
        Q50 = quantile(Q_ref.*((S_non(end-366:end)-S_ref)./m).^n,0.5);
        
        if beta_vec(k) == 1
            t_rec_mat(j,k) = -(log(Q90/Q50))/alpha_vec(j);
        else
            t_rec_mat(j,k) = (-Q90^(1-beta_vec(k)) + Q50^(1-beta_vec(k)))./...
                (alpha_vec(j)*(1-beta_vec(k))); % Stoelzle et al., 2013
        end
                
        %% plot (example reservoir)
%         colour_mat = (brewermap(10,'Set1')); %flip
%         f1 = figure('pos',[10 10 350 150]);
%         hold off
%         plot(t-t(end-3*365),P,'Linewidth',1.5,'color',colour_mat(2,:))
%         hold on
%         ylabel('Precipitation')
%         plot(t-t(end-3*365),Q_ref.*((S_non-S_ref)./m).^n,'-','Linewidth',1.5,'color',colour_mat(5,:))
%         xlabel('Time [days]'); ylabel('Flow [mm]')
%         leg = legend('Q_{in}','Q_{out}','location','ne');
%         set(leg,'Position',[0.79 0.68 0.13 0.27]);
%         leg.ItemTokenSize = [10,10];
%         legend boxoff        
%         xlim([t(end-3*365)-t(end-3*365) t(end)-t(end-3*365)])
%         ylim([0 1])
%         box off
        
        % save fig
%         set(f1,'Units','Inches');
%         position = get(f1,'Position');
%         set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[position(3),position(4)]);
%         fig_name_raw = strcat('nonlinear_numerical_reservoir');
%         fig_name = regexprep(fig_name_raw,'[^a-zA-Z0-9]','');
%         path_name = './Seasonal_signatures_paper_public/Images';
%         fig_path = strcat(path_name,'\',fig_name);
%         print(f1,fig_path,'-dpdf','-r600');
        
        %% display (example reservoir)
%         disp(mean(Q_ref.*((S_non(end-366:end)-S_ref)./m).^n))
%         disp(mean(Q_ref.*((maxSnonValue-S_ref)./m).^n) - mean(mean(Q_ref.*((S_non(end-366:end)-S_ref)./m).^n)))
%         disp(mean(Q_ref.*((S_non(end-366:end)-S_ref)./m).^n) - Q_ref.*((minSnonValue-S_ref)./m).^n)
%         disp(diff_max_non(j,k))
%         disp(diff_min_non(j,k))
    end
    
    A = alpha_vec(j)/sqrt(alpha_vec(j)^2 + w^2);
    phase_lag_days(j) = acos(A)/w; % analytical phase shift linear res.
    
end


%% plot phase shift vs dampening
colour_mat = (brewermap(9,'Set1')); %flip
f1 = figure('pos',[10 10 350 300]);
a_range = logspace(-6,-0,100);
A_theory = a_range./sqrt(a_range.^2 + w.^2);
acosA_theory = acos(A_theory)/2/pi*360;
hold on
phaseshift = reshape(diff_non,n_alpha*n_beta,1);
phaseshift_max = reshape(diff_max_non,n_alpha*n_beta,1);
phaseshift_min = reshape(diff_min_non,n_alpha*n_beta,1);
dampening = reshape(amplitude_change,n_alpha*n_beta,1);
scatter(dampening,phaseshift_max./365.*360,25,'filled','MarkerFaceColor',[0.3 0.3 0.3],'MarkerFaceAlpha',1.0)
scatter(dampening,phaseshift_min./365.*360,25,'filled','MarkerFaceColor',[0.7 0.7 0.7],'MarkerFaceAlpha',1.0)
scatter(dampening,phaseshift./365.*360,25,'filled','MarkerFaceColor',colour_mat(5,:))
plot(A_theory,acosA_theory,'k','linewidth',1.5)
% set(gca,'yscale','log')
ylabel('Phase shift [days]')
xlabel('Amplitude ratio [-]')
xlim([0 1.2])
ylim([0 120])
leg = legend('Phase shift maxima','Phase shift minima','Mean phase shift','Linear reservoir');
%grid on
leg.ItemTokenSize = [10,10];
legend boxoff

% save fig
set(f1,'Units','Inches');
position = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[position(3),position(4)]);
fig_name_raw = strcat('seasonal_signatures_nonlinear_reservoir');
fig_name = regexprep(fig_name_raw,'[^a-zA-Z0-9]','');
path_name = './Seasonal_signatures_paper_public/Images';
fig_path = strcat(path_name,'\',fig_name);
print(f1,fig_path,'-dpdf','-r600');

% diff_max_non = reshape(diff_max_non,n_alpha*n_beta,1);
% a_mat = reshape(a_mat,n_alpha*n_beta,1);
% n_mat = reshape(n_mat,n_alpha*n_beta,1);

%% plot peak lag
% figure
% grid on
% hold on
% plot(a_mat,ones(size(n_mat)).*period/4,'k --')
% for i = 1:length(a_mat(1,:))
%     %     colour_mat = (brewermap(length(a_mat(1,:)),'Spectral')); %RdGy PiYG
%     %     plot(a_mat(:,i),diff_max_non(:,i),'o','linewidth',2,'color',colour_mat(i,:))
%     scatter(a_mat(:,i),diff_non(:,i),50,n_mat(:,i),'filled') %1./t_rec_mat(:,i)
%     c = colorbar;
%     title(c,'n [-]')
% end
% plot(a_mat(:,1),phase_lag_days,'k -','linewidth',2)
% xlabel('a [f(n)]')
% ylabel('Phase lag [d]')
% title('Q = aS^n')
% ylim([0 200])
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% % legend('max. phase shift','peak difference',...
% %     'analytical phase lag','peak difference non-linear','location','best')
% % set(gca,'Xdir','reverse')
