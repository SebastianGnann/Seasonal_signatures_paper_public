function [] = plotSeasonalityCatchment(F,Q,w,ID,title_str)
%PLOTSEASONALITYCATCHMENT Plots forcing (e.g. P-PET) and streamflow and 
% their seasonal components.
%
% INPUT
% F: forcing, typically approximated by P - PET [mm]
% Q: streamflow [mm]
% w: angular frequency [1/days]
% ID: catchment ID
% title: catchment title
%
% OUTPUT
% plot
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

if nargin<5
    title_str = num2str(ID);
end

if nargin < 4
    error('Not enough input arguments.')
end

% calculate sine curve
[A_F,phi_F,offset_F,~,~] = fitSineCurve(F,w);
[A_Q,phi_Q,offset_Q,~,~] = fitSineCurve(Q,w);
F_sin = offset_F + A_F*sin(w.*F(:,1) + phi_F);
Q_sin = offset_Q + A_Q*sin(w.*Q(:,1) + phi_Q);

% note that the phase is with respect to a different reference datum -
% since we're only interested in the phase difference this does not matter
% except for plotting
% n = 3*365;
% [param_F, ~, ~, ~, ~, param_F_confInt] = TSDecomposePaperSeb(F, n);
% [param_Q, ~, ~, ~, ~, param_Q_confInt] = TSDecomposePaperSeb(Q, n);
% F_sin2 = mean(F(:,2)) + param_F(1)*sin(w.*(F(:,1)-F(1,1)) + param_F(3));
% Q_sin2 = mean(Q(:,2)) + param_Q(1)*sin(w.*(Q(:,1)-Q(1,1)) + param_Q(3));



n = 10;
colour_mat = (brewermap(n,'Set1')); %flip
fig1 = figure('Name',title_str,'NumberTitle','off','pos',[10 10 300 250]); % 

subplot(3,1,1)
title(title_str)
hold on
plot(F(:,1),movmean(F(:,2),30),'Linewidth',1,'color',colour_mat(9,:))
plot(F(:,1),F_sin,'Linewidth',1.5,'color',colour_mat(2,:))
%             xlabel('Time [d]')
ylabel('P - E_p [mm/d]')
tick5 = datenum(2000:2:2009,10,1);
tick5 = tick5 + 92;
set(gca, 'xtick', tick5);
datetick('x', 'yyyy', 'keepticks');
xlim([median(F(:,1)) max(F(:,1))])
% ylim([-10 round(max(P_eff(:,2)))]) %round(max(P(:,2)))
% ylim([-4 8]) %ylim([-5 40]) % ylim([-2 10]) %
%             xlim([P(round(end/2),1) P(end,1)])
%             title('(a)')
subplot(3,1,2)
hold on
plot(Q(:,1),movmean(Q(:,2),30),'Linewidth',1,'color',colour_mat(9,:))
plot(Q(:,1),Q_sin,'Linewidth',1.5,'color',colour_mat(5,:))
%             xlabel('Time [d]')
ylabel('Q [mm/d]')
tick5 = datenum(2000:2:2009,10,1);
tick5 = tick5 + 92;
set(gca, 'xtick', tick5);
datetick('x', 'yyyy', 'keepticks');
xlim([median(F(:,1)) max(F(:,1))])
% ylim([0 round(max(Q(:,2)))])
% ylim([-4 8]) %ylim([-5 40]) % ylim([-2 10]) %
%             xlim([P(round(end/2),1) P(end,1)])
%             title('(b)')
subplot(3,1,3)
hold on
plot(F(:,1),F_sin-mean(F_sin),'Linewidth',1.5,'color',colour_mat(2,:))
plot(Q(:,1),Q_sin-mean(Q_sin),'Linewidth',1.5,'color',colour_mat(5,:))
xlabel('Year')
ylabel('Seasonal [mm/d]')
tick5 = datenum(2000:2:2009,10,1);
tick5 = tick5 + 92;
set(gca, 'xtick', tick5);
datetick('x', 'yyyy', 'keepticks');
xlim([median(F(:,1)) max(F(:,1))])
%             xlim([P(round(end/2),1) P(end,1)])

% save fig
% NOTE RESOLUTION SEEMS TO BE WORSE WHEN USING THIS SHAPEFILE
set(fig1,'Units','Inches');
pos = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos(3),pos(4)]);
fig_name_raw = strcat('seasonal_components_','_',num2str(ID));
fig_name = regexprep(fig_name_raw,'[^a-zA-Z0-9]','');
path_name = './Seasonal_signatures_paper_public/Images';
fig_path = strcat(path_name,'\',fig_name);
print(fig1,fig_path,'-dpdf','-r500');
end

