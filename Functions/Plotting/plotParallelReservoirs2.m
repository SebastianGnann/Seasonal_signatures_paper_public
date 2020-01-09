function [] = plotParallelReservoirs2(w)
%PLOTPARALLELRESERVOIRS2 Plots dampening (amplification factor) vs. phase
%shift for two reservoirs in parallel.
%
% INPUT
% w: angular frequency
%
% OUTPUT
% plot
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

if nargin < 1
    error('Not enough input arguments.')
end

a_range = logspace(-5,0,100);
A_theory = a_range./sqrt(a_range.^2 + w.^2);
acosA_theory = acos(A_theory)./w;

[a_mat1, a_mat2] = meshgrid([1],[0.1 0.01 0.001]);
A_mat1 = a_mat1./sqrt(a_mat1.^2 + w.^2);
A_mat2 = a_mat2./sqrt(a_mat2.^2 + w.^2);
phi_mat1 = acos(A_mat1);
phi_mat2 = acos(A_mat2);
f_gw = [0:0.01:1]; 
A_theory_mat = sqrt(...
    ((1-f_gw).*A_mat1.*cos(phi_mat1)+f_gw.*A_mat2.*cos(phi_mat2)).^2 ...
   +((1-f_gw).*A_mat1.*sin(phi_mat1)+f_gw.*A_mat2.*sin(phi_mat2)).^2);
acosA_theory_mat = atan(...
    ((1-f_gw).*A_mat1.*sin(phi_mat1)+f_gw.*A_mat2.*sin(phi_mat2))./ ...
    ((1-f_gw).*A_mat1.*cos(phi_mat1)+f_gw.*A_mat2.*cos(phi_mat2)));

colour_mat = brewermap(10,'Set1');

f1 = figure('pos',[10 10 350 300]);
hold on
% grid on
p4 = fill([A_theory,fliplr(A_theory)],...
    [acosA_theory,fliplr(zeros(size(acosA_theory)))],[0.7 0.7 0.7],'facealpha',0.3,'Linestyle','None');
p2 = plot(A_theory,acosA_theory,'k','linewidth',1.5);
p3 = plot(A_theory,zeros(size(acosA_theory)),'k','linewidth',1.5,'color',[0.5 0.5 0.5]);
for i=1:length(A_theory_mat(:,1))
    p1(i) = plot(A_theory_mat(i,:),...
        acosA_theory_mat(i,:)./w,...
        '-','linewidth',1.5,'color',colour_mat(i,:));
end
annotation(f1,'textbox',[0.20 0.80 0.25 0.09],'String',['$0 \le p \le 1$'],'LineStyle','None','Interpreter','Latex')
annotation('textarrow',[0.4 0.33],[0.27 0.43],'String','$p$ increasing','HeadStyle','Plain',...
        'FontSize',10,'HeadWidth',6,'HeadLength',6,'Color',[0.2 0.2 0.2],'Interpreter','Latex')
ylim([1 120])
xlim([0 1.2])
ylabel('Phase shift [days]')
xlabel('Amplitude ratio [-]')
leg = legend([p1(:)', p2],'\tau_2 = 10 d','\tau_2 = 100 d','\tau_2 = 1000 d','Linear reservoir');
leg.ItemTokenSize = [10,10];
legend boxoff
title('(d)')

% % update cursor
% dcm_obj = datacursormode(figure(gcf));
% set(dcm_obj,'UpdateFcn',{@myupdatefcn,ID,[1:length(ID)]})

% save fig
set(f1,'Units','Inches');
position = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[position(3),position(4)]);
fig_name_raw = strcat('seasonal_signatures_parallel_reservoir_d');
fig_name = regexprep(fig_name_raw,'[^a-zA-Z0-9]','');
path_name = './Seasonal_signatures_paper_public/Images';
fig_path = strcat(path_name,'\',fig_name);
print(f1,fig_path,'-dpdf','-r600');

end
