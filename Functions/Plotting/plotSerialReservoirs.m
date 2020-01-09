function [] = plotSerialReservoirs(w)
%PLOTSERIALRESERVOIRS Plots dampening (amplification factor) vs. phase
%shift for two reservoirs in series.
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

a_range = logspace(-4,0,100);
A_theory = a_range./sqrt(a_range.^2 + w.^2);
acosA_theory = acos(A_theory)./w;

% example curves
a_range = logspace(-4,0,50);
[a_mat1, a_mat2] = meshgrid([0.1 0.01 0.001],a_range);
A_mat1 = a_mat1./sqrt(a_mat1.^2 + w.^2);
A_mat2 = a_mat2./sqrt(a_mat2.^2 + w.^2);
A_theory_mat = A_mat1.*A_mat2;
acosA_theory_mat = acos(A_mat1) + acos(A_mat2);

% outer envelope (gamma-distribution with parameter=2)
A_mat1_outer = a_range./sqrt(a_range.^2 + w.^2);
A_mat2_outer  = a_range./sqrt(a_range.^2 + w.^2);
A_theory_mat_outer  = A_mat1_outer.*A_mat2_outer;
acosA_theory_mat_outer  = acos(A_mat1_outer) + acos(A_mat2_outer);

colour_mat = brewermap(10,'Set1');

f1 = figure('pos',[10 10 350 300]);
hold on
% grid on
p4 = fill([A_theory,fliplr(A_theory_mat_outer)],...
    [acosA_theory,fliplr(acosA_theory_mat_outer./w)],[0.7 0.7 0.7],'facealpha',0.3,'Linestyle','None');
p2 = plot(A_theory,acosA_theory,'k','linewidth',1.5);
p3 = plot(A_theory_mat_outer,acosA_theory_mat_outer./w,'k','linewidth',1.5,'color',[0.5 0.5 0.5]);
for i=1:length(A_theory_mat(1,:))
    p1(i) = plot(A_theory_mat(:,i),...
        acosA_theory_mat(:,i)./w,...
        '-','linewidth',1.5,'color',colour_mat(i,:));
end
annotation('textarrow',[0.38 0.28],[0.24 0.24],'String','$\tau_2$ increasing','HeadStyle','Plain',...
        'FontSize',10,'HeadWidth',6,'HeadLength',6,'Color',[0.2 0.2 0.2],'Interpreter','Latex')
ylim([1 200])
xlim([0 1.2])
ylabel('Phase shift [days]')
xlabel('Amplitude ratio [-]')
leg = legend([p1(:)', p3, p2],'\tau_1 = 10 d','\tau_1 = 100 d','\tau_1 = 1000 d','\tau_1 = \tau_2','Linear reservoir');
leg.ItemTokenSize = [10,10];
legend boxoff

% % update cursor
% dcm_obj = datacursormode(figure(gcf));
% set(dcm_obj,'UpdateFcn',{@myupdatefcn,ID,[1:length(ID)]})

% save fig
set(f1,'Units','Inches');
position = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[position(3),position(4)]);
fig_name_raw = strcat('seasonal_signatures_serial_reservoir');
fig_name = regexprep(fig_name_raw,'[^a-zA-Z0-9]','');
path_name = './Seasonal_signatures_paper_public/Images';
fig_path = strcat(path_name,'\',fig_name);
print(f1,fig_path,'-dpdf','-r600');

end
