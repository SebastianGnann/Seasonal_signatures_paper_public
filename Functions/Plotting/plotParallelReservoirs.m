function [] = plotParallelReservoirs(w,panel)
%PLOTPARALLELRESERVOIRS Plots dampening (amplification factor) vs. phase
%shift for two reservoirs in parallel.
%
% INPUT
% w: angular frequency
% panel: which panel to plot
%
% OUTPUT
% plot
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

if nargin < 2
    error('Not enough input arguments.')
end

a_range = logspace(-5,0,100);
A_theory = a_range./sqrt(a_range.^2 + w.^2);
acosA_theory = acos(A_theory)./w;

a_range = logspace(-4,-1,50);
[a_mat1, a_mat2] = meshgrid([1 0.2 0.1],a_range);
A_mat1 = a_mat1./sqrt(a_mat1.^2 + w.^2);
A_mat2 = a_mat2./sqrt(a_mat2.^2 + w.^2);
phi_mat1 = acos(A_mat1);
phi_mat2 = acos(A_mat2);
if strcmp(panel,'(a)')
    f_gw = 0.3;
elseif strcmp(panel,'(b)')
    f_gw = 0.6;
elseif strcmp(panel,'(c)')
    f_gw = 0.9;
else
    error('Wrong input. Should be (a), (b), or (c).')
end
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
for i=1:length(A_theory_mat(1,:))
    p1(i) = plot(A_theory_mat(:,i),...
        acosA_theory_mat(:,i)./w,...
        '-','linewidth',1.5,'color',colour_mat(i,:));
end

ylim([1 120])
xlim([0 1.2])
ylabel('Phase shift [days]')
xlabel('Amplitude ratio [-]')
leg = legend([p1(:)', p2],'\tau_1 = 1 d','\tau_1 = 5 d','\tau_1 = 10 d','Linear reservoir');
leg.ItemTokenSize = [10,10];
legend boxoff
title(panel)

% % update cursor
% dcm_obj = datacursormode(figure(gcf));
% set(dcm_obj,'UpdateFcn',{@myupdatefcn,ID,[1:length(ID)]})

% save fig
set(f1,'Units','Inches');
position = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[position(3),position(4)]);
fig_name_raw = strcat('seasonal_signatures_parallel_reservoir',panel);
fig_name = regexprep(fig_name_raw,'[^a-zA-Z0-9]','');
path_name = './Seasonal_signatures_paper_public/Images';
fig_path = strcat(path_name,'\',fig_name);
print(f1,fig_path,'-dpdf','-r600');

end
