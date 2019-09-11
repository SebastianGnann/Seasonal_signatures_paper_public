function [] = plotSerialAndParallelReservoirs(w)
%PLOTSERIALANDPARALLELRESERVOIRS Plots dampening (amplification factor) vs. 
% phase shift for two reservoirs in series and two reservoirs in parallel.
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

% outer envelope (gamma-distribution with parameter=2)
A_mat1_outer = a_range./sqrt(a_range.^2 + w.^2);
A_mat2_outer  = a_range./sqrt(a_range.^2 + w.^2);
A_theory_mat_outer  = A_mat1_outer.*A_mat2_outer;
acosA_theory_mat_outer  = acos(A_mat1_outer) + acos(A_mat2_outer);

f1 = figure('pos',[10 10 350 300]);
hold on
% grid on
p4 = fill([A_theory,fliplr(A_theory_mat_outer)],...
    [acosA_theory,fliplr(acosA_theory_mat_outer./w)],[0.8 0.8 0.8],'facealpha',0.3,'Linestyle','None');
p5 = fill([A_theory,fliplr(A_theory)],...
    [acosA_theory,fliplr(zeros(size(acosA_theory)))],[0.6 0.6 0.6],'facealpha',0.3,'Linestyle','None');
p2 = plot(A_theory,acosA_theory,'k','linewidth',1.5,'color',[0.5 0.5 0.5]);
p3 = plot(A_theory_mat_outer,acosA_theory_mat_outer./w,'--','linewidth',1.5,'color',[0.5 0.5 0.5]);

ylim([1 140])
xlim([0 1.2])
ylabel('Phase shift [days]')
xlabel('Amplitude ratio [-]')
leg = legend([p2, p3, p4, p5],...
    'Single linear reservoir','Boundary serial reservoirs',...
    'Area two serial reservoirs','Area two parallel reservoirs');
set(leg,'Position',[0.44 0.7 0.46 0.22]);
leg.ItemTokenSize = [18,10];
legend boxoff

% % update cursor
% dcm_obj = datacursormode(figure(gcf));
% set(dcm_obj,'UpdateFcn',{@myupdatefcn,ID,[1:length(ID)]})

title('(b)')

% save fig
set(f1,'Units','Inches');
position = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[position(3),position(4)]);
fig_name_raw = strcat('seasonal_signatures_serial_and_parallel_reservoir');
fig_name = regexprep(fig_name_raw,'[^a-zA-Z0-9]','');
path_name = './Seasonal_signatures_paper_public/Images';
fig_path = strcat(path_name,'\',fig_name);
print(f1,fig_path,'-dpdf','-r600');

end
