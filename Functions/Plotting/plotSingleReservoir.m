function [] = plotSingleReservoir(w)
%PLOTSINGLERESERVOIR Plots dampening (amplification factor) vs. phase
% shift for a single reservoir.
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

% get theoretical curve
a_range = logspace(-4,0,100);
A_theory = a_range./sqrt(a_range.^2 + w.^2);
acosA_theory = acos(A_theory)./w;

a_marker1 = [0.1];
A_theory_marker1 = a_marker1./sqrt(a_marker1.^2 + w.^2);
acosA_marker1 = acos(A_theory_marker1)./w;

a_marker2 = [0.01];
A_theory_marker2 = a_marker2./sqrt(a_marker2.^2 + w.^2);
acosA_marker2 = acos(A_theory_marker2)./w;

a_marker3 = [0.001];
A_theory_marker3 = a_marker3./sqrt(a_marker3.^2 + w.^2);
acosA_marker3 = acos(A_theory_marker3)./w;

f1 = figure('pos',[10 10 350 300]);
hold on
% grid on
plot(A_theory,acosA_theory,'k','linewidth',1.5)
plot(A_theory_marker1,acosA_marker1,'o k','linewidth',1.)
plot(A_theory_marker2,acosA_marker2,'d k','linewidth',1.)
plot(A_theory_marker3,acosA_marker3,'^ k','linewidth',1.)
ylim([1 120])
xlim([0 1.2])
ylabel('Phase shift [days]')
xlabel('Amplitude ratio [-]')
leg = legend('Linear Reservoir','\tau = 10 d','\tau = 100 d','\tau = 1000 d');
leg.ItemTokenSize = [10,10];
legend boxoff

% % update cursor
% dcm_obj = datacursormode(figure(gcf));
% set(dcm_obj,'UpdateFcn',{@myupdatefcn,ID,[1:length(ID)]})

% save fig
set(f1,'Units','Inches');
position = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[position(3),position(4)]);
fig_name_raw = strcat('seasonal_signatures_single_reservoir');
fig_name = regexprep(fig_name_raw,'[^a-zA-Z0-9]','');
path_name = './Seasonal_signatures_paper_public/Images';
fig_path = strcat(path_name,'\',fig_name);
print(f1,fig_path,'-dpdf','-r600');

end
