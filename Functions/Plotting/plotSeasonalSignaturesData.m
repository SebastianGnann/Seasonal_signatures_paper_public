function [] = plotSeasonalSignaturesData(A,phi,w,varargin)
%PLOTSEASONALSIGNATURESDATA Plots scatter plot of amplitude ratio and phase 
% shift and plots lines for single reservoir and two reservoirs in series.
% Options:
%   - colour dots according to third variable
%   - plot confidence intervals
%   - various plotting options, e.g. axes limits
%   (- save plot as PDF)
%
% INPUT
% A: amplitude ratio
% phi: phase shift
% w: angular frequency, defines the lines for the reservoirs
% attribute: attribute to be coloured in, e.g. BFI
% attribute_name: name of attribute
% ID: catchment ID
% A_confInt: confidence interval 1y
% phi_confInt: confidence interval 1y
% x_limit: axis limit for x-axis (A)
% y_limit: axis limit for y-axis (phi)
% colour_scheme: name of colour scheme
% flip_colour_scheme: flip colour scheme?
% c_limits: limits of colour axis, e.g. [0 1]
% c_lower_limit_open: is the lower limit open? 
% c_upper_limit_open: is the upper limit open?
% figure_title: title of plot, e.g. '(a)'
% figure_name: name for saving, e.g. UK_BFI
%
% OUTPUT
% plot and saved figure
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

if nargin < 3
    error('Not enough input arguments.')
end

p = inputParser;

addRequired(p, 'xdata', ...
    @(A) isnumeric(A) && (size(A,1)==1 || size(A,2)==1))
addRequired(p, 'ydata', ...
    @(phi) isnumeric(phi) && (size(phi,1)==1 || size(phi,2)==1))
addRequired(p, 'w', ...
    @(w) isnumeric(w) && length(w)==1)

addParameter(p, 'attribute', NaN(size(A)), @(x) isnumeric(x) || islogical(x))
addParameter(p, 'attribute_name', @ischar)
addParameter(p, 'ID', NaN(size(A)), @isnumeric)
addParameter(p, 'A_confInt', NaN(size(A)), @isnumeric)
addParameter(p, 'phi_confInt', NaN(size(phi)), @isnumeric)
addParameter(p, 'x_limits', [0 1.2], @(x) isnumeric(x) && length(x)==2)
addParameter(p, 'y_limits', [0 140], @(x) isnumeric(x) && length(x)==2)
addParameter(p, 'colour_scheme', 'Spectral', @ischar)
addParameter(p, 'flip_colour_scheme', false, @islogical)
addParameter(p, 'c_limits', [], @(x) isnumeric(x) && length(x)==2)
addParameter(p, 'c_lower_limit_open', false, @islogical)
addParameter(p, 'c_upper_limit_open', false, @islogical)
addParameter(p, 'figure_title', '', @ischar)
addParameter(p, 'figure_name', '', @ischar)

parse(p, A, phi, w, varargin{:})

attribute = p.Results.attribute;
attribute_name = p.Results.attribute_name;
ID = p.Results.ID;
A_confInt = p.Results.A_confInt;
phi_confInt = p.Results.phi_confInt;
x_limits = p.Results.x_limits;
y_limits = p.Results.y_limits;
colour_scheme = p.Results.colour_scheme;
flip_colour_scheme = p.Results.flip_colour_scheme;
c_limits = p.Results.c_limits;
c_lower_limit_open = p.Results.c_lower_limit_open;
c_upper_limit_open = p.Results.c_upper_limit_open;
figure_title = p.Results.figure_title;
figure_name = p.Results.figure_name;

% create colormap
if flip_colour_scheme
    colour_mat = flip(brewermap(10,colour_scheme));
else
    colour_mat = brewermap(10,colour_scheme);
end

a_range = logspace(-5,1,100);
A_theory = a_range./sqrt(a_range.^2 + w.^2);
phi_theory = acos(A_theory)./w;

% outer envelope (gamma-distribution with parameter=2)
A_mat1_outer = a_range./sqrt(a_range.^2 + w.^2);
A_mat2_outer  = a_range./sqrt(a_range.^2 + w.^2);
A_theory_mat_outer  = A_mat1_outer.*A_mat2_outer;
phi_theory_mat_outer  = acos(A_mat1_outer) + acos(A_mat2_outer);

f1 = figure('Name',attribute_name,'NumberTitle','off','pos',[10 10 350 250]);
hold on
% grid on
plot(A_theory,phi_theory,'color',[0.7 0.7 0.7],'linewidth',1.5)
plot(A_theory_mat_outer,phi_theory_mat_outer./w,'--','linewidth',1.5,'color',[0.7 0.7 0.7]); % serial reservoirs
errorbar(A,phi,...
    -phi_confInt,phi_confInt,...
    -A_confInt,A_confInt,'.','color','black')
scatter(A,phi,25,attribute,'filled','MarkerEdgeColor','k')%
xlim(x_limits)
ylim(y_limits)
ylabel('Phase shift [days]')
xlabel('Amplitude ratio [-]')

colormap(colour_mat)
c = colorbar;
title(c,attribute_name)
x1=get(gca,'position');
c_pos=get(c,'Position');
c_pos(3)=10/400;
% x(1)=0.95;
set(c,'Position',c_pos)
set(gca,'position',x1)
if ~isempty(c_limits)
    caxis(c_limits)
    if c_lower_limit_open
        c.TickLabels{1} = ['<' c.TickLabels{1}];
    end
    if c_upper_limit_open
        c.TickLabels{end} = ['>' c.TickLabels{end}];
    end
end

title(figure_title)

% update cursor
dcm_obj = datacursormode(figure(f1));
set(dcm_obj,'UpdateFcn',{@myupdatefcn,ID,[1:length(ID)]})

% print correlations
disp(['Correlation between amplitude ratio and ',attribute_name,' (Pearson/Spearman):'])
disp([num2str(0.01*round(100*corr(A,attribute,'rows','complete','type','Pearson'))),'/',...
    num2str(0.01*round(100*corr(A,attribute,'rows','complete','type','Spearman'))),''])
disp(['Correlation between phase shift and ',attribute_name,' (Pearson/Spearman):'])
disp([num2str(0.01*round(100*corr(phi,attribute,'rows','complete','type','Pearson'))),'/',...
    num2str(0.01*round(100*corr(phi,attribute,'rows','complete','type','Spearman'))),''])

% save fig
set(f1,'Units','Inches');
position = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[position(3),position(4)]);
fig_name_raw = strcat('seasonal_signatures_',figure_name);
fig_name = regexprep(fig_name_raw,'[^a-zA-Z0-9]','');
path_name = './Seasonal_signatures_paper_public/Images';
fig_path = strcat(path_name,'\',fig_name);
print(f1,fig_path,'-dpdf','-r600');

end
