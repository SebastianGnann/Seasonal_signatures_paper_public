function [] = plotBivariate(x,y,varargin)
%PLOTBIVARIATE Plots scatter plot of two variables.
% Optional:
%   - plot correlation
%   - plot fit
%   - plot histograms
%
% INPUT
% x: variable 1
% y: variable 2
% x_name: string with name of x
% y_name: string with name of y
% ID: if data cursor shall show ID input the respective array
% x_limit: axis limit for x-axis
% y_limit: axis limit for y-axis
% show_corr: boolean to specify whether to show correlation
% show_fit: boolean to specify whether to plot linear regression line
% show_hist: boolean to specify whether to plot histogram
% figure_title: title of plot, e.g. '(a)'
% figure_name: name for saving, e.g. UK_P_Q
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

p = inputParser;

addRequired(p, 'xdata', ...
    @(x) isnumeric(x) && (size(x,1)==1 || size(x,2)==1))
addRequired(p, 'ydata', ...
    @(y) isnumeric(y) && (size(y,1)==1 || size(y,2)==1))

addParameter(p, 'x_name', @ischar)
addParameter(p, 'y_name', @ischar)
addParameter(p, 'ID', NaN(size(x)), @isnumeric)
addParameter(p, 'x_limits', [min(x) max(x)], @(z) isnumeric(z) && length(z)==2)
addParameter(p, 'y_limits', [min(y) max(y)], @(z) isnumeric(z) && length(z)==2)
addParameter(p, 'show_corr', false, @islogical)
addParameter(p, 'show_fit', false, @islogical)
addParameter(p, 'show_hist', false, @islogical)
addParameter(p, 'figure_title', '', @ischar)
addParameter(p, 'figure_name', '', @ischar)

parse(p, x, y, varargin{:})

x_name = p.Results.x_name;
y_name = p.Results.y_name;
ID = p.Results.ID;
x_limits = p.Results.x_limits;
y_limits = p.Results.y_limits;
show_corr = p.Results.show_corr;
show_fit = p.Results.show_fit;
show_hist = p.Results.show_hist;
figure_title = p.Results.figure_title;
figure_name = p.Results.figure_name;

% get rid of NaN points
xNaN = isnan(x); x(xNaN) = []; y(xNaN) = [];
yNaN = isnan(y); x(yNaN) = []; y(yNaN) = [];
index = [1:length(ID)]';
ID(xNaN) = []; ID(yNaN) = [];
index(xNaN) = []; index(yNaN) = [];

% nice colours
colour_mat = brewermap(10,'YlGnBu');
if show_hist
    fig1 = figure('Name',figure_name,'NumberTitle','off','pos',[10 10 500 500]);
else
    fig1 = figure('Name',figure_name,'NumberTitle','off','pos',[10 10 300 280]);
end

hold on
grid on
title(figure_title)

% plot
if show_hist
    scatterhist(x,y,'Color','k','NBins',20)
    hold on
else
    scatter(x,y,25,'k') %,'filled'
end

Pearson_cor = corr(x,y,'Type','Pearson','rows','complete'); %,'rows','complete'
Kendall_cor = corr(x,y,'Type','Kendall','rows','complete');
Spearman_cor = corr(x,y,'Type','Spearman','rows','complete');

% do linear regression if requested
if show_fit
    P = polyfit(x,y,1);
    yfit = polyval(P,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    
    set(0,'CurrentFigure',fig1(1));
    plot(x,yfit,'-','color',colour_mat(end/2,:),'linewidth',3);
    equation = sprintf(...
        'y = %.2f x + %.2f \nR^2 = %.2f \nPearson = %.2f \nKendall = %.2f \nSpearman = %.2f',...
        P(1),P(2),rsq,Pearson_cor,Kendall_cor,Spearman_cor);
    if show_corr
        if show_hist
            annotation('textbox',[0.025, 0.035, 0.25, 0.2],...
                'string',equation,'EdgeColor','none')
        else
            TextLocation(equation,'Location',[0.52 0.31 0.46 0.17]);
        end
    end
elseif show_corr
    set(0,'CurrentFigure',fig1(1));
    equation = sprintf(...
        'Pearson = %.2f \nKendall = %.2f \nSpearman = %.2f',...
        Pearson_cor,Kendall_cor,Spearman_cor);
    TextLocation(equation,'Location',[0.52 0.18 0.46 0.17]);
end

xlim(x_limits)
ylim(y_limits)
% axis equal
xlabel(x_name)
ylabel(y_name)

% update cursor
dcm_obj = datacursormode(figure(fig1(1)));
set(dcm_obj,'UpdateFcn',{@myupdatefcn,ID,index})

% save fig
set(fig1,'Units','Inches');
position = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[position(3),position(4)]);
fig_name_raw = strcat('scatter_bivariate_',figure_name);
fig_name = regexprep(fig_name_raw,'[^a-zA-Z0-9]','');
path_name = './Seasonal_signatures_paper_public/Images';
fig_path = strcat(path_name,'\',fig_name);
print(fig1,fig_path,'-dpdf','-r0');

end

function hOut = TextLocation(textString,varargin)

l = legend(textString,varargin{:});
t = annotation('textbox');
t.String = textString;
t.Position = l.Position;
delete(l);
t.LineStyle = 'None';

if nargout
    hOut = t;
end
end