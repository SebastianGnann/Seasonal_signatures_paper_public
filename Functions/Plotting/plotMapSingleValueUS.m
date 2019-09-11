function [] = plotMapSingleValueUS(lat,lon,z,varargin)
%PLOTMAPSINGLEVALUEUS Plots map with dots coloured according to a single
% attribute.
% Options:
%   - colour dots according to third variable
%   - plot confidence intervals
%   - various plotting options, e.g. axes limits
%   (- save plot as PDF)
%
% INPUT
% lat: latitude
% lon: longitude
% z: attribute to be coloured in, e.g. BFI
% attribute_name: name of attribute
% ID: catchment ID
% colour_scheme: name of colour scheme
% flip_colour_scheme: flip colour scheme?
% c_limits: limits of colour axis, e.g. [0 1]
% c_lower_limit_open: is the lower limit open? 
% c_upper_limit_open: is the upper limit open?
% figure_title: title of plot, e.g. '(a)'
% figure_name: name for saving, e.g. CAMELS_BFI
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

addRequired(p, 'latitude', ...
    @(lat) isnumeric(lat) && (size(lat,1)==1 || size(lat,2)==1))
addRequired(p, 'longitude', ...
    @(lon) isnumeric(lon) && (size(lon,1)==1 || size(lon,2)==1))
addRequired(p, 'attribute', ...
    @(z) isnumeric(z) || islogical(z))

addParameter(p, 'attribute_name', @ischar)
addParameter(p, 'ID', NaN(size(z)), @isnumeric)
addParameter(p, 'colour_scheme', 'Spectral', @ischar)
addParameter(p, 'flip_colour_scheme', false, @islogical)
addParameter(p, 'c_limits', [0 1], @(x) isnumeric(x) && length(x)==2)
addParameter(p, 'c_lower_limit_open', false, @islogical)
addParameter(p, 'c_upper_limit_open', false, @islogical)
addParameter(p, 'figure_title', '', @ischar)
addParameter(p, 'figure_name', '', @ischar)

parse(p, lat, lon, z, varargin{:})

attribute_name = p.Results.attribute_name;
ID = p.Results.ID;
colour_scheme = p.Results.colour_scheme;
flip_colour_scheme = p.Results.flip_colour_scheme;
c_limits = p.Results.c_limits;
c_lower_limit_open = p.Results.c_lower_limit_open;
c_upper_limit_open = p.Results.c_upper_limit_open;
figure_title = p.Results.figure_title;
figure_name = p.Results.figure_name;

% get rid of NaN points
% zNaN = isnan(z);
% latitude(zNaN) = [];
% longitude(zNaN) = [];
% z(zNaN) = [];
index = [1:length(ID)]';
index(ID<1000000) = [];
ID(ID<1000000) = [];
% ID(zNaN) = []; index(zNaN) = [];
% ID(zNaN) = []; index(zNaN) = [];

fig1 = figure('Name',attribute_name,'NumberTitle','off','pos',[10 10 600 300]); 
ax = axesm('MapProjection','mercator','MapLatLimit',[24 52],'MapLonLimit',[-130 -65]);
states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow(ax, states, ...
    'DisplayType','polygon','DefaultFaceColor','white','DefaultEdgeColor','black')
hold on
% grid on

% create colormap
if flip_colour_scheme
    colour_mat = flip(brewermap(10,colour_scheme)); %
else
    colour_mat = (brewermap(10,colour_scheme));
end

% plot
scatterm(lat,lon,25,z,'filled') 
scatterm(lat(isnan(z)),lon(isnan(z)),'x k','linewidth',1.2);
xlabel('Latitude [km]'); ylabel('Longitude [km]')
axis equal
colormap(colour_mat)
c=colorbar;
title(c,attribute_name)
x1=get(gca,'position');
x=[0.8 0.2 0.01 0.25];
set(c,'Position',x)
set(gca,'position',x1)
set(gca,'Visible','off')
caxis(c_limits)
if c_lower_limit_open
    c.TickLabels{1} = ['<' c.TickLabels{1}];
end
if c_upper_limit_open
    c.TickLabels{end} = ['>' c.TickLabels{end}];
end
title(figure_title)

% update cursor
dcm_obj = datacursormode(figure(fig1));
set(dcm_obj,'UpdateFcn',{@myupdatefcn,ID,index})

% save fig
set(fig1,'Units','Inches');
pos = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos(3),pos(4)]);
fig_name_raw = strcat('map_','_',figure_name);
fig_name = regexprep(fig_name_raw,'[^a-zA-Z0-9]','');
path_name = './Seasonal_signatures_paper_public/Images';
fig_path = strcat(path_name,'\',fig_name);
print(fig1,fig_path,'-dpdf','-r500');

end

