%% Evaluate MARRMoT model runs and plot signature spaces
%
%   - change settings in script to get different plots (which model and
%   which evaluation metric to plot)
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

% clc
% clear all
% close all

%% add directories for functions to path
% addpath(genpath('./Data_processing'));

if exist('./BrewerMap') == 7
    addpath(genpath('./BrewerMap'));
else
    error('BrewerMap toolbox needed. Can be download from https://github.com/DrosteEffect/BrewerMap and should be in a folder named BrewerMap in the same directory.')
end

if exist('./MARRMoT') == 7
    addpath(genpath('./MARRMoT'));
else
    error('MARRMoT toolbox needed. Can be download from https://github.com/wknoben/MARRMoT and should be in a folder named MARRMoT in the same directory.')
end

%% specify which model runs to evaluate

% check if results files are unzipped
% if exist('./Seasonal_signatures_paper_public/Data_and_results/Results_MARRMoT') ~= 7   
%    disp('Unzipping MARRMoT results...')
%    unzip('./Seasonal_signatures_paper_public/Data_and_results/Results_MARRMoT.zip','./Seasonal_signatures_paper_public/Data_and_results')
% end

% load catchment ID and aridity index
results_UK_struc = load('seasonal_signatures_UK.mat');
results_UK = results_UK_struc.seasonal_signatures_UK;
ID_UK = results_UK.ID;
AI_UK = results_UK.AI;

% model specifications
model_list = ["m_05_ihacres_7p_1s",...
    "m_07_gr4j_4p_2s"];

model_id = 1; % change model: 1 for IHACRES, 2 for GR4J
model_name = char(model_list(model_id));
numStore = str2double(model_name(end-1)); % number of stores
if model_name(end-5) ~= '_'
    numParam = str2double(model_name(end-5:end-4));
else
    numParam = str2double(model_name(end-4)); % number of parameters
end

fprintf('Model name: %s \n', model_name);

% chose which signature to plot
signature_name = 'BFI'; % 'BFI' 'Amplitude ratio' 'Phase shift [days]'
        
% number of parameter sets
n_samples = 20000; % 2000, 5000, 10000, 20000
parRange = feval([model_name,'_parameter_ranges']); % parameter ranges

% specify subset of catchments
load('inSubset.mat');

% prepare plots
f1 = figure('pos',[10 10 300 150],'visible','off'); %300 150
hax1 = axes;
hold on
f2 = figure('pos',[10 10 350 250],'visible','off');
hax2 = axes;
hold on

for i=1:length(ID_UK) 
    
    % check whether record is part of subset
    if  inSubset(i)
        
        % load catchment data
        ID = ID_UK(i);

        % load results
        str_results = ...
            strcat('./Seasonal_signatures_paper_public/Data_and_results/Results_MARRMoT/',...
            model_name,'_NR_',num2str(n_samples),'_ID_',num2str(ID),'.mat');
        load(str_results);
        
        %% plot results
        % fig 1
        % define signature to be plotted
        if strcmp(signature_name,'BFI')
            signature = GLUE_results.BFI_UKIH;
            signature_min = 0; signature_max = 1;
            if model_id == 1; title_name = '(a)'; end
            if model_id == 2; title_name = '(b)'; end
        elseif strcmp(signature_name,'Amplitude ratio')
            signature = GLUE_results.amplitude_ratio;
            signature_min = 0; signature_max = 1.2;
            if model_id == 1; title_name = '(c)'; end
            if model_id == 2; title_name = '(d)'; end
        elseif strcmp(signature_name,'Phase shift [days]')
            signature = GLUE_results.phase_shift;
            signature_min = 0; signature_max = 200;
            if model_id == 1; title_name = '(e)'; end
            if model_id == 2; title_name = '(f)'; end
        else
            disp('Wrong signature name.')
        end
        
        % estimate distribution via kernel density estimation
        [density,temp_val] = ksdensity(signature);
        colour_mat = flip(brewermap(10,'Spectral')); 
        % get approximate colours for aridity index
        col_fac = floor(AI_UK(i)./(1.1*max(AI_UK))*10);
        plot(hax1,temp_val,density,'LineWidth',1.5,'color',colour_mat(col_fac,:))
        
        % fig 2
        % scatter plot 
        scatter(hax2,GLUE_results.amplitude_ratio(GLUE_results.KGE>-.41),...
            GLUE_results.phase_shift(GLUE_results.KGE>-.41),0.005,...
            GLUE_results.BFI_UKIH(GLUE_results.KGE>-.41))
        
    end
    
end

% fig 1
plotHelper1(f1,signature_min,signature_max,signature_name,title_name,model_name)

% fig2
plotHelper2(f2,model_id,model_name)

% get colourbar 
% plotHelperColourbar()

function plotHelper1(...
    f1,signature_min,signature_max,signature_name,title_name,model_name)
% helper function for amplitude vs. phase shift scatter plot
figure(f1)
xlim([signature_min signature_max])
xlabel(signature_name);
ylabel('Probability density')
set(gca,'yticklabel',[])
title(title_name)
colour_mat = flip(brewermap(10,'Spectral')); 
colormap(colour_mat)

set(f1,'Units','Inches');
position = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[position(3),position(4)]);
fig_name_raw = strcat('Model_results_pdf_',signature_name,model_name);
fig_name = regexprep(fig_name_raw,'[^a-zA-Z0-9]','');
path_name = './Seasonal_signatures_paper_public/Images';
fig_path = strcat(path_name,'\',fig_name);
print(f1,fig_path,'-dpdf','-r600');

end

function plotHelper2(f2,model_id,model_name)
% helper function for amplitude vs. phase shift scatter plot

% plot adjustments
figure(f2)
ylim([0 200])
xlim([0 1.2])
xlabel('Amplitude ratio [-]')
ylabel('Phase shift [days]')
c = colorbar;
colour_mat = (brewermap(10,'YlGnBu')); 
colormap(colour_mat)
x1=get(gca,'position');
x=get(c,'Position');
x(3)=10/400;
% x(1)=0.95;
set(c,'Position',x)
set(gca,'position',x1)
caxis([0 1])
title(c,'BFI')

% plot lines for single reservoir and two reservoirs in parallel
w = 2*pi/365;
a_range = logspace(-5,1,100);
% single reservoir
A_theory = a_range./sqrt(a_range.^2 + w.^2);
phi_theory = acos(A_theory)./(w);
% outer envelope serial reservoirs (gamma-distribution with parameter=2)
A_mat1_serial = a_range./sqrt(a_range.^2 + w.^2);
A_mat2_serial  = a_range./sqrt(a_range.^2 + w.^2);
A_theory_mat_serial  = A_mat1_serial.*A_mat2_serial;
phi_theory_mat_serial  = acos(A_mat1_serial) + acos(A_mat2_serial);
plot(A_theory,phi_theory,'color',[0.7 0.7 0.7],'linewidth',1.5) 
plot(A_theory_mat_serial,phi_theory_mat_serial./w,'--',...
    'linewidth',1.5,'color',[0.7 0.7 0.7]); 

if model_id == 1
    title('(a)')
elseif model_id == 2
    title('(b)')
end

set(f2,'Units','Inches');
position = get(f2,'Position');
set(f2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[position(3),position(4)]);
fig_name_raw = strcat('Model_results_seasonal_signatures_',model_name);
fig_name = regexprep(fig_name_raw,'[^a-zA-Z0-9]','');
path_name = './Seasonal_signatures_paper_public/Images';
fig_path = strcat(path_name,'\',fig_name);
print(f2,fig_path,'-dpdf','-r300');

end

function plotHelperColourbar()
% plot colourbar - needs to be adjusted and cropped manually
figure
c = colorbar('southoutside'); 
colour_mat = flip(brewermap(10,'Spectral')); 
colormap(colour_mat)
title(c,'E_p/P')
set(gca,'ColorScale','log')
caxis([0.2 2])
c.Ticks = [0.25 0.5 1 2];
c.TickLabels = [0.25 0.5 1 2];
end