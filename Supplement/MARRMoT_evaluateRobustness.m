%% Evaluate MARRMoT model runs and test whether the results are robust
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

% clc
% clear all
% close all

%% add directories for functions to path

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
if exist('./Seasonal_signatures_paper/Data_and_results/Results_MARRMoT') ~= 7
    disp('Unzipping MARRMoT results...')
    unzip('./Seasonal_signatures_paper/Data_and_results/Results_MARRMoT.zip','./Seasonal_signatures_paper/Data_and_results')
end

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
for signature_name = ["KGE", "BFI", "Amplitude ratio", "Phase shift"]
    % signature_name = 'KGE'; % 'KGE' 'BFI' 'Amplitude ratio' 'Phase shift'
    
    % number of parameter sets
    n_samples_list = [2000 5000 10000 20000];
    
    results_cell = cell(4,40);
    
    for k = 1:length(n_samples_list)
        
        n_samples = n_samples_list(k);
        
        % specify subset of catchments
        load('inSubset.mat');
        
        i_subset = 0;
        
        for i=1:length(ID_UK) % loop over catchments
            
            % check whether record is part of subset
            if  inSubset(i)
                
                i_subset = i_subset + 1;
                
                % load catchment data
                ID = ID_UK(i);
                
                % load results (no BC3 folder)
                str_res = strcat('./Seasonal_signatures_paper_public/Data_and_results/Results_MARRMoT/',...
                    model_name,'_NR_',num2str(n_samples),'_ID_',num2str(ID),'.mat');
                load(str_res);
                
                % define signature to be plotted
                if strcmp(signature_name,'BFI')
                    signature = MC_results.BFI_UKIH;
                    title_name = '(b)';
                    signature_min = 0.0; signature_max = 1.0;
                elseif strcmp(signature_name,'Amplitude ratio')
                    signature = MC_results.amplitude_ratio;
                    title_name = '(c)';
                    signature_min = 0; signature_max = 1.5;
                elseif strcmp(signature_name,'Phase shift')
                    signature = MC_results.phase_shift;
                    title_name = '(d)';
                    signature_min = 0; signature_max = 200;
                elseif strcmp(signature_name,'KGE')
                    signature = MC_results.KGE;
                    title_name = '(a)';
                    signature_min = -0.5; signature_max = 1;
                else
                    disp('Wrong signature name.')
                end
                
                % store results
                rem = boolean(MC_results.amplitude_ratio>0.01 & ...
                    MC_results.amplitude_ratio<1.2 & ...
                    MC_results.phase_shift<200);
                signature(~rem) = NaN;
                results_cell{k,i_subset} = signature;
                
            end
            
        end
        
    end
    
    plotHelper(...
        results_cell,model_name,title_name,signature_name,signature_min,signature_max)
    
end

function plotHelper(...
    results_cell,model_name,title_name,signature_name,signature_min,signature_max)
% helper function for box-whisker plots

f1 = figure('pos',[10 10 300 200]);
hold on
t1 = reshape([results_cell{1,:}],40*2000,1);
t2 = reshape([results_cell{2,:}],40*5000,1);
t3 = reshape([results_cell{3,:}],40*10000,1);
t4 = reshape([results_cell{4,:}],40*20000,1);
t_help = [2000*ones(40*2000,1); 5000*ones(40*5000,1); 10000*ones(40*10000,1); 20000*ones(40*20000,1)];
boxplot([t1; t2; t3; t4],t_help,...
    'Notch','off','Labels',{'n = 2000','n = 5000','n = 10000','n = 20000'},...
    'Colors','k','OutlierSize',0.2,'Symbol','r o','Widths',0.3)
ylim([signature_min signature_max])
ylabel(signature_name)
title(title_name)

%% save figs
set(f1,'Units','Inches');
position = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[position(3),position(4)]);
fig_name_raw = strcat('Model_results_boxplot',signature_name,model_name);
fig_name = regexprep(fig_name_raw,'[^a-zA-Z0-9]','');
path_name = './Seasonal_signatures_paper_public/Images';
fig_path = strcat(path_name,'\',fig_name);
print(f1,fig_path,'-dpdf','-r600');

end