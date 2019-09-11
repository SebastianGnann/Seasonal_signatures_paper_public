%% Create a data file containing selected catchment attributes
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

% clc
% clear all
% close all

%% add directories for functions to path

if exist('./Data_processing') == 7
    addpath(genpath('./Data_processing')); % my personal data store - 
    % not necessary to run the code. Only the data-files are needed
end

%% load data and metadata
if exist('data_UK_struc') % check if data is already loaded
else
    data_UK_struc = load('./Data_processing/Data/UK_data_new.mat'); 
end
if exist('data_CAMELS_struc') % check if data is already loaded
else
    data_CAMELS_struc = load('./Data_processing/Data/CAMELS_data_new.mat'); 
end

data_UK = data_UK_struc.UK_data;
data_CAMELS = data_CAMELS_struc.CAMELS_data;

%% save relevant metadata
UK_metadata.ID = data_UK.ID;
UK_metadata.AI = data_UK.AI;
UK_metadata.isBenchmark = data_UK.isBenchmark;

CAMELS_metadata.ID = data_CAMELS.ID;
CAMELS_metadata.AI = data_CAMELS.AI;

save('./Seasonal_signatures_paper/Data_and_results/metadata_UK.mat','UK_metadata')
save('./Seasonal_signatures_paper/Data_and_results/metadata_CAMELS.mat','CAMELS_metadata')
