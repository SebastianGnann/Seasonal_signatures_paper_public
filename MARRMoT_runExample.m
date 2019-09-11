%% run MARRMoT - example file
% Run MARRMoT in a Monte Carlo sampling scheme for n-thousand parameter 
% sets. Parameters are sampled using Latin Hypercube sampling 
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

if exist('./MARRMoT/Octave') == 7
    error('You are currently using MARRMoTs Octave version. Please REMOVE the folder ./MARRMoT/Octave.')
end

%% load catchments
if exist('data_Example_struc') % check if data is already loaded
else
    Example_data_struc = load('./Seasonal_signatures_paper_public/Data_and_results/Example_data_new.mat');
end

Example_data = Example_data_struc.Example_data_new;
ID_Example = Example_data.ID;
P_Example = Example_data.P_daily; % [Matlab_date data [mm]]
PET_Example = Example_data.PET_daily; % [Matlab_date data [mm]]
Q_Example = Example_data.Q_daily; % [Matlab_date data [mm]]
T_Example = PET_Example; % dummy T time series - not used in models without snow routine

%% trim time series (from Oct 1989 to Sep 2009)
date_start = [1999,10,1]; % change to 1999 for example catchments
date_start = datenum(date_start);
date_end = [2009,09,30];
date_end = datenum(date_end);

P_Example = getSubPeriod(date_start,date_end,P_Example,ID_Example);
PET_Example = getSubPeriod(date_start,date_end,PET_Example,ID_Example);
Q_Example = getSubPeriod(date_start,date_end,Q_Example,ID_Example);
T_Example = getSubPeriod(date_start,date_end,T_Example,ID_Example);

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

% Latin hypercube sampling
% intitialise random seed?
rng(15)
n_samples = 20; % 2000 5000 10000 20000
sample_points = lhsdesign(n_samples,numParam);
parRange = feval([model_name,'_parameter_ranges']); % parameter ranges
% scale Latin hypercube with parameter ranges
sample_parameters = sample_points.*(parRange(:,2)-parRange(:,1))'+parRange(:,1)';

for i=1:length(ID_Example) % loop over catchments
    
    tic
    
    % load catchment data
    ID = ID_Example(i);
    P = P_Example{i};
    PET = PET_Example{i};
    Q = Q_Example{i};
    T = T_Example{i};
    
    % initialise vectors to store results
    NSE = NaN((n_samples),1);
    KGE = NaN((n_samples),1);
    KGE_Pool = NaN((n_samples),1);
    Q_mean = NaN((n_samples),1);
    Q_median = NaN((n_samples),1);
    Q_var = NaN((n_samples),1);
    Q_skewness = NaN((n_samples),1);
    AET_mean = NaN((n_samples),1);
    BFI_LH = NaN((n_samples),1);
    BFI_UKIH = NaN((n_samples),1);
    RR = NaN((n_samples),1);
    PQ_elasticity = NaN((n_samples),1);
    Q5 = NaN((n_samples),1);
    Q95 = NaN((n_samples),1);
    FDC_midslope = NaN((n_samples),1);
    HFD_mean = NaN((n_samples),1);
    HFI_mean = NaN((n_samples),1);
    high_Q_freq = NaN((n_samples),1);
    high_Q_dur = NaN((n_samples),1);
    low_Q_freq = NaN((n_samples),1);
    low_Q_dur = NaN((n_samples),1);
    no_Q_freq = NaN((n_samples),1);
    no_Q_dur = NaN((n_samples),1);
    Q_AC1 = NaN((n_samples),1);
    rising_limb_density = NaN((n_samples),1);
    phase_shift = NaN((n_samples),1);
    amplitude_ratio = NaN((n_samples),1);
    
    fprintf('Catchment nr (max 1192): %d \n', i);
    
    % run MARRMoT model with all parameter sets
    for j = 1:n_samples % use parfor to parallelise code
        disp(j)
        % warming up
        initial_storages = zeros(numStore,1);
        initial_storages = calculateInitialStorages(...
            P,T,PET,Q,....
            model_name,sample_parameters(j,:),initial_storages);
        % run model
        [Qmod,AETmod] = ...
            runMARRMoT(P,T,PET,Q,...
            model_name,sample_parameters(j,:),initial_storages,...
            false,false);
        
        % calculate signatures using modelled streamflow
        try
            NSE(j) = of_NSE(Q(:,2),Qmod);
            KGE(j) = of_KGE(Q(:,2),Qmod);
            KGE_Pool(j) = of_KGE_Pool(Q(:,2),Qmod);
            Q_mean(j) = nanmean(Qmod);
            Q_median(j) = nanmedian(Qmod);
            Q_var(j) = nanvar(Qmod);
            Q_skewness(j) = mean(Qmod)/median(Qmod);
            AET_mean(j) = nanmean(AETmod);
            BFI_LH(j) = calc_BFI_LyneHollick(Qmod);
            BFI_UKIH(j) = calc_BFI_UKIH(Qmod);
            RR(j) = mean(Qmod)/mean(P(:,2));
            PQ_elasticity(j) = sig_streamflow_elasticity(Qmod,P(:,2),Q(:,1));
            Q5(j) = sig_x_percentile(Qmod,5);
            Q95(j) = sig_x_percentile(Qmod,95);
            FDC_midslope(j) = sig_FDC_midslope(Qmod);
            HFD_mean(j) = sig_HFD_mean(Qmod,Q(:,1));
            HFI_mean(j) = sig_HFI_mean(Qmod,Q(:,1));
            high_Q_freq(j) = sig_high_Q_freq(Qmod);
            high_Q_dur(j) = sig_high_Q_dur(Qmod);
            low_Q_freq(j) = sig_low_Q_freq(Qmod);
            low_Q_dur(j) = sig_low_Q_dur(Qmod);
            no_Q_freq(j) = sig_no_Q_freq(Qmod);
            no_Q_dur(j) = sig_no_Q_dur(Qmod);
            Q_AC1_temp = autocorr(Qmod);
            Q_AC1(j) = Q_AC1_temp(2);
            rising_limb_density(j) = sig_rising_limb_density(Qmod);
            [amplitude_ratio(j),phase_shift(j)] = sig_seasonal(Qmod,P(:,2),PET(:,2),Q(:,1));
        catch
            fprintf('Error in calculating signatures. Parameter set %d ',j)
        end
    end
    
    % create results struc
    %         clear MC_results
    MC_results.NSE = NSE;
    MC_results.KGE = KGE;
    MC_results.KGE_Pool = KGE_Pool;
    MC_results.Q_mean = Q_mean;
    MC_results.Q_median = Q_median;
    MC_results.Q_var = Q_var;
    MC_results.Q_skewness = Q_skewness;
    MC_results.AET_mean = AET_mean;
    MC_results.BFI_LH = BFI_LH;
    MC_results.BFI_UKIH = BFI_UKIH;
    MC_results.RR = RR;
    MC_results.PQ_elasticity = PQ_elasticity;
    MC_results.Q5 = Q5;
    MC_results.Q95 = Q95;
    MC_results.FDC_midslope = FDC_midslope;
    MC_results.HFD_mean = HFD_mean;
    MC_results.HFI_mean = HFI_mean;
    MC_results.high_Q_freq = high_Q_freq;
    MC_results.high_Q_dur = high_Q_dur;
    MC_results.low_Q_freq = low_Q_freq;
    MC_results.low_Q_dur = low_Q_dur;
    MC_results.no_Q_freq = no_Q_freq;
    MC_results.no_Q_dur = no_Q_dur;
    MC_results.Q_AC1 = Q_AC1;
    MC_results.rising_limb_density = rising_limb_density;
    MC_results.amplitude_ratio = amplitude_ratio;
    MC_results.phase_shift = phase_shift;
    
    MC_results.parameter_set = sample_parameters;
    % TODO: new one for every model or always the same x-thousand ones?
    MC_results.ID = ID;
    MC_results.n_samples = n_samples;
    
    % save results
    str_run = strcat('./Seasonal_signatures_paper_public/Data_and_results/Results_MARRMoT_Example/',model_name,'_NR_',num2str(n_samples),'_ID_',num2str(ID),'.mat');
    save(str_run,'MC_results')
    
    toc
    
end
