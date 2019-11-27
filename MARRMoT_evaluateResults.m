%% Evaluate MARRMoT model runs and plot signature spaces
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
if exist('./Seasonal_signatures_paper_public/Data_and_results/Results_MARRMoT') ~= 7
    disp('Unzipping MARRMoT results...')
    unzip('./Seasonal_signatures_paper_public/Data_and_results/Results_MARRMoT.zip','./Seasonal_signatures_paper/Data_and_results')
end

% load catchment ID and aridity index
results_UK_struc = load('seasonal_signatures_UK.mat');
results_UK = results_UK_struc.seasonal_signatures_UK;
ID_UK = results_UK.ID;
I_m_UK = results_UK.I_m;

% model specifications
model_list = ["m_05_ihacres_7p_1s",...
    "m_07_gr4j_4p_2s"];

% prepare fig1
f1 = figure('pos',[10 10 700 750],'visible','on'); %300 150

for model_id = 1:2
    
    % model_id = 2; % change model: 1 for IHACRES, 2 for GR4J
    model_name = char(model_list(model_id));
    numStore = str2double(model_name(end-1)); % number of stores
    if model_name(end-5) ~= '_'
        numParam = str2double(model_name(end-5:end-4));
    else
        numParam = str2double(model_name(end-4)); % number of parameters
    end
    
    fprintf('Model name: %s \n', model_name);
    
    % chose which signature to plot
    signature_name_list = ["Amplitude ratio", "Phase shift [days]", "BFI"]';
    
    for i_name = 1:3
        
        signature_name = signature_name_list(i_name);
        
        % number of parameter sets
        n_samples = 20000; % 2000, 5000, 10000, 20000
        parRange = feval([model_name,'_parameter_ranges']); % parameter ranges
        
        % specify subset of catchments
        load('inSubset.mat');
        
        % prepare fig 2
        if i_name > 1
        else
            f2 = figure('pos',[10 10 350 250],'visible','on');
            hax2 = axes;
            hold on
        end
        
        A_plot = [];
        phi_plot = [];
        BFI_plot = [];
        
        rng(11)
        % set random seed
        random_loop = randperm(length(ID_UK));
        
        max_density = 0;
        n_catchments = length(results_UK.ID(boolean(inSubset)));
        sig_min = NaN(n_catchments,1);
        sig_max = NaN(n_catchments,1);
        sig_med = NaN(n_catchments,1);
        sig_q25 = NaN(n_catchments,1);
        sig_q75 = NaN(n_catchments,1);
        sig_q01 = NaN(n_catchments,1);
        sig_q99 = NaN(n_catchments,1);
        sig_obs = NaN(n_catchments,1);
        sig_ID = NaN(n_catchments,1);
        sig_mat = NaN(n_catchments,n_samples);
        I_m = NaN(n_catchments,1);
        col_fac = NaN(n_catchments,1);
        count = 0;
        
        for i=random_loop
            
            try
                % check whether record is part of subset
                if  inSubset(i)
                    
                    % load catchment data
                    count = count+1;
                    sig_ID(count) = ID_UK(i);
                    
                    % load results
                    str_results = ...
                        strcat('./Seasonal_signatures_paper_public/Data_and_results/Results_MARRMoT/',...
                        model_name,'_NR_',num2str(n_samples),'_ID_',num2str(sig_ID(count)),'.mat');
                    load(str_results);
                    
                    %% plot results
                    % fig 1
                    % define signature to be plotted
                    if strcmp(signature_name,'BFI')
                        signature = MC_results.BFI_UKIH;
                        signature_min = 0; signature_max = 1;
                        sig_obs(count) = results_UK.BFI(i);
                        if model_id == 1; title_name = '(c)'; end
                        if model_id == 2; title_name = '(f)'; end
                    elseif strcmp(signature_name,'Amplitude ratio')
                        signature = MC_results.amplitude_ratio;
                        signature_min = 0; signature_max = 1.2;
                        sig_obs(count) = results_UK.A_vec(i);
                        if model_id == 1; title_name = '(a)'; end
                        if model_id == 2; title_name = '(d)'; end
                    elseif strcmp(signature_name,'Phase shift [days]')
                        signature = MC_results.phase_shift;
                        signature_min = 0; signature_max = 200;
                        sig_obs(count) = results_UK.phi_vec(i);
                        if model_id == 1; title_name = '(b)'; end
                        if model_id == 2; title_name = '(e)'; end
                    else
                        disp('Wrong signature name.')
                    end
                    
                    % estimate distribution via kernel density estimation
                    %         [density,temp_val] = ksdensity(signature);
                    %         max_density = max(max_density, max(density));
                    %         colour_mat = (brewermap(10,'Spectral'));
                    % get approximate colours for aridity index
                    %         col_fac = ceil((I_m_UK(i)+1)*5);
                    %         plot(hax1,temp_val,density,'LineWidth',1.5,'color',colour_mat(col_fac,:))
                    
                    % remove runs that are "really bad"
                    rem = boolean(MC_results.amplitude_ratio>0.01 & ...
                        MC_results.amplitude_ratio<1.2 & ...
                        MC_results.phase_shift<200);
                    % & MC_results.KGE>-.41
                    % & (MC_results.Q_mean./(results_UK.Q_mean(i)/365))>0.1);
                    
                    % get min, max and observed value
                    sig_mat(count,:) = signature;
                    signature = signature(rem);
                    sig_min(count) = min(signature);
                    sig_max(count) = max(signature);
                    sig_med(count) = median(signature);
                    sig_q25(count) = quantile(signature,0.25);
                    sig_q75(count) = quantile(signature,0.75);
                    sig_q01(count) = quantile(signature,0.01);
                    sig_q99(count) = quantile(signature,0.99);
                    I_m(count) = I_m_UK(i);
                    col_fac(count) = ceil((I_m(count)+1)*5);
                    
                    % tmp = MC_results.parameter_set(:,4);
                    A_plot = [A_plot; MC_results.amplitude_ratio(rem)];
                    phi_plot = [phi_plot; MC_results.phase_shift(rem)];
                    BFI_plot = [BFI_plot; MC_results.BFI_UKIH(rem)];
                    
                end
            catch
            end
            
        end
        
        % fig 2
        if i_name>1
        else
            shuffle_UK = randperm(length(A_plot))'; % shuffle plotting order
            scatter(hax2,A_plot(shuffle_UK),phi_plot(shuffle_UK),0.5,...
                BFI_plot(shuffle_UK),'filled')
            
            scatter(hax2,results_UK.A_vec(boolean(inSubset)),results_UK.phi_vec(boolean(inSubset)),...
                25,results_UK.BFI(boolean(inSubset)),'filled','markeredgecolor','k')
        end
        
        % fig 1
        if model_id == 1
            nr_plot = 0+i_name;
        else
            nr_plot = 3+i_name;
        end
        
        figure(f1)
        hax = subplot(2,3,nr_plot);
        hold on
        plotHelper1(...
            hax,signature_min,signature_max,signature_name,title_name,model_name,...
            sig_min,sig_max,sig_med,sig_q25,sig_q75,sig_q01,sig_q99,sig_obs,I_m,col_fac)
                
        % fig2
        if i_name>1
        else
            plotHelper2(f2,model_id,model_name)
        end
                
    end
end

% get colourbar and save fig
plotHelperColourbar(f1)

function plotHelper1(...
    hax,signature_min,signature_max,signature_name,title_name,model_name,...
    sig_min,sig_max,sig_med,sig_q25,sig_q75,sig_q01,sig_q99,sig_obs,I_m,col_fac)
% helper function for amplitude vs. phase shift scatter plot
% figure(f1)
colour_mat = (brewermap(10,'RdBu'));
[~,sorted_i] = sort(I_m); %sig_obs'
ind = 0;
isOutside = 0;
for i = sorted_i'
    try
        ind = ind+1;
        sig_range_lex = sig_min(i):(sig_q01(i)-sig_min(i))/1000:sig_q01(i);
        sig_range_low = sig_q01(i):(sig_q25(i)-sig_q01(i))/1000:sig_q25(i);
        sig_range_mid = sig_q25(i):(sig_q75(i)-sig_q25(i))/1000:sig_q75(i);
        sig_range_upp = sig_q75(i):(sig_q99(i)-sig_q75(i))/1000:sig_q99(i);
        sig_range_uex = sig_q99(i):(sig_max(i)-sig_q99(i))/1000:sig_max(i);
        sig_i_lex = ind*ones(size(sig_range_lex));
        sig_i_low = ind*ones(size(sig_range_low));
        sig_i_mid = ind*ones(size(sig_range_mid));
        sig_i_upp = ind*ones(size(sig_range_upp));
        sig_i_uex = ind*ones(size(sig_range_uex));
        scatter(hax,sig_obs(i),ind,10,'markeredgecolor',colour_mat(col_fac(i),:),...
            'linewidth',1.0)
        plot(hax,sig_range_lex,sig_i_lex,':','color',colour_mat(col_fac(i),:),'linewidth',1.0)
        plot(hax,sig_range_low,sig_i_low,'-','color',colour_mat(col_fac(i),:),'linewidth',1.0)
        plot(hax,sig_range_mid,sig_i_mid,'-','color',colour_mat(col_fac(i),:),'linewidth',2)
        plot(hax,sig_range_upp,sig_i_upp,'-','color',colour_mat(col_fac(i),:),'linewidth',1.0)
        plot(hax,sig_range_uex,sig_i_uex,':','color',colour_mat(col_fac(i),:),'linewidth',1.0)
        %         plot(sig_med(i),ind,'+','color',colour_mat(col_fac(i),:),...
        %             'linewidth',2)
        %         plot(sig_q25(i),ind,'+','color',colour_mat(col_fac(i),:),...
        %             'linewidth',1)
        %         plot(sig_q75(i),ind,'+','color',colour_mat(col_fac(i),:),...
        %             'linewidth',1)
        if sig_obs(i)<sig_min(i) || sig_obs(i)>sig_max(i)
            isOutside = isOutside+1;
        else
            scatter(hax,sig_obs(i),ind,10,'filled','markeredgecolor',colour_mat(col_fac(i),:),...
                'markerfacecolor',colour_mat(col_fac(i),:),'linewidth',1.0)
        end
    catch
    end
end

disp('Catchments outside range: ')
disp(isOutside)

xlim([signature_min signature_max])
xlabel(signature_name);
ylabel('Catchments')
ylim([0 ind+1])
set(gca,'yticklabel',[])
set(gca,'ytick',[])
title(title_name)

if strcmp(model_name,'m_07_gr4j_4p_2s')
    hax.Position = hax.Position + [0 0.03 0 0];
end
% colour_mat = flip(brewermap(10,'Spectral'));
% colormap(colour_mat)

% set(f1,'Units','Inches');
% position = get(f1,'Position');
% set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[position(3),position(4)]);
% fig_name_raw = strcat('Model_results_range_',signature_name,model_name);
% fig_name = regexprep(fig_name_raw,'[^a-zA-Z0-9]','');
% path_name = './Seasonal_signatures_paper/Images';
% fig_path = strcat(path_name,'\',fig_name);
% print(f1,fig_path,'-dpdf','-r900');

end

function plotHelper2(f2,model_id,model_name)
% helper function for amplitude vs. phase shift scatter plot

% plot adjustments
figure(f2)
ylim([0 140])
xlim([0 1.2])
xlabel('Amplitude ratio [-]')
ylabel('Phase shift [days]')
c = colorbar;
colour_mat = (brewermap(10,'YlGnBu'));
% colour_mat = (brewermap(10,'Spectral'));
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
    title('(a) IHACRES')
elseif model_id == 2
    title('(b) GR4J')
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

function plotHelperColourbar(f1)
% plot colourbar and model names
figure(f1)
c = colorbar('southoutside');
c.Position = [0.4115 0.03 0.21 0.01];
c.Label.Rotation = 45;
colour_mat = (brewermap(10,'RdBu'));
colormap(colour_mat)
title(c,'Moisture index')
% set(gca,'ColorScale','log')
caxis([-1 1])
c.Ticks = [-1 0 1];
c.TickLabels = [-1 0 1];

an1 = annotation('textarrow',[0.05 0.05],[0.816 0.816],'String',"IHACRES",...
    'HeadStyle','none','LineStyle','none','TextRotation',90,'FontSize',16);
an2 = annotation('textarrow',[0.05 0.05],[0.348 0.348],'String',"GR4J",...
    'HeadStyle','none','LineStyle','none','TextRotation',90,'FontSize',16);

set(f1,'Units','Inches');
position = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[position(3),position(4)]);
fig_name_raw = strcat('Model_results_range');
fig_name = regexprep(fig_name_raw,'[^a-zA-Z0-9]','');
path_name = './Seasonal_signatures_paper_public/Images';
fig_path = strcat(path_name,'\',fig_name);
print(f1,fig_path,'-dpdf','-r600');

end