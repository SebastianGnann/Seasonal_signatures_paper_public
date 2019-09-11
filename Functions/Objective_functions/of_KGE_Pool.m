function [val,c,w] = of_KGE_Pool(obs,sim,varargin)
% OF_KGE_POOL Calculates non-parametric Kling-Gupta Efficiency of simulated 
% streamflow (Pool et al, 2018). Ignores time steps with -999 values.
%
% INPUT
% obs       - time series of observations       [nx1]
% sim       - time series of simulations        [nx1]
% varargin  - optional weights of components    [3x1]
%
% OUTPUT
% val       - objective function value          [1x1]
% c         - components [r,alpha,beta]         [3x1]
% w         - weights    [wr,wa,wb]             [3x1]
%
% References
% Sandra Pool, Marc Vis & Jan Seibert (2018) Evaluating model performance:
% towards a non-parametric variant of the Kling-Gupta efficiency, 
% Hydrological Sciences Journal, 63:13-14, 1941-1953, 
% DOI: 10.1080/02626667.2018.1552002
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk
%
% Acknowledgements: Wouter Knoben

%% check inputs and set defaults
if nargin < 2
    error('Not enough input arguments')
elseif nargin > 3
    error('Too many inputs.')    
    
elseif nargin == 2
    w = [1,1,1]; % no weights specified, use defaults
    
elseif nargin == 3
    if size(varargin{1}) == [1,3] | size(varargin{1}) == [3,1] % check weights variable for size
        w = varargin{1}; % apply weights if size = correct
    else
        error('Weights should be a 3x1 or 1x3 vector.') % or throw error
    end
end    

% check time series size and rotate one if needed
if checkTimeseriesSize(obs,sim) == 0
    error('Time series not of equal size.')
    
elseif checkTimeseriesSize(obs,sim) == 2
    sim = sim'; % 2 indicates that obs and sim are the same size but have different orientations
end

%% check for missing values
% -999 is used to denote missing values in observed data, but this is later
% scaled by area. Therefore we check for all negative values, and ignore those.
idx = find(obs >= 0);   

%% calculate components
c(1) = corr(obs(idx),sim(idx),'Type','Spearman'); % r: rank correlation
FDC_sim = sort(sim/sum(sim)); % calculate normalised FDCs
FDC_obs = sort(obs/sum(obs));
c(2) = 1-0.5*sum(abs(FDC_sim - FDC_obs)); % alpha: normalised FDC
c(3) = mean(sim(idx))/mean(obs(idx)); % beta: bias 

%% calculate value
val = 1-sqrt((w(1)*(c(1)-1))^2 + (w(2)*(c(2)-1))^2 + (w(3)*(c(3)-1))^2); % weighted KGE

end