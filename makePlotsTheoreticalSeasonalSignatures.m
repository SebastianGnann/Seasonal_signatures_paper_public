%% Plot theoretical seasonal signature plots
%
%   - single linear reservoir
%   - linear reservoirs in series
%   - linear reservoirs in parallel with different parameters
%   - a combination thereof
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

%% specify period
period_1y = 365;
w_1y = 2*pi/period_1y;

%% plot figures
plotSingleReservoir(w_1y);
plotSerialReservoirs(w_1y);
plotParallelReservoirs(w_1y,'(a)');
plotParallelReservoirs(w_1y,'(b)');
plotParallelReservoirs(w_1y,'(c)');
plotParallelReservoirs2(w_1y);
plotSerialAndParallelReservoirs(w_1y);
