function [X] = getSubPeriod(date_start,date_end,X,ID)
%GETSUBPERIOD Extract subperiod from time series, e.g. from 1 Oct 1989 to
% 30 Sept 1999
%
% INPUT
% date_start: start date as Matlab date integer
% date_end: end date as Matlab date integer
% X: time series (time,value)
% ID: catchment ID list
%
% OUTPUT
% trimmed time series (time,value)
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

% figure
for i = 1:length(ID)
    X_temp = X{i};
    date_temp = X_temp(:,1);
    index_start = find(date_temp == date_start);
    index_end = find(date_temp == date_end);
    X{i} = [X_temp(index_start:index_end,1) X_temp(index_start:index_end,2)];
end

end

