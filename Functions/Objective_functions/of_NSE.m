function [NS] = of_NSE(Qobs, Qmod)
%OF_NSE Returns Nash-Sutcliffe Efficiency Index given simulated and
% observed values (Nash and Sutcliffe, 1970).
%
% INPUT
% Qobs: observed streamflow
% Qmod: modelled streamflow
%
% OUTPUT
% NS: Nash-Sutcliffe Efficiency Index
%
% References
% Nash, J.E. and Sutcliffe, J.V., 1970. River flow forecasting through 
% conceptual models part I—A discussion of principles. Journal of 
% Hydrology, 10(3), pp.282-290.
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

Num=sum((Qobs-Qmod).^2);
QobsAv=mean(Qobs);
Den=sum((Qobs-QobsAv).^2);
NS=1-Num/Den;

end

