function [BFI_LyneHollick] = calc_BFI_LyneHollick(Q_tot)
%CALC_BFI_LYNEHOLLICK: Calculate baseflow and BFI using Lyne-Hollick filter
%(Lyne and Hollick, 1979)
%
% INPUT
% Q_tot: total streamflow 
%
% OUTPUT
% BFI calculated using Lyne-Hollick filter
%
% References
% Lyne, V. and Hollick, M., 1979, September. Stochastic time-variable 
% rainfall-runoff modelling. In Institute of Engineers Australia National
% Conference (Vol. 1979, pp. 89-93). Barton, Australia: Institute of 
% Engineers Australia.
%
% ---
% 
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

% remove NaN values
if any(isnan(Q_tot))
    disp('NaN values are assigned the median value for baseflow calculation')
end
Q_tot(isnan(Q_tot)) = nanmedian(Q_tot);

% calculate baseflow by applying RDF forwards, backwards, and forwards again
a_LyneHollick = 0.925; % Nathan and McMahon, 1990
B_LyneHollick = LyneHollickFilter(Q_tot, a_LyneHollick);
B_LyneHollick = LyneHollickFilter(flip(B_LyneHollick), a_LyneHollick);
B_LyneHollick = LyneHollickFilter(flip(B_LyneHollick), a_LyneHollick);

% calculate BFI
BFI_LyneHollick = sig_BFI(B_LyneHollick,Q_tot);

end