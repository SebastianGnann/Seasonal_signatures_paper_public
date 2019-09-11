function [EI] = BudykoCurve(AI)
%BUDYKOCURVE "Original" Budyko curve (Budyko, 1974)
%
% INPUT
% AI: aridity index, PET/P
%
% OUTPUT
% EI: evaporative index, Ep/P
%
% References
% Budyko, M.I., Miller, D.H. and Miller, D.H., 1974. Climate and life 
% (Vol. 508). New York: Academic press.
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

EI = sqrt(AI.*tanh(1./AI).*(1 - exp(-AI)));

end

