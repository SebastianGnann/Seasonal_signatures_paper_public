function [A,phi,k,X_res,X_hat] = fitSineCurve(X,w)
%FITSINECURVE fit sine curve to time series
%
% INPUT
% X: time series (time,data)
% w: angular frequency
%
% OUTPUT
% A0: amplitude
% phi: phase shift
% k: offset
% X_res
% X_hat
%
% ---
%
% Sebastian Gnann (2019)
% sebastian.gnann@bristol.ac.uk

% create matrix
M = ones(length(X(:,1)),3);
M(:,2) = cos(w*X(:,1));
M(:,3) = sin(w*X(:,1));

% change NaN values to median for linear regression
X(isnan(X(:,2)),2) = nanmedian(X(:,2));

% solve equation system
b = M\X(:,2);

% get estimated sine curve parameters
phi = atan2(b(2),b(3)); % get unambigous value for phi
A = sqrt(b(2)^2+b(3)^2);
k = b(1);

% get estimated sine curve
X_hat = A*sin(w*X(:,1) + phi) + k;
X_res = X;
X_res(:,2) = X(:,2) - X_hat;

end

