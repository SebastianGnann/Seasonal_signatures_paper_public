function [B] = LyneHollickFilter(Q,a)
%LYNEHOLLICKFILTER Digital filter as shown in Lyne and Hollick (1979).
%
% INPUT
% Q: streamflow
% a: filter parameter
%
% OUTPUT
% B: baseflow
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

n = length(Q);
f = NaN(n,1);
f(1) = Q(1) - min(Q); % initial condition, e.g. Su et al., 2016

for i=2:1:n
    f(i) = a*f(i-1) + (1+a)/2*(Q(i) - Q(i-1));
    if f(i)<0
        f(i) = 0;
    end
end

B = Q - f;
%B(B>Q) = Q(B>Q); % during or after filtering?

end

