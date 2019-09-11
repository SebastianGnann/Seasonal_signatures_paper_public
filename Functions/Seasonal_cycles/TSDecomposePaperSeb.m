%==========================================================================
%   function [f, gof, output, Xr, yX] = TSDecompose(X,n)
%==========================================================================
%
%   Function to identify the characteristics of an underlying seasonality
%   in a given dataset X.  X needs to be an (n x 2) matrix with column 1 a
%   list of matlab date numbers, and column 2 a set of daily observations
%   that contain an underlying seasonal cycle.
%
%==========================================================================
%
%   - Original code by Nicholas Howden
%	- Sebastian: added confidence intervals and changed a few small things
%   
%==========================================================================
%
%   uses function nanxcov: 
%   http://www.swissfluxnet.ch/eddycalc/html/nanxcov.html
%
%==========================================================================


function [paramsS, paramsE, Xr, yX, I, paramsS_confInt] = TSDecomposePaperSeb(X, nS)

    %   extract the vector of date numbers from the input
    d = X(:,1);
    
    %   find out the MATLAB first day number of the record
    d1 = d(1,1);
    
    %   extract the datevector for the first day
    d1v = datevec(d(1,1));
    
    %   get the MATLAB datenumber for the first day of the year
    ds = datenum(d1v(1),10,1); 
    
    %   find the phase required to correct the reference series
    %   to the start of the observed record in days
    pdd = d1 - ds;
    
    %   construct dummy variables for the reference time series
    x = d - ds;
    y = sin(2*pi*x/365);

    %   find the unbiased autocovariance function from the cross-covariance
    xc = nanxcov(X(:,2),y,'unbiased');
    start = (length(xc) - 1)/2 + 1;
    xc = xc(start:length(xc));
    
    
    %   find the variance of the time series
    nv = nanvar(X(:,2));
    
    %   convert this into the acf
    xc = xc./(nv*sqrt(2));
    
    %   fit the seasonal cycle to the cross-covariance function
    xcs = xc(1:nS);
    xxc = 0:1:nS-1;

    %   fit the seasonal curve to the dependence structure
    [f, gof, output] = fit(xxc',-xcs,'sin1');

    %   get the coefficients of the equation
    paramsS = coeffvalues(f);
    paramsS_confInt = reshape(confint(f),6,1);
    
    %   calculate the seaonal pattern and then the residuals
    %   (de-seasonalised time series)
    yX = 2*paramsS(1)*nv*sqrt(2)*(1+sin(2.*pi.*x./365+paramsS(3)+pi/2));
    Xr = X(:,2) - yX;
        
    % adjust parameters
    paramsS(1) = 2*paramsS(1)*(nv*sqrt(2));
    paramsS(3) = paramsS(3)+pi/2;
    paramsS_confInt(1:2) = 2*paramsS_confInt(1:2)*(nv*sqrt(2));
    paramsS_confInt(5:6) = paramsS_confInt(5:6)+pi/2;

    %   now look at the autocorrelation in the residuals after
    %   de-seasonalising the flow time series
    %   find the autocovariance of the remaining signal
    acv = nanxcov(Xr,Xr,'unbiased');
    start2 = (length(acv) - 1)/2 +1;
    acv = acv(start2:length(acv));

    %   find the acf
    acf = acv ./ nanvar(X(:,2));
    
    %   find the first negative value in the ACV
    I = find(acf < 0.01,1,'first');
    
    %   fix I if it is too small
    if I < 5; I = 5; end
    
    %   truncate this to the early lags and define a timebase
    acfs = acf(1:I-1);
    %acvs2 = acv2(1:n2);
    x = 0:1:I-2;
    
    %   fit the two exponential decay model to the residuals
    [fexp2_1, gof_exp2_1, output_exp2_1] = fit(x', acfs, 'exp2');
    
    paramsE = coeffvalues(fexp2_1);    
    
    %   check that the two exponentials are listed in the correct places in
    %   the matrix
    
    if paramsE(2) < paramsE(4)
    
        paramsE = [paramsE(3), paramsE(4), paramsE(1), paramsE(2)];
              
    end
    
      %plot all the outcomes of the decomposition

      %{
      h = figure('Name','Flow Decomposition','NumberTitle','off')
  
    subplot(2,2,[1 2])    
        plot(X(:,1),X(:,2),'b-')
        hold on
        plot(X(:,1),yX,'r-')
        grid on
        datetick('x')
        set(gca, 'FontName', 'Times')
        set(gca, 'FontSize', 16)
        xlabel('date','FontWeight','Bold')
        ylabel('streamflow (mm / day)','FontWeight','Bold')
    
    subplot(2,2,3)
        plot(f, xxc,-xcs)
        grid on
        set(gca, 'FontName', 'Times')
        set(gca, 'FontSize', 16)
        xlabel('lag','FontWeight','Bold')
        ylabel('cross-covariance','FontWeight','Bold')    
    
    subplot(2,2,4)

        plot(x', acfs, 'r.')
        hold on
        plot(fexp2_1, 'k')
        grid on
        set(gca, 'FontName', 'Times')
        set(gca, 'FontSize', 16)
        xlabel('lag','FontWeight','Bold')
        ylabel('covariance','FontWeight','Bold')
     %}
end