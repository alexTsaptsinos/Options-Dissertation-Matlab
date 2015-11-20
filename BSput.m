function [BSprice] = BSput(strike,maturity,interest,volatility,initialPrice);
%% BLACK-SCHOLES CALCULATOR FOR EUROPEAN PUT
% This function will take in parameters for a european put option and will
% spit out the BS solution.

% Now calculate the BS explicit value
d1 = (log(initialPrice/strike) + (interest + volatility^2/2)*maturity)/(volatility*sqrt(maturity));
d2 = (log(initialPrice/strike) + (interest - volatility^2/2)*maturity)/(volatility*sqrt(maturity));
d2 = normcdf(-d2)*strike*exp(-interest*maturity);
d1 = normcdf(-d1)*initialPrice;

BSprice = d2 - d1;

