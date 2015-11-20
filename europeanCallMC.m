%% Monte Carlo for a European Option
% Even though we know real-value through BS we run MC simulations for a cal
% call European option and compare values. We let the underlying stock
% price have a GBM.
% dS(t)/S(t) = rdt + sdW(t)
% This has an explicit solution - the logarithm of the stock price is
% normally distributed - the stock price has a lognormal distribution.

%% SET UP FIXED VARIABLES
strikePrice = 100;
maturity = 0.5;
interestRate = 0.06;
volatility = 0.4;
initialPrice = 80;
loops = 100000;

%% BEGIN LOOPS
Cis = [];
for i = 1:loops;
    Z = randn;
    S = initialPrice*exp((interestRate - volatility^2/2)*maturity + volatility*sqrt(maturity)*Z);
    payoff = max(S - strikePrice,0);
    Ctemp = exp(-interestRate*maturity)*payoff;
    Cis = [Cis Ctemp];
end

% Final values
C = mean(Cis)
s_C = std(Cis)

% Now calculate the BS explicit value
tempValue = (log(initialPrice/strikePrice) + (interestRate + volatility^2/2)*maturity)/(volatility*sqrt(maturity));
tempValue = initialPrice*normcdf(tempValue);
tempValue2 = (log(initialPrice/strikePrice) + (interestRate - volatility^2/2)*maturity)/(volatility*sqrt(maturity));
tempValue2 = exp(-interestRate*maturity)*strikePrice*normcdf(tempValue2);
BCvalue = tempValue - tempValue2