%% SINGLE ASSET AMERICAN PUT
% We use methods from Anderson and Broadie to find the value of an American put option
% with a single underlying asset

%% Set up variables
strike = 100;
interest = 0.06;
maturity = 0.5;
volatility = 0.4;
initialPrice = 80;
loops = 5000;
timesteps = 100;

%% First find European value
europeanValue = BSput(strike,maturity,interest,volatility,initialPrice)