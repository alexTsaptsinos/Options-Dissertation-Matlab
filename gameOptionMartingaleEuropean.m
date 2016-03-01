function [europeanValue,martApproximation,martStdError,martRelativeStdError,martApproximation2,martStdError2,martRelativeStdError2,totaltime] = gameOptionMartingale(S0);
%% LSM Lower Bound Simulation with Antithetic Variate & Control Variate
% This simulation will calculate an upper bound of a game option with put
% payoff.
%clc,clear all;
tic
%% Set up variables
K = 100; % strike price
r = 0.06; % interest
T = 0.5; % maturity
s = 0.4; % volatility (sigma)
%S0 = 36; % initial price
N = 5*10^4; % sample path pairs for valuation
d = 500; % number of timesteps
M = 4; % number of basis functions
dt = T/d; % size of each timestep
delta = 5; % penalty payoff
trueValue = 2.54;

%% First find European value
europeanValue = BSput(K,T,r,s,S0)

%% Generate new sample paths for lower bound
% Generate all the new sample paths in a matrix S of size (timesteps + 1) x
% loops, so each column corresponds to a different path
S = zeros(d+1,2*N);

% the first entry in each row will be the initial price
S(1,:) = S0;
for i = 2:d+1;
    Z = randn(1,N);
    Z = [Z,-Z]; % create antithetic pairs
    S(i,:) = S(i-1,:).*exp((r - s^2/2)*dt + s*Z*sqrt(dt));
end

%% Calculate the payoff matrix
% h is a matrix of size timesteps x loops, so each column corresponds to
% the payoffs along a path at each time NOT DISCOUNTED BACK TO TIME 0. Note
% that time 0 is not included so when matching with S there will be one
% rows difference.
h = max(K-S(2:d+1,:),0);
g = h + delta; % penalty matrix
g(d,:) = h(d,:); % take off penalty at maturity

%% Now build discounted European option martingale
mart = zeros(d,2*N); % no time 0
mart(d,:) = h(d,:);

for i=1:d-1
    i
    % check for which paths we need to run sub simulations
    currentTimePaths = S(i+1,:);
    mart(i,:) = exp(-r*dt*i)*BSput(K,(d-i)*dt,r,s,currentTimePaths);
        
end

for i = 1:d
    h(i,:) = exp(-r*dt*i).*h(i,:);
    g(i,:) = exp(-r*dt*i).*g(i,:);
end

clear S;
%% Calculate Bounds, sup then inf!
% now need to calculate R(s,t) and M(s,t) for each time point t = 1,...,d
% and s = 1,...,d
%Rt = zeros(d,2*N);
Rs = zeros(d,2*N);
for s = 1:d
    s
%     martTemp = mart;
%     Rt(1:s,:) = h(1:s,:); % equals t times as t<= s
    
    
%     if s~= d
%         for j = s+1:d
%             Rt(j,:) = g(s,:);  % now t > s, so is stopped at s
%             martTemp(j,:) = martTemp(s,:);
%         end
%     end
    extraDiff = g(s,:) - mart(s,:);
    tempDiff = [h(1:s,:) - mart(1:s,:);extraDiff];
    %tempDiff = [tempDiff;extraDiff];
    tempTMax = max(tempDiff);
    Rs(s,:) = tempTMax;
    clear extraDiff tempDiff tempTMax;
end

minimums = min(Rs);

martApproximation = mean(minimums) + europeanValue
martStdError = std(minimums)/sqrt(2*N)
martRelativeStdError = abs(martStdError/martApproximation)*100
martRelativeError = abs((trueValue - martApproximation)./trueValue).*100

% %% Now do the other way round, inf then sup!
% Rt = zeros(d,2*N);
% %Rs = zeros(d,2*N);
% for t = 1:d
%     %martTemp = mart;
% %     if t ~= 1        
% %         Rs(1:t-1,:) = g(1:t-1,:); % equals s times
% %     end
%     %size(Rt(s+1:d,:))
%     %size(g(s,:))
% %     for j = t:d
% %         Rs(j,:) = h(t,:);  % now s > t, so is stopped at t
% %         martTemp(j,:) = martTemp(t,:);
% %     end
%     
%     extraDiff = h(t,:) - mart(t,:);
%     tempDiff = g(1:t-1,:) - mart(1:t-1,:);
%     tempDiff = [tempDiff;extraDiff];
%     tempSMin = min(tempDiff);
%     Rt(t,:) = tempSMin;
% end
% 
% maximums = max(Rt);
% 
% martApproximation2 = mean(maximums) + europeanValue
% martStdError2 = std(maximums)/sqrt(2*N)
% martRelativeStdError2 = abs(martStdError2/martApproximation2)*100;

    
totaltime = toc



end

