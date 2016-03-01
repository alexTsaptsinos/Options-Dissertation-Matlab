function [europeanValue,martApproximation,martStdError,martRelativeStdError,martApproximation2,martStdError2,martRelativeStdError2,totaltime] = gameOptionDualDual(S0);
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
N = 1*10^5; % sample paths pairs for coefficients
N2 = 1*10^5; % sample paths pairs for valuation
N3 = 5*10^3; % sample path pairs for lower bound valuation
N4 = 5*10^3; % subloops
d = 100; % number of timesteps
M = 4; % number of basis functions
dt = T/d; % size of each timestep
delta = 5; % penalty payoff

%% Get regression coefficients
%[beta,approximation,approximationStdError] = gameOptionCoefficients(K,r,T,s,S0,N,d,M,delta);
[controlApproximation,europeanValue,controlStdError,beta] = gameOptionLSM(K,r,T,s,S0,N,N2,d,M,delta);

%% Generate new sample paths for lower bound
% Generate all the new sample paths in a matrix S of size (timesteps + 1) x
% loops, so each column corresponds to a different path
S = zeros(d+1,2*N3);

% the first entry in each row will be the initial price
S(1,:) = S0;
for i = 2:d+1;
    Z = randn(1,N3);
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

%% Compute the continuation matrix
C = zeros(d,2*N3); % continuation matrix. no time zero, but time d
for i = 1:d-1
    % CHOOSE either basis functions or choice ones
    D = generateBasisFunctions(S(i+1,:),M);
    %D = generateChoiceFunctions(S(i+1,:),M,K,(d-i)*dt,r,s);
    
    C(i,:) = D*beta(i+1,:)';
end

%% Build the indicator matrix which will tell us when to exercise
I = ((h >= C) & (h>0)) | (g <= C); % so I is d x N2
I(d,:) = 1; % writer always cancels at maturity
%J = ~I;
V = min(g,max(h,C)); % no time 0
clear D Z;

% now build martingale
mart = zeros(d,2*N3); % no time 0
mart(1,:) = exp(-r*dt).*V(1,:);

for i=1:d-1
    i
    % check for which paths we need to run sub simulations
    timePaths = S(i+1,:);
    %indexs = I(i,:);
    relaventTimePaths = timePaths(I(i,:));
    subloopsNeeded = length(relaventTimePaths);
    
    % where no exercise just put in LtBt
    diff = zeros(1,2*N3);
    tempV1 = V(i+1,:);
    tempV2 = V(i,:);
    diff(~I(i,:)) = exp(-r*dt*(i+1)).*tempV1(~I(i,:)) - exp(-r*dt*i).*tempV2(~I(i,:));

    means = zeros(1,subloopsNeeded);
    
    for n=1:subloopsNeeded
    
        Z = randn(1,N4);
        Z = [Z,-Z];
        subS = relaventTimePaths(1,n).*exp((r-s^2/2)*dt + s*Z*sqrt(dt));
        subH = max(K-subS,0);
        subG = subH + delta;
        if i==d-1
            means(1,n) = mean(subH);
        else
            %subD = generateChoiceFunctions(subS,M,K,(d-i-1)*dt,r,s);
            subD = generateBasisFunctions(subS,M);
            subC = (subD*beta(i+2,:)')';
            subV = exp(-r*dt*(i+1)).*min(subG,max(subH,subC));
            means(1,n) = mean(subV);
        end
    end
    
    diff(I(i,:)) = exp(-r*dt*(i+1)).*tempV1(I(i,:))- means;

    mart(i+1,:) = mart(i,:) + diff;

        
end
clear tempV1 tempV2 diff subD subC subH subV relaventTimePaths timePaths;

for i = 1:d
    h(i,:) = exp(-r*dt*i).*h(i,:);
    g(i,:) = exp(-r*dt*i).*g(i,:);
end


%% Calculate Bounds, sup then inf!
% now need to calculate R(s,t) and M(s,t) for each time point t = 1,...,d
% and s = 1,...,d
Rt = zeros(d,2*N3);
Rs = zeros(d,2*N3);
for s = 1:d
    martTemp = mart;
    Rt(1:s,:) = h(1:s,:); % equals t times as t<= s
    %size(Rt(s+1:d,:))
    %size(g(s,:))
    if s~= d
        for j = s+1:d
            Rt(j,:) = g(s,:);  % now t > s, so is stopped at s
            martTemp(j,:) = martTemp(s,:);
        end
    end
    tempDiff = Rt - martTemp;
    tempTMax = max(tempDiff);
    Rs(s,:) = tempTMax;
end

minimums = min(Rs);

controlApproximation
martApproximation = mean(minimums) + controlApproximation
martStdError = std(minimums)/sqrt(2*N3)
martRelativeStdError = abs(martStdError/martApproximation)*100;

%% Now do the other way round, inf then sup!
Rt = zeros(d,2*N3);
Rs = zeros(d,2*N3);
for t = 1:d
    martTemp = mart;
    if t ~= 1        
        Rs(1:t-1,:) = g(1:t-1,:); % equals s times
    end
    %size(Rt(s+1:d,:))
    %size(g(s,:))
    for j = t:d
        Rs(j,:) = h(s,:);  % now s > t, so is stopped at t
        martTemp(j,:) = martTemp(s,:);
    end
    
    tempDiff = Rs - martTemp;
    tempSMin = min(tempDiff);
    Rt(t,:) = tempSMin;
end

maximums = max(Rt);

martApproximation2 = mean(maximums) + controlApproximation
martStdError2 = std(maximums)/sqrt(2*N3)
martRelativeStdError2 = abs(martStdError2/martApproximation2)*100;

    
totaltime = toc



end

