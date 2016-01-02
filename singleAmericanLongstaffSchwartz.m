function [LSMlowerbound,europeanValue,LSMstdError,totaltime,relativeStdError] = singleAmericanLongstaffSchwartz(S0);
%% LSM Lower Bound Simulation
% This simulation will calculate a lower bound of an American Put on a
% single asset following a geometric brownian motion
%clc,clear all;
tic
%% Set up variables
K = 40; % strike price
r = 0.06; % interest
T = 1; % maturity
s = 0.2; % volatility (sigma)
%S0 = 36; % initial price
N = 1*10^5; % sample paths pairs for coefficients
N2 = 1*10^5; % sample paths pairs for valuation
d = 50; % number of timesteps
M = 4; % number of basis functions
dt = T/d; % size of each timestep

%% First find European value
europeanValue = BSput(K,T,r,s,S0)

%% Get regression coefficients
[beta,lowerBound,lowerBoundStdError] = LSMregressioncoefficients(K,r,T,s,S0,N,d,M);
coefficientsTime = toc;

%% Generate new sample paths for lower bound
% Generate all the new sample paths in a matrix S of size (timesteps + 1) x
% loops, so each column corresponds to a different path
S = zeros(d+1,N2);

% the first entry in each row will be the initial price
S(1,:) = S0;

for i = 2:d+1;
    Z = randn(1,N2);
    S(i,:) = S(i-1,:).*exp((r - s^2/2)*dt + s*Z*sqrt(dt));
end

%% Calculate the payoff matrix
% h is a matrix of size timesteps x loops, so each column corresponds to
% the payoffs along a path at each time NOT DISCOUNTED BACK TO TIME 0. Note
% that time 0 is not included so when matching with S there will be one
% rows difference.
h = max(K-S(2:d+1,:),0);
%h = K-S(2:d+1,:);
% %set the final payoff
% h(d,:) = exp(-r*T)*max(K - S(d+1,:),0);
% 
% for i = 1:d
%     h(i,:) = exp(-r*dt*i)*max(K - S(i+1,:),0);
% end


%% Compute the continuation matrix
C = zeros(d,N2); % continuation matrix. no time zero, but time d
for i = 1:d-1
    % CHOOSE either basis functions or choice ones
    %D = generateBasisFunctions(S(i+1,:),M);
    D = generateChoiceFunctions(S(i+1,:),M,K,(d-i)*dt,r,s);
    
    C(i,:) = D*beta(i+1,:)';
end

controlVariate = zeros(1,N2); % stores the control variate values
Y = zeros(1,N2);
for i = 1:N2
    indx = find(h(:,i) >= C(:,i) & h(:,i)>0,1);
    if indx
        Y(i) = exp(-r*dt*indx)*h(indx,i);
        controlVariate(i) = exp(-r*dt*indx)*BSput(K,(d-indx)*dt,r,s,S(indx+1,i));
    end
    %Y(i) = h(find(h(:,i) >= C(:,i),1),i);
end

LSMlowerbound = mean(Y)
LSMstdError = std(Y)/sqrt(N2)
relativeStdError = abs(LSMstdError/LSMlowerbound)*100;
totaltime = toc
end

