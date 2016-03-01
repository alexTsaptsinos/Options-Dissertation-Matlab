function [controlApproximation,europeanValue,controlStdError,relativeStdError,beta,totaltime] = gameOptionLSM(K,r,T,s,S0,N,N2,d,M,delta)
%% LSM Lower Bound Simulation with Antithetic Variate & Control Variate
% This simulation will calculate a lower bound of an American Put on a
% single asset following a geometric brownian motion
tic

%% Set up variables
%K = 100; % strike price
%r = 0.06; % interest
%T = 0.5; % maturity
%s = 0.4; % volatility (sigma)
%S0 = 36; % initial price
%N = 1*10^6; % sample paths pairs for coefficients
%N2 = 5*10^4; % sample paths pairs for valuation
%d = 32; % number of timesteps
%M = 4; % number of basis functions
dt = T/d; % size of each timestep
%delta = 5; % penalty payoff

%% First find European value
europeanValue = BSput(K,T,r,s,S0)

%% Get regression coefficients
[beta] = gameOptionCoefficients(K,r,T,s,S0,N,d,M,delta);

%% Generate new sample paths for lower bound
% Generate all the new sample paths in a matrix S of size (timesteps + 1) x
% loops, so each column corresponds to a different path
S = zeros(d+1,2*N2);

% the first entry in each row will be the initial price
S(1,:) = S0;

for i = 2:d+1;
    Z = randn(1,N2);
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
C = zeros(d,2*N2); % continuation matrix. no time zero, but time d
for i = 1:d-1
    % CHOOSE either basis functions or choice ones
    %D = generateBasisFunctions(S(i+1,:),M);
    D = generateChoiceFunctions(S(i+1,:),M,K,(d-i)*dt,r,s);
    
    C(i,:) = D*beta(i+1,:)';
end

controlVariate = zeros(1,2*N2); % stores the control variate values
Y = zeros(1,2*N2);
for i = 1:2*N2
    indx1 = find(h(:,i) >= C(:,i) & h(:,i)>0,1);
    indx2 = find(g(:,i) <= C(:,i),1);
    if isempty(indx1) == 1
        indx1 = d+1;
    end
    
    if isempty(indx2) == 1
        indx2 = d; % exercise at maturity always
    end
    %indx = min(indx1,indx2);
    if indx1 <= indx2
        Y(1,i) = exp(-r*dt*indx1)*h(indx1,i);
    else
        Y(1,i) = exp(-r*dt*indx2)*g(indx2,i);
    end
    %controlVariate(1,i) = exp(-r*dt*indx)*BSput(K,(d-indx)*dt,r,s,S(indx+1,i));
    %Y(i) = h(find(h(:,i) >= C(:,i),1),i);
end

meanControl = mean(controlVariate);
meanNormal = mean(Y);
covarianceEstimate = mean(Y.*controlVariate);
varianceEstimate = mean(controlVariate.^2);
controlVariateConstant = -(covarianceEstimate - meanNormal*meanControl)/(varianceEstimate - meanControl^2);
Z = Y + controlVariateConstant*(controlVariate - europeanValue);

% control variate is fucking up so just use normal
controlApproximation = mean(Y)
controlStdError = std(Y)/sqrt(2*N2)
relativeStdError = abs(controlStdError/controlApproximation)*100;
%controlApproximation = mean(Z)
%controlStdError = std(Z)/sqrt(2*N2)
totaltime = toc
end

