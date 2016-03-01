%% Belomestny, Bender and Schoenmakers algorithm
% This simulation will calculate an upper and lower bound of an American
% Put on a single asset following a geometric brownian motion
clear all, clc;
tic
%% Set up variables
K = 40; % strike price
r = 0.06; % interest
T = 1; % maturity
s = 0.2; % volatility (sigma)
S0 = 36; % initial price
N = 1*10^5; % sample paths for upper bound
N1 = 1*10^2;
d = 50; % number of timesteps
M = 4; % number of basis functions
dt = T/d; % size of each timestep

%% First find European value
%europeanValue = BSput(K,T,r,s,S0);

%% Calculate Lower Bound and Regression Coefficients
% Utilise the LSM method to find a lower bound and give us regression
% coefficients which we can then use to define an exercise policy.
% Remember beta includes 0 but not d
[controlLowerBound,europeanValue,controlStdError,totaltime,relativeStdError,beta] = singleAmericanLSMAntithetic(S0);
lbtime = toc;

%% Generate sample paths
% Generate all the new sample paths in a matrix S of size (timesteps +
% 1) x loops, so each column corresponds to a different path
S = zeros(d+1,2*N);

% the first entry in each row will be the initial price
S(1,:) = S0;
% Simulate N1 independent samples of the Bronwian increments
bMIntervals = sqrt(dt)*randn(d,2*N);

for i = 2:d+1;
    Z = bMIntervals(i-1,:);
    %Z = [Z,-Z]; % create antithetic pairs
    S(i,:) = S(i-1,:).*exp((r - s^2/2)*dt + s*Z);
end

%% Calculate the European option value at each time step for each sample path
% does not contain time zero as can't exercise then anyway
europeanValues = zeros(d,2*N);
for i = 1:d
    europeanValues(i,:) = BSput(K,(d-i)*dt,r,s,S(i+1,:));
end

%% Calculate the payoff matrix
% h is a matrix of size timesteps x loops, so each column corresponds to
% the payoffs along a path at each time. Note that time 0 is not included
% so when matching with S there will be one rows difference.
h = max(K-S(2:d+1,:),0);


%% Build the value matrix
C = zeros(d,2*N);
for i = 1:d-1
    % at each time (not time 0)
    subS = S(i+1,:); % path values at that time
    D = generateChoiceFunctions(subS,M,K,(d-i)*dt,r,s);
    C(i,:) = D*beta(i+1,:)';
end
V = max(h,C);


% Now we do another least squares regression to find Z
betaZ = zeros(d,M); % matrix to store the regression coefficients includes time zero but not time d
martingale = zeros(d,1); % no times 0

for i = 1:d
    
    D = generateChoiceFunctions(S(i+1,:),M,K,(d-i)*dt,r,s);
    
    Y = (bMIntervals(i,:)./dt).*V(i,:).*exp(-r*i*dt);
    Y=Y';

    betaSub = D\Y;
    betaZ(i,:) = betaSub;
end

% now have regression we calcuate the martingale
bMIntervals = sqrt(dt)*randn(d,2*N);
for i = 2:d+1;
    Z = bMIntervals(i-1,:);
    %Z = [Z,-Z]; % create antithetic pairs
    S(i,:) = S(i-1,:).*exp((r - s^2/2)*dt + s*Z);
end
h = max(K-S(2:d+1,:),0);

martAlmost = zeros(d,2*N);
for i = 1:d
    
    D = generateChoiceFunctions(S(i+1,:),M,K,(d-i)*dt,r,s);
    z = D*betaZ(i,:)';
    
    martTemp = z'.*bMIntervals(i,:);
    martAlmost(i,:) = martTemp;

end

mart = cumsum(martAlmost);

diff = h-mart;

maximums = max(diff);

upperEstimate = mean(maximums)





