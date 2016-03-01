%% Anderson Broadie Simulation with variate reduction techniques attempt 2
% This simulation will calculate an upper and lower bound of an American
% Put on a single asset following a geometric brownian motion
clear, clc;
tic
%% Set up variables
K = 40; % strike price
r = 0.06; % interest
T = 1; % maturity
s = 0.2; % volatility (sigma)
S0 = 36; % initial price
N = 5*10^3; % sample paths for upper bound
d = 32;%0; % number of timesteps
%N1 = 1*10^6;%2*10^6; % number of sample paths for lower bound
N2 = 500; % number of subpath loops at continuation
%N3 = 10000; % number of subpath loops at exercise
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

for i = 2:d+1;
    Z = randn(1,N);
    Z = [Z,-Z]; % create antithetic pairs
    S(i,:) = S(i-1,:).*exp((r - s^2/2)*dt + s*Z*sqrt(dt));
end

%% Calculate the payoff matrix
% h is a matrix of size timesteps x loops, so each column corresponds to
% the payoffs along a path at each time. Note that time 0 is not included
% so when matching with S there will be one rows difference.
h = max(K-S(2:d+1,:),0);

%% Calculate the European option value at each time step for each sample path
% does not contain time zero as can't exercise then anyway
europeanValues = zeros(d,2*N);
for i = 1:d
    europeanValues(i,:) = BSput(K,(d-i)*dt,r,s,S(i+1,:));
end


%% Build the indicator matrix which will tell us when to exercise
C = zeros(d,2*N); % no time zero
for i = 1:d-1
    % at each time (not time 0)
    subS = S(i+1,:); % path values at that time
    D = generateChoiceFunctions(subS,M,K,(d-i)*dt,r,s);
    %D = generateBasisFunctions(subS,M);
    C(i,:) = D*beta(i+1,:)';
end
I = (h >= C) & (h>0); % so I is d x N
I(d,:) = h(d,:) > 0;
V = max(h,C);

% now build martingale
mart = zeros(d,2*N); % no time 0

for i=1:d
    i
    subS = zeros(2*N2,2*N);
    subV = zeros(2*N2,2*N);
    subC = zeros(2*N2,2*N);
    subH = zeros(2*N2,2*N);
    
    for n=1:2*N2
    
        Z = randn(1,N);
        Z = [Z,-Z];
        subS(n,:) = S(i,:).*exp((r-s^2/2)*dt + s*Z*sqrt(dt));
        subH(n,:) = max(K-subS(n,:),0);
        if i==d
            subV(n,:) = subH(n,:);
        else
            subD = generateChoiceFunctions(subS(n,:),M,K,(d-i)*dt,r,s);
            %subD = generateBasisFunctions(subS(n,:),M);
            subC(n,:) = subD*beta(i+1,:)';
            subV(n,:) = exp(-r*dt*i).*max(subH(n,:),subC(n,:));
        end
    end
    
    %sum(subV,1)/(2*N2));
    diff = exp(-r*dt*i).*V(i,:) - mean(subV);
    if i==1
        mart(i,:) = diff;%exp(-r*dt*i).*V(i,:);
    else
        mart(i,:) = mart(i-1,:) + diff;
    end
    

        
end

for i = 1:d
    h(i,:) = exp(-r*dt*i).*h(i,:);
end

diff = h - mart;
maximums = max(diff);

upperBound = mean(maximums) %+ controlLowerBound
upperStdError = std(maximums)/sqrt(2*N)
upperRelativeStdError = abs(upperStdError/upperBound)*100;

% Construct CI
alpha = 0.05;
z = norminv(1-alpha/2);
CIlower = controlLowerBound - z*controlStdError;
CIupper = upperBound + z*sqrt(controlStdError^2 + upperStdError^2);
CI = [CIlower,CIupper]

endtime = toc


