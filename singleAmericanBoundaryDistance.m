%% Anderson Broadie Simulation with boundary distance grouping
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
N = 5*10^4; % sample paths for upper bound
d = 32;%0; % number of timesteps
%N1 = 1*10^6;%2*10^6; % number of sample paths for lower bound
N2 = 1*10^3; % number of subpath loops at continuation
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
% europeanValues = zeros(d,2*N);
% for i = 1:d
%     europeanValues(i,:) = BSput(K,(d-i)*dt,r,s,S(i+1,:));
% end


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
%J = ~I;
V = max(h,C);
clear D Z;

% now build martingale
mart = zeros(d,2*N); % no time 0
mart(1,:) = exp(-r*dt).*V(1,:);

for i=1:d-1
    i
    % check for which paths we need to run sub simulations
    timePaths = S(i+1,:);
    %indexs = I(i,:);
    relaventTimePaths = timePaths(I(i,:));
    subloopsNeeded = length(relaventTimePaths);
    
    % where no exercise just put in LtBt
    diff = zeros(1,2*N);
    tempV1 = V(i+1,:);
    tempV2 = V(i,:);
    diff(~I(i,:)) = exp(-r*dt*(i+1)).*tempV1(~I(i,:)) - exp(-r*dt*i).*tempV2(~I(i,:));
    means = zeros(1,subloopsNeeded);
    
    for n=1:subloopsNeeded
    
        Z = randn(1,N2);
        Z = [Z,-Z];
        subS = relaventTimePaths(1,n).*exp((r-s^2/2)*dt + s*Z*sqrt(dt));
        subH = max(K-subS,0);
        if i==d-1
            means(1,n) = mean(subH);
        else
            subD = generateChoiceFunctions(subS,M,K,(d-i-1)*dt,r,s);
            %subD = generateBasisFunctions(subS(n,:),M);
            subC = (subD*beta(i+2,:)')';
            subV = exp(-r*dt*(i+1)).*max(subH,subC);
            means(1,n) = mean(subV);
        end
    end
    
    %sum(subV,1)/(2*N2));
    diff(I(i,:)) = exp(-r*dt*(i+1)).*tempV1(I(i,:))- means;
%     if i==1
%         mart(i,:) = diff;%exp(-r*dt*i).*V(i,:);
%     else
%         mart(i,:) = mart(i-1,:) + diff;
%     end
    mart(i+1,:) = mart(i,:) + diff;

        
end
clear tempV1 tempV2 diff subD subC subH subV relaventTimePaths timePaths;

for i = 1:d
    h(i,:) = exp(-r*dt*i).*h(i,:);
end

diff = h - mart;
maximums = max(diff);

controlVariate = zeros(1,2*N2); % stores the control variate values
Y = max(diff);
for i = 1:2*N
    indx = find(h(:,i) >= C(:,i) & h(:,i)>0,1);
    if indx
        controlVariate(i) = exp(-r*dt*indx)*BSput(K,(d-indx)*dt,r,s,S(indx+1,i));
    end
    %Y(i) = h(find(h(:,i) >= C(:,i),1),i);
end

meanControl = mean(controlVariate);
meanNormal = mean(Y);
covarianceEstimate = mean(Y.*controlVariate);
varianceEstimate = mean(controlVariate.^2);
controlVariateConstant = -(covarianceEstimate - meanNormal*meanControl)/(varianceEstimate - meanControl^2);
Z = Y + controlVariateConstant*(controlVariate - europeanValue);

testStd = std(maximums)/sqrt(2*N)
testMean = mean(maximums) + controlLowerBound
upperBound = mean(Z) + controlLowerBound
upperStdError = std(Z)/sqrt(2*N)
upperRelativeStdError = abs(upperStdError/upperBound)*100;

% Construct CI
alpha = 0.05;
z = norminv(1-alpha/2);
CIlower = controlLowerBound - z*controlStdError;
CIupper = upperBound + z*sqrt(controlStdError^2 + upperStdError^2);
CI = [CIlower,CIupper]

endtime = toc


