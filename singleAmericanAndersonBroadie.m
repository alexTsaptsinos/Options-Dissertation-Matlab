%% Anderson Broad Simulation
% This simulation will calculate an upper and lower bound of an American
% Put on a single asset following a geometric brownian motion
tic
%% Set up variables
K = 100; % strike price
r = 0.06; % interest
T = 0.5; % maturity
s = 0.4; % volatility (sigma)
S0 = 100; % initial price
N = 1500; % sample paths for upper bound
d = 50; % number of timesteps
N1 = 2*10^6; % number of sample paths for lower bound
N2 = 10000; % number of subpath loops at continuation
N3 = 10000; % number of subpath loops at exercise
dt = T/d; % size of each timestep

%% First find European value
europeanValue = BSput(K,T,r,s,S0)

%% Calculate Lower Bound
% launch N1 paths that stop according to the stopping rule
stoppedValues = zeros(N1,1);
for j = 1:N1
    S = 110;
    i=0;
    while (S > 85) && (i < d)
        i = i+1;
        Z = randn;
        S = S*exp((r-s^2/2)*dt + s*Z*sqrt(dt));
    end
    if i == d
        % didn't exercise
        stoppedValues(j,1) = 0;
    else
        % exercised
        stoppedDiscountedExerciseValue = exp(-r*i*dt)*max(K - S,0);
        stoppedValues(j,1) = stoppedDiscountedExerciseValue;
    end
end

LowerBound = mean(stoppedValues)
LBStandardError = std(stoppedValues)/sqrt(N1)
lbtime = toc

%% Generate sample paths
% Generate all the new sample paths in a matrix S of size (timesteps +
% 1) x loops, so each column corresponds to a different path
S = zeros(d+1,N);

% the first entry in each row will be the initial price
S(1,:) = S0;

for i = 2:d+1;
    Z = randn(1,N);
    S(i,:) = S(i-1,:).*exp((r - s^2/2)*dt + s*Z*sqrt(dt));
end

%% Calculate the payoff matrix
% h is a matrix of size timesteps x loops, so each column corresponds to
% the payoffs along a path at each time. Note that time 0 is not included
% so when matching with S there will be one rows difference.
h = zeros(d,N);
% set the final payoff
h(d,:) = exp(-r*T)*max(K - S(d+1,:),0);

for i = 1:d
    time = i*dt;
    h(i,:) = exp(-r*time)*max(K - S(i+1,:),0);
end


%% Ok now lets go through the algorithm
% We implement a simple stopping rule - we stop if the price of the option
% goes under 85


%% Calculate all the discounted lower bound values + expectations needed
maximums = zeros(N,1); % to store the max difference for each path
for j = 1:N
    S1 = S(:,j); % extract the subpath we want
    I_t = S1 < 85;
    L = zeros(d,1); % matrix to hold the discounted lower bound values
    E = zeros(d,1); % matrix to hold the expected values at exercise points
    for k = 1:d
        % look at the indicator at time k
        if I_t(k) == 0
            % continuation
            % launch N2 subpaths starting from Sk stopped according to tau_k to calculate Lk/Bk
            subS = zeros(d-k+1,N2);
            subS(1,:) = S1(k+1);
            for i = 2:d-k+1
                Z = randn(1,N2);
                subS(i,:) = subS(i-1,:).*exp((r - s^2/2)*dt + s*Z*sqrt(dt));
            end
            
            % now we have the subpaths find the first time we stop and
            % calculate the average
            stoppedValues = zeros(N2,1); % to store the stopped discounted values
            for i = 1:N2
                row = find(subS(:,i) < 85,1); % index of first stopping time
                if isempty(row) == 1
                    stoppedValues(i,1) = 0;
                else
                    stoppedValue = subS(row,i);
                    stoppedDiscountedExerciseValue = exp(-r*(row+k-1)*dt)*max(K - stoppedValue,0);
                    stoppedValues(i,1) = stoppedDiscountedExerciseValue;
                end
            end
            LkBk = mean(stoppedValues);
            L(k,j) = LkBk;
        else
            % exercise
            % first set Lk/Bk = hk/Bk
            LkBk = exp(-r*k)*S1(k+1);
            L(k,j) = LkBk;
            
            % now launch N3 subpaths starting from Sk stopped according to
            % tau_(k+1) to calculate Ek[L(k+1)/B(k+1)]
            subS = zeros(d-k+1,N3);
            subS(1,:) = S1(k+1);
            for i = 2:d-k+1
                Z = randn(1,N3);
                subS(i,:) = subS(i-1,:).*exp((r - s^2/2)*dt + s*Z*sqrt(dt));
            end
            
            % now we have the subpaths find the first time we stop past time k
            stoppedValues = zeros(N3,1); % to store the stopped discounted values
            for i = 1:N3
                row = find(subS(:,i) < 85,2); % index of second stopping time, ie first stopping time after time k
                if length(row) < 2
                    stoppedValues(i,1) = 0;
                else
                    row = row(2);
                    stoppedValue = subS(row,i);
                    stoppedDiscountedExerciseValue = exp(-r*(row+k-1)*dt)*max(K - stoppedValue,0);
                    stoppedValues(i,1) = stoppedDiscountedExerciseValue;
                end
            end
            Ek = mean(stoppedValues); % MC approx of Ek[L(k+1)/B(k+1)]
            E(k,j) = Ek;
        end
        
    end
    % Now have calculated all the Lk/Bk and the reqd Ek[L(k+1)/B(k+1)]
    
    % Now we can calculate the martingale
    mart = zeros(d,1);
    mart(1,1) = L(1,1);
    for  k = 1:d-1
        if I_t(k) == 0
            % continuation
            mart(k+1,1) = mart(k,1) + L(k+1,1) - L(k,1);
        else
            % exercise
            mart(k+1,1) = mart(k,1) + L(k+1,1) - L(k,1) - E(k,1) + h(k,j);
        end
        
    end
    diff = h(:,j) - mart;

    maxDiff = max(diff);
    maximums(j,1) = maxDiff; 
end

UpperBound = LowerBound + mean(maximums)
UBStandardError = std(maximums)/sqrt(N)


endtime = toc















