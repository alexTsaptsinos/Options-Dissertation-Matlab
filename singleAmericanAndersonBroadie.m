%% Anderson Broadie Simulation
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
N = 10; % sample paths for upper bound
d = 50; % number of timesteps
%N1 = 1*10^6;%2*10^6; % number of sample paths for lower bound
N2 = 10000; % number of subpath loops at continuation
N3 = 10000; % number of subpath loops at exercise
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
% set the final payoff
% h(d,:) = max(K - S(d+1,:),0);
%
% for i = 1:d
%     time = i*dt;
%     h(i,:) = max(K - S(i+1,:),0);
% end

%% Build the indicator matrix which will tell us when to exercise
C = zeros(d,N);
for i = 1:d-1
    % at each time (not time 0)
    subS = S(i+1,:); % path values at that time
    D = generateBasisFunctions(subS,M);
    C(i,:) = D*beta(i+1,:)';
end
I = h >= C; % so I is d x N
I(d,:) = h(d,:) > 0;

%% Calculate all the discounted lower bound values + expectations needed
maximums = zeros(N,1); % to store the max difference for each path
for j = 1:N
    S1 = S(:,j); % extract the subpath we want
    %j
    L = zeros(d,1); % matrix to hold the discounted lower bound values
    E = zeros(d,1); % matrix to hold the expected values at exercise points
    for k = 1:d-1
        % look at the indicator at time k
        if I(k,j) == 0
            % continuation
            % launch N2 subpaths starting from Sk stopped according to our exercise policy to calculate Lk/Bk
            subS = zeros(d-k,N2);
            
            % take the first step
            startVector = S(k+1,j)*ones(1,N2);
            Z = randn(1,N2);
            subS(1,:) = startVector.*exp((r - s^2/2)*dt + s*Z*sqrt(dt));
            
            for i = 2:d-k;
                Z = randn(1,N2);
                subS(i,:) = subS(i-1,:).*exp((r - s^2/2)*dt + s*Z*sqrt(dt));
            end
            
            % So have generated the subpaths. Now want to look along each
            % subpath and stop it at the first time the exercise value is
            % greater than the continuation value.
            subh = max(K-subS,0);
            
            if k == d-1
                subI = subh > 0;
                [sel c]= max(subI~=0, [], 1 );
                idx = sub2ind(size(subh),c,[1:length(subh)]);
                subEx = subh(idx); % exercise values
                
                subEx2 = (exp(-r*c*dt).*subEx)'; % discounted exercise values
                L(k,1) = mean(subEx2);
                
            else
                subC = zeros(d-k,N2);
                for i = 1:d-k-1
                    subD = generateBasisFunctions(subS(i,:),M);
                    subC(i,:) = subD*beta(i+1,:)';
                end
                subC(d-k,:) = subh(d-k,:);
                
                subI = subh > subC;
                [sel c]= max(subI~=0, [], 1 );
                idx = sub2ind(size(subh),c,[1:length(subh)]);
                subEx = subh(idx); % exercise values
                
                subEx2 = (exp(-r*c*dt).*subEx)'; % discounted exercise values
                L(k,1) = mean(subEx2);
                
            end
            
            %% STUCK HERE
            
            
            
            
            
            
            
            %             stoppedValues = zeros(N2,1);
            %             for l = 1:N2
            %                 subS = S1(k+1);
            %                 % CHANGE THIS SO NOT LOOPING BUT JUST LIKE IN LSMREGRESSION
            %                 % SO THAT WE LOOK AT ALL AT ONCE IN A MATRIX
            %                 for i = k+1:d
            %                     Z = randn;
            %                     subS = subS*exp((r-s^2/2)*dt + s*Z*sqrt(dt));
            %                     tempExValue = max(K - subS,0);
            %
            %                     if i == d
            %                         if tempExValue > 0
            %                             % exercise
            %                             stoppedValues(l,1) = exp(-r*i*dt)*tempExValue;
            %                         else
            %                             % don't exercise
            %                             stoppedValues(l,1) = 0;
            %                         end
            %                     else
            %                         tempD = generateBasisFunctions(subS,M);
            %                         tempContValue = tempD*beta(i,:)';
            %                         if tempExValue >= tempContValue
            %                             % exercise
            %                             stoppedValues(l,1) = exp(-r*i*dt)*tempExValue;
            %                             break
            %                         end
            %                     end
            %                 end
            %
            %             end
            %             L(k,1) = mean(stoppedValues);
            
        else
            % exercise
            % first set Lk/Bk = hk/Bk
            L(k,1) = exp(-r*k*dt)*h(k,j);
            if k < d
                % now launch N3 subpaths starting from Sk stopped according to
                % tau_(k+1) to calculate Ek[L(k+1)/B(k+1)]
                subS = zeros(d-k,N3);
                
                % take the first step
                startVector = S(k+1,j)*ones(1,N3);
                Z = randn(1,N3);
                subS(1,:) = startVector.*exp((r - s^2/2)*dt + s*Z*sqrt(dt));
                
                for i = 2:d-k;
                    Z = randn(1,N3);
                    subS(i,:) = subS(i-1,:).*exp((r - s^2/2)*dt + s*Z*sqrt(dt));
                end
                
                % So have generated the subpaths. Now want to look along each
                % subpath and stop it at the first time the exercise value is
                % greater than the continuation value.
                subh = max(K-subS,0);
                
                if k == d-1
                    subI = subh > 0;
                    [sel c]= max(subI~=0, [], 1 );
                    idx = sub2ind(size(subh),c,[1:length(subh)]);
                    subEx = subh(idx); % exercise values
                    
                    subEx2 = (exp(-r*c*dt).*subEx)'; % discounted exercise values
                    E(k,1) = mean(subEx2);
                    
                else
                    subC = zeros(d-k,N2);
                    for i = 1:d-1-k
                        subD = generateBasisFunctions(subS(i,:),M);
                        subC(i,:) = subD*beta(i+1,:)';
                    end
                    subC(d-k,:) = subh(d-k,:);
                    
                    subI = subh > subC;
                    [sel c]= max(subI~=0, [], 1 );
                    idx = sub2ind(size(subh),c,[1:length(subh)]);
                    subEx = subh(idx); % exercise values
                    
                    subEx2 = (exp(-r*c*dt).*subEx)'; % discounted exercise values
                    E(k,1) = mean(subEx2);
                    
                end
                
                
                
                
                %                 stoppedValues = zeros(N3,1);
                %                 for l = 1:N3
                %                     subS = S1(k+1);
                %                     for i = k+1:d
                %                         Z = randn;
                %                         subS = subS*exp((r-s^2/2)*dt + s*Z*sqrt(dt));
                %                         tempExValue = max(K - subS,0);
                %
                %                         if i == d
                %                             if tempExValue > 0
                %                                 % exercise
                %                                 stoppedValues(l,1) = exp(-r*i*dt)*tempExValue;
                %                             else
                %                                 % don't exercise
                %                                 stoppedValues(l,1) = 0;
                %                             end
                %                         else
                %                             tempD = generateBasisFunctions(subS,M);
                %                             tempContValue = tempD*beta(i,:)';
                %                             if tempExValue >= tempContValue
                %                                 % exercise
                %                                 stoppedValues(l,1) = exp(-r*i*dt)*tempExValue;
                %                                 break
                %                             end
                %                         end
                %                     end
                %                 end
                %                 E(k,1) = mean(stoppedValues); % MC approx of Ek[L(k+1)/B(k+1)]
            end
        end
        
    end
    % Just need to calculate at k = d.
    if I(d,j) == 0
        % continue
        L(d,1) = 0;
    else
        % exercise
        L(d,1) = exp(-r*d*dt)*h(d,j);

    end
    % We have calculated all the Lk/Bk and the reqd Ek[L(k+1)/B(k+1)]
    
    % Now we can calculate the martingale
    mart = zeros(d,1);
    mart(1,1) = L(1,1);
    for  k = 1:d-1
        if I(k,j) == 0
            % continuation
            mart(k+1,1) = mart(k,1) + L(k+1,1) - L(k,1);
        else
            % exercise
            mart(k+1,1) = mart(k,1) + L(k+1,1) - L(k,1) - E(k,1) + exp(-r*k*dt)*h(k,j);
        end
        
    end
    diff = exp(-r*k*dt).*h(:,j) - mart;
    
    maxDiff = max(diff);
    maximums(j,1) = maxDiff;
end

upperBound = lowerBound + mean(maximums)
upperBoundStdError = std(maximums)/sqrt(N)

% Construct CI
alpha = 0.05;
z = norminv(1-alpha/2);
CIlower = lowerBound - z*lowerBoundStdError;
CIupper = upperBound + z*sqrt(lowerBoundStdError^2 + upperBoundStdError^2);
CI = [CIlower,CIupper]

endtime = toc















