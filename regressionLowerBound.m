%% LSM Method
% Function to calculate the regression coefficients using the LSM method

%% Set up variables
strike = 100;
interest = 0.06;
maturity = 0.5;
volatility = 0.4;
initialPrice = 80;
loops = 10000; % number of trajectories
timesteps = 100;
timestep = maturity/timesteps;
numberBasisFunctions = 5;

[regressionCoefficients] = LSMregressioncoefficients(strike,interest,maturity,volatility,initialPrice,loops,timesteps,numberBasisFunctions);

%% Estimate option Value using optimal stopping rule via 2nd round of MC

% Calculate new trajectories
for i = 2:timesteps+1;
    Z = randn(1,loops);
    samplePaths(i,:) = samplePaths(i-1,:).*exp((interest - volatility^2/2)*timestep + volatility*Z*sqrt(timestep));
end

% Calculate new payoff matrix h
% final value
h(timesteps,:) = exp(-interest*maturity)*max(strike - samplePaths(timesteps+1,:),0);

for t = 1:timesteps
    time = t*timestep;
    h(t,:) = exp(-interest*time)*max(strike - samplePaths(t+1,:),0);
end

% Now calculate the optimal stopping value
valueFunctions = zeros(timesteps,loops);
continuationValues2 = zeros(timesteps,loops);

for n = 1:timesteps
    if (n==timesteps)
        valueFunctions(n,:) = h(n,:);
    else
        basisFunctions = generateBasisFunctions(samplePaths(n+1,:),numberBasisFunctions);
        continuationValues2(n,:) = (regressionCoefficients(:,n)')*basisFunctions;
        valueFunctions(n,:) = max(h(n,:),continuationValues2(n,:));
    end
    
end

stoppingRule = zeros(1,loops); % the ith entry containing the optimal stopping value for the ith trajectory

for n=1:loops
    stoppingRule(n) = h(find(h(:,n)>=continuationValues2(:,n),1),n);
end

optionValue = sum(stoppingRule)/loops
standardError = std(stoppingRule)/sqrt(loops)
variance = (sum(stoppingRule.^2)/loops - optionValue^2)/loops
    
    