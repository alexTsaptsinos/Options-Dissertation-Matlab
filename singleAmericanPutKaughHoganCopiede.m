%% SINGLE ASSET AMERICAN PUT
% We use methods from Anderson and Broadie to find the value of an American put option
% with a single underlying asset

%% Set up variables
strike = 100;
interest = 0.06;
maturity = 0.5;
volatility = 0.4;
initialPrice = 80;
loops = 10000;
timesteps = 100;
subpathLoops = 200;
timestep = maturity/timesteps;
numberBasisFunctions = 5;


%% First find European value
europeanValue = BSput(strike,maturity,interest,volatility,initialPrice)

%% Load regression coefficients
[regressionCoefficients] = LSMregressioncoefficients(strike,interest,maturity,volatility,initialPrice,loops,timesteps,numberBasisFunctions);

%% Generate sample paths
% Generate all the sample paths in a matrix of size loops x (timesteps +
% 1), so each column corresponds to a different path
samplePaths = zeros(timesteps+1,loops);

% the first entry in each row will be the initial price
samplePaths(1,:) = initialPrice;

for i = 2:timesteps+1;
    Z = randn(1,loops);
    samplePaths(i,:) = samplePaths(i-1,:).*exp((interest - volatility^2/2)*timestep + volatility*Z*sqrt(timestep));
end

%% Calculate the payoff matrix
h = zeros(timesteps,loops);
% set the final payoff
h(timesteps,:) = exp(-interest*maturity)*max(strike - samplePaths(timesteps+1,:),0);

for t = 1:timesteps
    time = t*timestep;
    h(t,:) = exp(-interest*time)*max(strike - samplePaths(t+1,:),0);
end

%% Calculate approximating martingale

valueFunctions = zeros(timesteps,loops);
continuationValues = zeros(timesteps,loops);
martingales = zeros(timesteps,loops);
dualityValues = zeros(timesteps,loops);

for i = 1:timesteps
    time = i*timestep;
    
    if (i==timesteps) 
        valueFunctions(i,:) = h(i,:);
    else
        basisFunctions = generateBasisFunctions(samplePaths(i+1,:),numberBasisFunctions);
        continuationValues(i,:) = (regressionCoefficients(:,i)')*basisFunctions;
        valueFunctions(i,:) = max(h(i,:),continuationValues(i,:));
    end
    
    % now we need the sub paths
    stockSubpaths = zeros(subpathLoops,loops);
    valueSubpaths = zeros(subpathLoops,loops);
    continuationSubpaths = zeros(subpathLoops,loops);
    payoffSubpaths = zeros(subpathLoops,loops);
    
    for j = 1:subpathLoops
        Z = randn(1,loops);
        stockSubpaths(j,:) = samplePaths(i,:).*exp((interest - volatility^2/2)*timestep + volatility*Z*sqrt(timestep));
        payoffSubpaths(j,:) = exp(-interest*time)*max(strike - stockSubpaths(j,:),0);
        
        if (i == timesteps)
            valueSubpaths(j,:) = payoffSubpaths(j,:);
        else
            basisFunctions = generateBasisFunctions(stockSubpaths(j,:),numberBasisFunctions);
            continuationSubpaths(j,:) = (regressionCoefficients(:,i)')*basisFunctions;
            valueSubpaths(j,:) = max(payoffSubpaths(j,:),continuationSubpaths(j,:));
        end
    end
    
    martingaleDifference = valueFunctions(i,:) - sum(valueSubpaths,1)/subpathLoops;
    if (i==1)
        martingales(i,:) = martingaleDifference;
    else
        martingales(i,:) = martingales(i-1,:) + martingaleDifference;
    end
    
    dualityValues = h(i,:) - martingales(i,:);
    
end

dualityVector = max(dualityValues,[],1);

optionUpperBound = sum(dualityVector)/loops
    































