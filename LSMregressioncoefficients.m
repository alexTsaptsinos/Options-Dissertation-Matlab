function [regressionCoefficients] = LSMregressioncoefficients(strike,interest,maturity,volatility,initialPrice,loops,timesteps,numberBasisFunctions);
%% LSM Method
% Function to calculate the regression coefficients using the LSM method

timestep = maturity/timesteps;

regressionCoefficients = zeros(numberBasisFunctions,timesteps);
continuationValues = zeros(timesteps,loops);

%% Generate sample paths
% First we generate all the sample paths in a matrix of size loops x (timesteps +
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

%% Calculate regression coefficients
% Vmatrix stores the value of the option for each trajectory
Vmatrix = zeros(1,loops);
% set the final option value
Vmatrix = h(timesteps,:);

% Now use backward induction to find the regression coefficients,
% continuation values and option values

for j = timesteps-1:-1:1
    
    % now need the basis functions, grab them
    basisFunctions = generateBasisFunctions(samplePaths(j+1,:),numberBasisFunctions);
    
    sum1 = 0;
    sum2 = 0;
    
    for k = 1:loops
        sum1 = sum1 + basisFunctions(:,k)*basisFunctions(:,k)';
        sum2 = sum2 + basisFunctions(:,k)*Vmatrix(k);
    end
    
    regressionCoefficients(:,j) = sum1\sum2;
    continuationValues(j,:) = (regressionCoefficients(:,j)')*basisFunctions;
    Vmatrix = max(h(j,:),continuationValues(j,:));
end
    
    