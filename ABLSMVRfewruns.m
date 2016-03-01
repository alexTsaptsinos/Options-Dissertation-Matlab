% Script to run a few LSM algorithms at once
clear all,clc;

S0 = [36,38,40,42,44];
lengthS = length(S0);
FDMvalues = [4.478,3.250,2.314,1.617,1.110];

europeanValues = zeros(1,lengthS);
lowerBounds = zeros(1,lengthS);
upperBounds = zeros(1,lengthS);
upperStdErrors = zeros(1,lengthS);
upperRelativeStdErrors = zeros(1,lengthS);
times = zeros(1,lengthS);
confidenceIntervals = zeros(2,lengthS);

for i = 1:lengthS
    % Choose which approx using for continuation value
    [controlLowerBound,europeanValue,upperBound,upperStdError,upperRelativeStdError,upperTime,CI] = singleAmericanAB3(S0(i));
   % [controlLowerBound,europeanValue,upperBound,upperStdError,upperRelativeStdError,upperTime,CI] = singleAmericanAB4(S0(i));
    
    europeanValues(i) = europeanValue;
    lowerBounds(i) = controlLowerBound;
    upperBounds(i) = upperBound;
    upperStdErrors(i) = upperStdError;
    upperRelativeStdErrors(i) = upperRelativeStdError;
    times(i) = upperTime;
    confidenceIntervals(:,i) = CI;
end

relativeErrors = abs((FDMvalues - upperBounds)./FDMvalues).*100;


europeanValues
FDMvalues
lowerBounds
upperBounds
relativeErrors
upperStdErrors
upperRelativeStdErrors
confidenceIntervals
times


