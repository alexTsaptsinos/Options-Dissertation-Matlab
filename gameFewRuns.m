% Script to run a few LSM algorithms at once
clear all,clc;

%S0 = [80,90,100,110,120];
S0 = 120;
lengthS = length(S0);
%trueValues = [20.6,12.4,5.00,3.64,2.54];
trueValues = 2.54;

europeanValues = zeros(1,lengthS);
estimatedValues = zeros(1,lengthS);
stdErrors = zeros(1,lengthS);
relativeStdErrors = zeros(1,lengthS);
times = zeros(1,lengthS);

for i = 1:lengthS
    [controlApproximation,europeanValue,controlStdError,relativeStdError,beta,totaltime] = gameOptionLSM(100,0.06,0.5,0.4,S0(i),10^5,10^5,642,4,5);
    
    europeanValues(i) = europeanValue;
    estimatedValues(i) = controlApproximation;
    stdErrors(i) = controlStdError;
    relativeStdErrors(i) = relativeStdError;
    times(i) = totaltime;
end

relativeErrors = abs((trueValues - estimatedValues)./trueValues).*100;


europeanValues
trueValues
estimatedValues
relativeErrors
stdErrors
relativeStdErrors
times


