% Script to run a few LSM algorithms at once
clear all,clc;

S0 = [36,38,40,42,44];
lengthS = length(S0);
FDMvalues = [4.478,3.250,2.314,1.617,1.110];

europeanValues = zeros(1,lengthS);
simulatedValues = zeros(1,lengthS);
stdErrors = zeros(1,lengthS);
times = zeros(1,lengthS);
relativeStdErros = zeros(1,lengthS);

for i = 1:lengthS
    % Choose whether standard or variance reduction method
    %[LSMlowerbound,europeanValue,LSMstdError,totaltime,relativeStdError] = singleAmericanLongstaffSchwartz(S0(i));
    [LSMlowerbound,europeanValue,LSMstdError,totaltime,relativeStdError,beta] = singleAmericanLSMAntithetic(S0(i));
    
    europeanValues(i) = europeanValue;
    simulatedValues(i) = LSMlowerbound;
    stdErrors(i) = LSMstdError;
    times(i) = totaltime;
    relativeStdErrors(i) = relativeStdError;
end

relativeErrors = abs((FDMvalues - simulatedValues)./FDMvalues).*100;


europeanValues
FDMvalues
simulatedValues
stdErrors
relativeErrors
relativeStdErrors
times


