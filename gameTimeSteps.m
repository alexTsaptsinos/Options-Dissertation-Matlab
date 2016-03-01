% Script to run a few LSM algorithms at once
clear all,clc;

d = [50,100,250,500,1000];
lengthd = length(d);
trueValue = 12.4*ones(1,lengthd);

europeanValues = zeros(1,lengthd);
lowerBounds = zeros(1,lengthd);
upperBounds = zeros(1,lengthd);

for i = 1:lengthd
    %[controlApproximation,europeanValue,controlStdError,relativeStdError,beta,totaltime] = gameOptionLSM(100,0.06,0.5,0.4,80,10^5,10^5,d(i),4,5);
    [europeanValue,lowerBound,lowerStdError,lowerRelativeStdError,upperBound,upperStdError,upperRelativeStdError,totaltime] = gameOptionBounds(90,d(i));
    
    europeanValues(i) = europeanValue;
    lowerBounds(i) = lowerBound;
    upperBounds(i) = upperBound;
end


europeanValues
lowerBounds
upperBounds

figure(1)
hold off
plot(d,lowerBounds);
hold on
plot(d,upperBounds);
plot(d,trueValue,'--');
set(gca,'FontSize',20);
%title({'$S_{0} = 120$'},'Interpreter','Latex');


