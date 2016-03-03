% Script to run a few LSM algorithms at once
clear all,clc;

d = [50,100,250,500,1000];
lengthd = length(d);
trueValue = 20.6*ones(1,lengthd);

europeanValues = zeros(1,lengthd);
simValues = zeros(1,lengthd);
lowerBounds = zeros(1,lengthd);
upperBounds = zeros(1,lengthd);

for i = 1:lengthd
    %[controlApproximation,europeanValue,controlStdError,relativeStdError,beta,totaltime] = gameOptionLSM(100,0.06,0.5,0.4,80,10^5,10^5,d(i),4,5);
    %[europeanValue,lowerBound,lowerStdError,lowerRelativeStdError,upperBound,upperStdError,upperRelativeStdError,totaltime] = gameOptionBounds(90,d(i));
    %[europeanValue,martApproximation,martStdError,martRelativeStdError,martApproximation2,martStdError2,martRelativeStdError2,totaltime] = gameOptionMartingale(80,d(i));
    [europeanValue,lowerBound,lowerStdError,lowerRelativeStdError,upperBound,upperStdError,upperRelativeStdError,martApproximation,totaltime] = gameOptionBoundsPath(80,d(i));
    
    europeanValues(i) = europeanValue;
    simValues(i) = martApproximation;
    lowerBounds(i) = lowerBound;
    upperBounds(i) = upperBound;
end


europeanValues
simValues
lowerBounds
upperBounds

figure(1)
hold off
scatter(d,simValues,'x','LineWidth',3);
hold on
plot(d,upperBounds);
plot(d,lowerBounds);
plot(d,trueValue,'--');
set(gca,'FontSize',20);
%title({'$S_{0} = 120$'},'Interpreter','Latex');


