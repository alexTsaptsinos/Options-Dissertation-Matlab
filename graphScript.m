% Script to produce the graphs
clf
estimatedValues80 = [20.8299, 20.7433, 20.6716, 20.5977, 20.5456];
estimatedValues90 = [12.8824,12.6971,12.4805,12.4250,12.3393];
estimatedValues100 = [6.3442,5.9596,5.6340,5.4480,5.3222];
estimatedValues110 = [4.0715,3.9513,3.8397,3.7937,3.7501];
estimatedValues120 = [2.8279,2.7487,2.6692,2.6354,2.6080];

estimatedValues = [estimatedValues80;estimatedValues90;estimatedValues100;estimatedValues110;estimatedValues120];
d = [50,100,250,500,1000];
trueValues = [20.6,12.4,5,3.64,2.54]';
trueValues = repmat(trueValues,1,length(d));

figure(1)
hold all

for i=1:length(d)
    relativeErrors = abs((trueValues(i,:) - estimatedValues(i,:))./trueValues(i,:)).*100;
    plot(d,relativeErrors);
    
end

h = legend({'$S_{0}=80$','$S_{0}=90$','$S_{0}=100$','$S_{0}=110$','$S_{0}=120$'},'FontSize',20);
set(h,'Interpreter','latex');
legend('boxoff');
set(gca,'FontSize',20);

% i=5;
% figure(1)
% hold off
% plot(d,estimatedValues(i,:));
% hold on
% plot(d,trueValues(i,:),'--');
% set(gca,'FontSize',20);