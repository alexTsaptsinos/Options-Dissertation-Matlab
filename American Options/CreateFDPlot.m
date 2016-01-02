function [h,kdx,idx,Z] = CreateFDPlot(S0,K,T,s,t,P)

% CreateFDPlot - Use MATLAB graphics commands to create a custom plot to
% illustrate the finite difference approach to solving the price problem

idx = find(s < S0); idx = idx(end); a = S0-s(idx); b = s(idx+1)-S0;
Z = 1/(a+b)*[b a]; % Interpolation vector
% create a Logical index for points close to the initial asset price
kdx = (s >= S0-10) & (s <= S0+10); 

% Create figure
spos = get(0,'screensize');
fpos = [spos(3)/2-300, spos(4)/2-300 600 600];
figure('position',fpos);

subplot(2,1,1);
h = plot(s(kdx),P(end,kdx),'.-',S0,Z*P(end,(idx:idx+1))','ro'); 
ax = gca;
set(ax,'ylim',[0,K-S0+10]);
h(3) = title(['t = ',num2str(T)],'fontsize',18); 
xlabel('Share price');
ylabel('Option price'); 
grid on;
shg

subplot(2,1,2);
h(4) = surf(s,t,P,'edgecolor','none'); 
shading interp;
hold on; 
line(s,t(end)*ones(size(s)),P(end,:),'color','k','linewidth',2);
line(s(1)*ones(size(t)),t,P(:,1),'color','k','linewidth',2);
line(s(end)*ones(size(t)),t,P(:,end),'color','k','linewidth',2);
xlabel('Share price'); 
ylabel('Time'); 
zlabel('Option price'); 
grid on
h(5) = title(['t = ',num2str(T)],'fontsize',18);
