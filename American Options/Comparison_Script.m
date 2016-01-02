%% Compare the various methods
% Try up to 100 timesteps and compare the results in terms of time taken
% and how they agree

S0 = 36; K = 40; r = 0.06; T = 1; sigma = 0.2;

Timings = zeros(98,3);

Results = zeros(98,3);

M = 10000; % Number of grid spacings/MC paths for LSM and FD

for N = 3:20
    tic;
    Results(N-1,1) = AmericanOptCRR(S0,K,r,T,sigma,N);
    Timings(N-1,1) = toc;
    tic;
    Results(N-1,2) = AmericanOptFD(S0,K,r,T,sigma,N,M);
    Timings(N-1,2) = toc;
    tic;
    Results(N-1,3) = AmericanOptLSM(S0,K,r,T,sigma,N,M);
    Timings(N-1,3) = toc;
end

%% Plot the results of this
subplot(2,1,1);
plot(Results);
grid
title('Option price','fontsize',14);
xlabel('Number of timesteps','fontsize',14);
ylabel('Option price','fontsize',14);
legend('CRR','FD','LSM','location','SE');

subplot(2,1,2);
plot(Timings);
grid
title('Time taken to compute','fontsize',14);
xlabel('Number of timesteps','fontsize',14);
ylabel('Method timings','fontsize',14);
legend('CRR','FD','LSM','location','NW');

   
%%  How do they compare over a surface?
[Price,Pcrr,Scrr,Tcrr] = AmericanOptCRR(S0,K,r,T,sigma,100);
[Price,Pfd,Sfd,Tfd] = AmericanOptFD(S0,K,r,T,sigma,100,100);
[Price,Plsm,Slsm,Tlsm] = AmericanOptLSM(S0,K,r,T,sigma,100,100);

figure;
surf(Sfd,Tfd,Pfd); shading interp
line(Slsm,Tlsm,Plsm,'linestyle','none','marker','.','color','r');
line(Scrr,Tcrr,Pcrr,'linestyle','none','marker','.','color','b');

data.CRR.P = Pcrr;
data.CRR.S = Scrr;
data.CRR.T = Tcrr;

data.FD.P = Pfd;
data.FD.S = Sfd;
data.FD.T = Tfd;

data.LSM.P = Plsm;
data.LSM.S = Slsm;
data.LSM.T = Tlsm;

save AMERICAN_OPTION_DATA data;
