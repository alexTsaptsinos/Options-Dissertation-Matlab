%% Price via Finite differences
% Solve the Black Scholes PDE by a finite difference approach. We take the
% Black Scholes PDE:
%
% $$ \frac{\partial V}{\partial t}+\frac{1}{2}\sigma^2s^2\frac{\partial^2 V}{\partial s^2}+rs\frac{\partial V}{\partial s}-rV = 0 $$
%
% And discretise this to a grid on N point in time and M point in asset
% price:
%
% $$ V_{i+1,j}-V_{i,j}+rj\delta t \frac{V_{i,j+1}-V{i,j}}{2}+\frac{1}{2}\sigma^2j^2\delta t (V_{i,j+1}+V_{i,j-1}-2V_{i,j}) = 0 $$
%
% $$ i \leq N, j \leq M $$
%
% This can be rewritten as 
%
% $$ a_jV_{i,j-1}+b_jV_{i,j}+c_jV_{i,j+1} = V_{i+1,j} $$
%
% $$ a_j = \frac{1}{2}\delta t(rj-\sigma^2 j^2) $$
%
% $$ b_j = 1+\sigma^2 j^2 \delta t +r \delta t $$
%
% $$ c_j = -\frac{1}{2}\delta t(rj+\sigma^2 j^2) $$
%
% This becomes:
%
% $$ V_{i,j} = D\backslash V^*_{i+1,j} $$
%
% for D a matrix constructed from the _a_, _b_ and _c_ values, and V* the
% next row of the V matrix augmented with boundary conditions

%% Parameters
% Set up parameters for the option. American option price will depend on
% the current price of the underlying asset _S0_, the strike price _K_, the
% risk free interest rate _r_, the time to maturity of the option _T_ and
% the volatility _sigma_

S0 = 36;    % Initial price
K = 40;     % Strike price
r = 0.06;   % Interest rate
sigma = 0.2;% Volatility
T = 1;      % Maturity (years)

%% Finite difference parameters
% In addition to the above, to price the option via finite differences, we 
% need to specify some extra parameters. These are the maximum share price
% that we want to consider in the algorithm (this should be taken to be
% large enough so that the option has effectively 0 value at this price
% throughout its life), and the grid size - the grid being taken in time
% and share price. For the purposes of this file we take 
%
% $$ S_{max} = 2e^{rT}max(S0,K) $$
%
% Twice the interest adjusted of the maximum of the initial price and the
% strike price at maturity. Along with the grid size we could vary this
% parameter to investigate the effect this choice has on the final value of
% the share price.

%%
% Given the parameters above now construct the matrices required for the
% finite diffference operation. Here we form _D_ the matrix we need to
% "divide" by as a sparse tridiagonal matrix

N = 50; % Number of timesteps
M = 50; % number of Share prices to look at

[D,a,c,t,dt,s,dS] = SetupFDMatrix(S0,K,r,T,sigma,N,M);

spy(D); shg

%% Boundary conditions
% We now need to set up a matrix for finite differences so that it has the
% boundary conditions that we want to apply in place. These are:
%
% $$ P(t,S_{max}) = 0 $$
%
% $$ P(T,S) = max(K-S,0) $$
%
% $$ P(t,0) = K $$
%

P = NaN*ones(N+1,M+1);          % Pricing Matrix (t,S)
P(end,:) = max(K-(0:M)*dS,0);   % Value of option at maturity
P(:,1) = K;                     % Value of option when stock price is 0
P(:,end) = 0;                   % Value of option when S = Smax

%% Plot initial conditions
% Create a plot to visualise what the current state of the PDE simulation 
% is.

[handles,kdx,idx,Z] = CreateFDPlot(S0,K,T,s,t,P);

%% Finite difference solver
% We now solve the finite difference problem. This is done in a for loop
% which starts and time _t = T_, and steps back through time to reach 
% _t = 0_, from where we can read off the option price. 

for ii = N:-1:1
    y = P(ii+1,2:end-1)'+[-a(1)*K; zeros(M-3,1); -c(end)*0];
    x = D\y; % Solve FD-equation
    P(ii,2:end-1) = max(x,K-s(2:end-1)'); % Take into account early exercise option
    
    % Update plot.
    set(handles(1),'ydata',P(ii,kdx));
    set(handles(2),'ydata',Z*P(ii,idx:idx+1)');
    set(handles([3,5]),'string',['t = ',num2str(T-(N+1-ii)*dt)]); 
    set(handles(4),'zdata',P); shg
    drawnow
    pause(.1);
end
    
%% Price the option
% Extract the final price and plot the resulting surface, and plot a line
% corresponding to the initial asset price.

P_FD = Z*P(1,idx:idx+1)'; % Option price

spos = get(0,'screensize');
figure('position',[spos(3)/2-400 spos(4)/2-300 800 600]);

sf = surf(s,t,P); shading interp; light
Zdata = P(:,idx:idx+1)*Z';
Ydata = t;
Xdata = S0*ones(size(t));
line(Xdata,Ydata,Zdata,'color','r','linewidth',2);
xlabel('Share price','fontsize',16);
ylabel('Time','fontsize',16); 
zlabel('Option Price','fontsize',16);
title('Price of American Put Option by Finite Differences','fontsize',16);
rotate3d on;

%% Compare with CRR method
% The Cox Ross Rubenstein tree based method also gives us prices for the
% option over time and share price, we can now ask how this compares with
% the results from the finite difference method.

load CRRData;
hold on;

crr = plot3(Scrr,Tcrr,Pcrr,'.'); shg
title('Comparison of pricing methods','fontsize',16);
legend([sf,crr],'Finite differences','CRR tree');

%% Reference:
%
% Chapter 18, _Options, Futures and Other Derivatives_, John Hull, *Prentice-Hall*, (5th edition 2003, 744 pages)