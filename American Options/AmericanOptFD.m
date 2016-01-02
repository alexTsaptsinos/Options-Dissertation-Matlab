function [P_FD,P,s,t] = AmericanOptFD(S0,K,r,T,sigma,N,M,type)

%AmericanOptFD - Price an american option via finite differences
%
%  Returns the price of an American option computed using finite
%  difference method applied to the BlackScoles PDE. 
%
% Inputs:
%
%   S0      Initial asset price
%   K       Strike Price
%   r       Interest rate
%   T       Time to maturity of option
%   sigma   Volatility of underlying asset
%   N       Number of points in time grid to use (minimum is 3, default is 50)
%   M       Number of points in asset price grid to use (minimum is 3, default is 50)
%   type    True (default) for a put, false for a call

if nargin < 6 || isempty(N), N = 50; elseif N < 3, error('N has to be at least 3'); end
if nargin < 7 || isempty(M), M = 50; elseif M < 3, error('M has to be at least 3'); end
if nargin < 8, type = true; end

% create time grid

t = linspace(0,T,N+1);
dt = T/N; % Time step

% Share price grid
Smax = 2*max(S0,K)*exp(r*T); % Maximum price considered
dS = Smax/(M);
s = 0:dS:Smax;

% Now find points either side of the initial price so that we can calculate
% the price of the option via interploation
idx = find(s < S0); idx = idx(end); a = S0-s(idx); b = s(idx+1)-S0;
Z = 1/(a+b)*[a b]; % Interpolation vector

% Set up a pricing matrix to hold the values we compute
P = NaN*ones(N+1,M+1); % Pricing Matrix (t,S)

% Boundary condition
if type
    P(end,:) = max(K-(0:M)*dS,0); % Value of option at maturity - Put
else
    P(end,:) = max((0:M)*dS-K,0); % Value of option at maturity - Call
end
P(:,1) = K; % Value of option when stock price is 0)
P(:,end) = 0; % Value of option when S = Smax

% Create matrix for finite difference calculations
J = (1:M-1)';
a = r/2*J*dt-1/2*sigma^2*J.^2*dt;
b = 1+sigma^2*J.^2*dt+r*dt;
c = -r/2*J*dt-1/2*sigma^2*J.^2*dt;
D = spdiags([[a(2:end);0] b [0;c(1:end-1)]],[-1 0 1],M-1,M-1);

% Finite difference solver
for ii = N:-1:1
    y = P(ii+1,2:end-1)'+[-a(1)*K; zeros(M-3,1); -c(end)*0];
    x = D\y; % Value of the option
    if type
        P(ii,2:end-1) = max(x,K-s(2:end-1)'); % Put
    else
        P(ii,2:end-1) = max(x,s(2:end-1)'-K); % Call
    end
end
% Extract the final price
P_FD = Z*P(1,idx:idx+1)';

end