function [Price,P,S,Time] = AmericanOptCRR(S0,K,r,T,sigma,N,type)

%AmericanOptCRR - Price an american option via Cox-Ross-Rubenstein tree
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
%   N       Number of points in time grid to use (minimum is 2, default is 50)
%   type    True (default) for a put, false for a call

if nargin < 6 || isempty(N), N = 50; end
if nargin < 7, type = true; end

dt = T/N;

u = exp(sigma*sqrt(dt)); d = 1/u;                   
a = exp(r*dt); p = (a-d)/(u-d);           

% Create final Returns on the tree
S{N+1} = S0*u^N*d.^(0:2:2*N);
if type
    % Put option
    P{N+1} = max(K-S{N+1},0);
else
    P{N+1} = max(S{N+1}-K,0);
end
Time{N+1} = T*ones(1,N+1);
% Now move back through time and calculate the expected return at previous
% nodes on the tree. Compare this with the immediate return. Exercise the
% option if the immediate return is greater than the expected return

for ii = N:-1:1
    Q = zeros(1,ii);
    V = zeros(1,ii);
    for jj = 1:ii
        % Share price at current node
        V(jj) = S0*u^(ii-1)*d^(2*(jj-1));   
        % Expected value of option due if we continue to hold
        E = p*P{ii+1}(jj)/a+(1-p)*P{ii+1}(jj+1)/a;
        % Value of early exercise
        if type
            % Put option
            I = max(K-V(jj),0);
        else
            I = max(V(jj-K),0);
        end
        % Value of option at this Node
        Q(jj) = max(E,I);
    end
    S{ii} = V;
    P{ii} = Q;
    Time{ii} = ii*dt*ones(size(S{ii}));
end

Price = P{1};
P = [P{:}];
S = [S{:}];
Time = [Time{:}];