function [FDMValue] = singleAmericanPutFDM(K,r,T,sigma,S0,N,d);
% K,r,T,s,S0 same as usual in option price. d is number of timesteps
% N is number of sammple paths

% Pricing Amercian Option Using FDM Penalty Method  

Smax = 10*S0;

S = linspace(0,Smax,N+1)'; 

dS = S(2)-S(1);
dt = T/d;

lambda = 0.5*sigma^2*S.^2/dS^2;
gamma = 0.5*r*S/dS;

a = lambda - gamma;
b = - 2*lambda - r;


c = lambda + gamma;

a(end) = - 2*gamma(end);
b(end) = 2*gamma(end) - r;

L = sparse_tri(a,b,c);
I = speye(length(S));

u = max(0,K-S);
payoff=u;
cnt = zeros(1,d);
ue = zeros(1,d+1);
ue(1) = K;
rho = 10000;

for n = 1:d
    rhs = (I+0.5*dt*L)*u;

    u_old = u;
    conv = 0;
    count = 0;

    while ~conv
        P = max(0,K-S-u);
        PP = spdiags( -sign(P), 0, N+1,N+1);
        du = (I-0.5*dt*L-rho*dt*PP) \ (rhs-(I-0.5*dt*L)*u+rho*dt*P);
        u = u + du;
        conv = max(abs(du))<1e-10;
        count = count + 1;
    end
    
    cnt(n) = count;
    ue(n+1) = S(min(find(u > (K-S))));
end

FDMValue = u(1 + N*S0/Smax);


function A = sparse_tri(a,b,c)

A = spdiags([c b a],[-1 0 1],length(a),length(a))';
