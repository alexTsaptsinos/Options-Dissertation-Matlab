function [D,a,c,t,dt,s,dS] = SetupFDMatrix(S0,K,r,T,sigma,N,M)

% SetupFDMatrix - Create a matrix which can be used to solve the
% discretised version of the Black-Scholes equation
%
% Inputs:
%   S0 - Initial asset price
%   K - Strike price
%   r - interest rate
%   T - maturity
%   sigma - volatility
%   N - number of steps in time grid
%   M - Number of steps in share price grid
%
% Ref: Chapter 18, Options, Futures and Other Derivatives, John Hull, 
% Prentice-Hall, (5th edition 2003, 744 pages)

t = linspace(0,T,N+1); dt = T/N;
Smax = 2*exp(r*T)*max(S0,K); % Maximum price considered
dS = Smax/(M); 
s = 0:dS:Smax;

J = (1:M-1)';
a = r/2*J*dt-1/2*sigma^2*J.^2*dt;
b = 1+sigma^2*J.^2*dt+r*dt;
c = -r/2*J*dt-1/2*sigma^2*J.^2*dt;

D = spdiags([[a(2:end);0] b [0;c(1:end-1)]],[-1 0 1],M-1,M-1);
