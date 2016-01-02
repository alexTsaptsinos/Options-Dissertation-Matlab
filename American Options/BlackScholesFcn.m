function [C,P] = BlackScholesFcn(A,E,r,T,s)

% A - Asset price
% E - strike price
% r - risk free interest rate
% T - maturity
% s - volatility

d1 = (log(A./E)+(r-s.^2/2).*T)./(s.*sqrt(T));
d2 = d1-s.*sqrt(T);

C = A.*normcdf(d1)-E*exp(-r*T)*normcdf(d2);
P = E*exp(-r*T)*normcdf(-d2)-A.*normcdf(-d1);