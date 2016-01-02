function Y = choicePolynomials(k,S,K,maturity,r,s);
% With each entry of S, will calculate the kth Laguerre polynomial. K is
% the strike price and time is the current time.

if (k==0)
    Y = K;
elseif (k==1)
    Y = S;
elseif (k==2)
    % equivalent European option price for put option expiring at time T
   Y = BSput(K,maturity,r,s,S);
else
   Y = S.*BSput(K,maturity,r,s,S);
end

end