function Y = Hermite(k,S);
% With each entry of S, will calculate the kth Laguerre polynomial

if k==0
    Y = 1;
elseif k==1
    Y = 2.*S;
else
    Y = 2.*S.*Hermite(k-1,S) - 2*(k-1).*Hermite(k-2,S);
end

end
% if (k==0)
%     Y = 1;
% elseif (k==1)
%     Y = 1-S;
% else
%     Y = (1/k)*((2*k-1-S).*Laguerre(k-1,S)-(k-1)*Laguerre(k-2,S));
% end
% 
% end