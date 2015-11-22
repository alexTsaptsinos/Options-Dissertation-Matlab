function Y = Laguerre(S,k);

if (k==0)
    Y = 1;
elseif (k==1)
    Y = 1-S;
else
    Y = (1/k)*((2*k-1-S).*Laguerre(S,k-1)-(k-1)*Laguerre(S,k-2));
end

end