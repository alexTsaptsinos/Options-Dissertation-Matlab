function basisFunctions = generateChoiceFunctions(S,M,K,maturity,r,s);

% M is the number of basis functions, S is input vector. K is strike price,
% only needed for choice polynomials and maturity is time to maturity also only
% needed for choice polynomials.
% Use weighted Laguerre polynomials
basisFunctions = zeros(length(S),M);

%% Choice Polynomials
for i = 1:M
    basisFunctions(:,i) = choicePolynomials(i-1,S,K,maturity,r,s);
end


end