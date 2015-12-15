function basisFunctions = generateBasisFunctions(S,J);

% J is the number of basis functions, S is input vector
% Use weighted Laguerre polynomials
basisFunctions = zeros(J,length(S));

for i = 1:J
    basisFunctions(i,:) = exp(-S./2).*Laguerre(S,i-1);
end
end