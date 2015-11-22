function basisFunctions = generateBasisFunctions(S,J);

% N is the number of basis functions, S is input matrix
% Use weighted Laguerre polynomials
basisFunctions = zeros(J,length(S));

for i = 1:J
    basisFunctions(i,:) = exp(-S./2).*Laguerre(S,i-1);
end
end