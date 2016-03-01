function basisFunctions = generateBasisFunctions(S,M);

% M is the number of basis functions, S is input vector.

basisFunctions = zeros(length(S),M);

% Must make a choice of basis functions

% %% Laguerre
for i = 1:M
    basisFunctions(:,i) = Laguerre(i-1,S);
end

% %% Weighted Laguerre
% for i = 1:M
%    basisFunctions(:,i) = exp(-S/2).*Laguerre(i-1,S);
% end
% 
%% Hermite
% for i = 1:M
%     basisFunctions(:,i) = Hermite(i-1,S);
% end



end