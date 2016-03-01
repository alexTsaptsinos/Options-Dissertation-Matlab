function [beta] = gameOptionCoefficients(K,r,T,s,S0,N,d,M,delta)
%% LSM Method
% Function to calculate the regression coefficients using the LSM method.
% These regression coefficients can then be used to provide an exercise
% strategy to be used in Anderon/Broadie simulation.
dt = T/d;

%% Generate sample paths
% First we generate all the sample paths in a matrix of size (timesteps +
% 1) x loops, so each column corresponds to a different path
S = zeros(d+1,2*N);

% the first entry in each row will be the initial price
S(1,:) = S0;

for i = 2:d+1;
    Z = randn(1,N);
    Z = [Z,-Z]; % create antithetic pairs
    S(i,:) = S(i-1,:).*exp((r - s^2/2)*dt + s*Z*sqrt(dt));
end

%% Calculate the payoff matrix
% h is a matrix of size timesteps x loops, so each column corresponds to
% the payoffs along a path at each time NOT DISCOUNTED BACK TO TIME 0. Note
% that time 0 is not included so when matching with S there will be one
% rows difference.
h = max(K-S(2:d+1,:),0);

%% Calculate regression coefficients
C = zeros(1,2*N); % continuation values for each path
beta = zeros(d,M); % matrix to store the regression coefficients includes time zero but not time d
% We will work backwards from maturity. First at time d
C(1,:) = exp(-r*dt).*h(d,:); % first set the continuation value as the final time exercise value

for i = d-1:-1:1
    % METHOD 1
    % At time i find where the stock is in the money, ie when below K
    %S(i+1,:);
%         I = S(i+1,:) < K;
%         Iindex = find(S(i+1,:)<K);
%         X = S(i+1,I)'; % Same X as in LS paper
%     
%         % Now we want to work out Y. Lets get the relevant continuation values
%         % ie only for paths in the money
%         Csub = C(1,I);
%     
%         Y = Csub';
%         %scatter(X,Y)
%     
%         % CHOOSE either basis functions or choice ones
%         %D = generateBasisFunctions(X,M); % design matrix
%         D = generateChoiceFunctions(X,M,K,(d-i)*dt,r,s);
%     
%     
%         betaSub = D\Y; % regression coefficients
%         beta(i+1,:) = betaSub;
%     
%         % Work out the continuation
%         % values for each of the paths in the money and compare them to the
%         % exercise value.
%         contValue = D*betaSub;
%         exerciseValue = max(K - X,0);
%         penaltyValue = exerciseValue + delta;
%         Iex = find(exerciseValue > 0 & contValue <=  exerciseValue)';
%         Ipen = find(contValue >= penaltyValue)';
%     
%         C(1,Iindex(Iex)) = exerciseValue(Iex);
%         C(1,Iindex(Ipen)) = penaltyValue(Ipen);
%     
%         % Continuation values are discounted one more time step.
%         C(1,:) = exp(-r*dt).*C(1,:);
    
    % METHOD 2
%     I = S(i+1,:) < K;
%     Iindex = find(S(i+1,:)<K);
%     X = S(i+1,I)';
%     Csub = C(1,I);
%     Y = Csub';
%     D = generateBasisFunctions(X,M); % design matrix
%    % D = generateChoiceFunctions(X,M,K,(d-i)*dt,r,s);
% 
%     betaSub = D\Y;
%     beta(i+1,:) = betaSub;
%     
%     contValue = D*betaSub;
%     exerciseValue = max(K-X,0);
%     penaltyValue = exerciseValue + delta;
%     C(1,I) = min(penaltyValue,max(exerciseValue,contValue));
%     C(1,:) = exp(-r*dt).*C(1,:);
    
    
    % METHOD 3
    X = S(i+1,:)';
    Y = C(1,:)';
    %D = generateBasisFunctions(X,M); % design matrix
    D = generateChoiceFunctions(X,M,K,(d-i)*dt,r,s);

    
    betaSub = D\Y;
    beta(i+1,:) = betaSub;
    
    contValue = D*betaSub;
    exerciseValue = max(K-X,0);
    penaltyValue = exerciseValue + delta;
    C(1,:) = exp(-r*dt).*min(penaltyValue,max(exerciseValue,contValue));

    
end

% The continuation values have been computed now all the way to time zero.
% All that remains is to compute the value of the option.

approximation = mean(C(1,:))
approximationStdError = std(C(1,:))/sqrt(2*N);
end
