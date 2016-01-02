function [beta,lowerBound,lowerBoundStdError] = LSMregressioncoefficientsOld(K,r,T,s,S0,N,d,M)
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
C = zeros(d,2*N); % cash flow matrix (does not include time 0)
beta = zeros(d,M); % matrix to store the regression coefficients includes time zero but not time d
% We will work backwards from maturity. First at time d
C(d,:) = h(d,:);

for i = d-1:-1:1
    % At time i find where the stock is in the money, ie when below K
    I = S(i+1,:) < K;
    Iindex = find(S(i+1,:)<K);
    X = S(i+1,I)'; % Same X as in LS paper
    
    % Now we want to work out Y. Lets get the relevant part of the cash
    % flow matricex. ie from time i+1 onwards and only for paths in the
    % money
    Csub = C(i+1:d,I);
    [sel c]= max(Csub~=0, [], 1 );
    % sel is indicator of whether there is a non zero entry in the column and c contains
    % the index of the first non zero row (it puts 1 if there is no non
    % zero entries as well, but this does not matter since it will still extract zero entries)
    % We can extract the cash flow values
    idx = sub2ind(size(Csub),c,[1:length(X)]);
    Y1 = Csub(idx); % cash flow values
    
    Y = (exp(-r*c*dt).*Y1)'; % discounted cash flow values
    %scatter(X,Y)
    
    D = generateBasisFunctions(X,M); % design matrix
    
    % SELF-CALCULATED METHOD
    
    %     phi = zeros(M,M);
    %     phi2 = zeros(M,1);
    %     for j = 1:length(Y)
    %         temp = D(j,:)'*D(j,:);
    %         phi = phi + temp;
    %
    %         temp2 = D(j,:)'*Y(j);
    %         phi2 = phi2 + temp2;
    %     end
    %     phi = phi/length(Y);
    %     phi2 = phi2/length(Y);
    %     betaSub = phi2'/phi';
    
    % AUTOMATIC MATLAB
    betaSub = D\Y; % regression coefficients
    beta(i+1,:) = betaSub;
    
    % Work out the continuation
    % values for each of the paths in the money and compare them to the
    % exercise value.
    contValue = D*betaSub;
    exerciseValue = K - X;
    Icont = find(contValue <=  exerciseValue)';
    
    % Update the cash flow matrix on the paths that exercise is indicated
    if length(Icont) > 0
        for j = 1:length(Icont)
            index = Icont(j);
            newIndex = Iindex(index);
            newCashFlow = exerciseValue(index);
            C(i,newIndex) = newCashFlow;
            C(i+1:d,newIndex) = 0;
            
        end
    end
end

% Now we have C updated for all times > 0. Finally we just need to
% calculate the value of the option by discounting each non zero entry in C
% back to 0 and averaging over all the paths.

[sel c]= max(C~=0, [], 1 );
idx = sub2ind(size(C),c,[1:length(C)]);
Y1 = C(idx); % cash flow values
Y = (exp(-r*c*dt).*Y1); % discounted cash flow values

lowerBound = mean(Y)
lowerBoundStdError = std(Y)/sqrt(2*N)    
    