%% SINGLE ASSET AMERICAN PUT
% We use methods from Rogers to find the value of an American put option
% with a single underlying asset using the martingale of the Doob-Meyer
% decomposition of the Snell envelope.
tic
%% Set up variables
strike = 100;
interest = 0.06;
maturity = 0.5;
volatility = 0.4;
initialPrice = 120;
loops = 5000;
timesteps = 100;

%% First find European value
europeanValue = BSput(strike,maturity,interest,volatility,initialPrice)

%% Sample small number of paths for optimisation
optimisationTimeSteps = 40;
optimisationTimeStep = maturity/optimisationTimeSteps;
optimisationLoops = 300;
lambdas = 0.9:0.001:1.1;
lambdaloops = length(lambdas)
expectations = [];
for p = 1:lambdaloops
    testLambda = lambdas(p);
    supremums = [];
    for M = 1:optimisationLoops
        Snew = initialPrice;
        time = 0;
        supremum = 0;
        
        for k = 0:optimisationTimeSteps;
            % Z_t is the discounted exercise value of the option
            Z_t = exp((-interest)*time)*max(strike - Snew,0);
            
            % M_t is the discounted values of the corresponding European
            % put only calculate if in the money
            tempdiff = strike - Snew;
            if tempdiff > 0
                if maturity - time < 0.0000001
                     BSprice = max(strike - Snew,0);
                 else
                     [BSprice] = BSput(strike,maturity - time,interest,volatility,Snew);
                end
                 M_t = exp((-interest)*time)*BSprice;
            else
                M_t = 0;
            end
            
            difference = Z_t - testLambda*(M_t - europeanValue);
            
            if difference > supremum
                supremum = difference;
            end
            time = time + optimisationTimeStep;
            Z = randn;
            %Snew is the new price of the underlying asset for the next
            %loop
            Snew = Snew*exp((interest - volatility^2/2)*optimisationTimeStep + volatility*Z*sqrt(optimisationTimeStep));
        end
        
        supremums = [supremums supremum];
    end
    newExpectation = mean(supremums);
    expectations = [expectations newExpectation];
    p
end

% Now we have all the expectations for the different values of lambda,
% let's find the smallest
[M,I] = min(expectations);
lambda = lambdas(I)


%% Finished optimisation, now find value
results = [];
standardErrors = [];
for t = 1:2
    % We do it twice with different timesteps so we can apply richardson
    % extrapolation
    if t == 1
        timesteps = timesteps;
    else
        timesteps = timesteps/2;
    end
    timestep = maturity/timesteps;
    
    
    
    supremums = [];
    for i = 1:loops
        Snew = initialPrice;
        time = 0;
        supremum = 0;
        startedEuropean = 0;
        
        for j = 0:timesteps
            
            % Z_t is the discounted exercise value of the option
            Z_t = exp((-interest)*time)*max(strike - Snew,0);
            
            % M_t is the discounted values of the corresponding European
            % put only calculate if in the money
            tempdiff = strike - Snew;
            if tempdiff > 0
                if maturity - time < 0.0000001
                     BSprice = max(strike - Snew,0);
                 else
                     [BSprice] = BSput(strike,maturity - time,interest,volatility,Snew);
                end
                 M_t = exp((-interest)*time)*BSprice;
            else
                M_t = 0;
            end
                 
                 
%             if startedEuropean == 0
%                 % check to start
%                 tempdiff = strike - Snew;
%                 if tempdiff > 0
%                     startedEuropean = 1;
%                     europeanStartTime = time;
%                 end
%             end
%             
%             
%             if startedEuropean == 1
%                 % european already started
%                 if maturity - time < 0.0000001
%                     BSprice = max(strike - Snew,0);
%                 else
%                     [BSprice] = BSput(strike,maturity - time,interest,volatility,Snew);
%                 end
%                 M_t = exp((-interest)*time)*BSprice;
%             else
%                 M_t = 0;
%             end
 
            difference = Z_t - lambda*M_t;
            
            if difference > supremum
                supremum = difference;
            end
            %Svalues = [Svalues Snew];
            time = time + timestep;
            Z = randn;
            %Snew is the new price of the underlying asset for the next
            %loop
            Snew = Snew*exp((interest - volatility^2/2)*timestep + volatility*Z*sqrt(timestep));
            
        end
        supremums = [supremums supremum];
    end
    
    averageValue = mean(supremums)+ lambda*europeanValue;
    stdeviation = std(supremums);
    standarderror = stdeviation/sqrt(loops);
    standardErrors = [standardErrors standarderror];
    results = [results averageValue];
    
end


% Richardson extrapolate finally
finalValue = (4*results(2) - results(1))/3
standardErrorFinal = (4*standardErrors(2) - standardErrors(1))/3
timeFinal = toc