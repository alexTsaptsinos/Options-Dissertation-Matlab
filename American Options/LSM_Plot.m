function LSM_Plot(S0,K,r,T,sigma,N,M)

% Function to demonstrate the Longstaff-Schwartz Least Squares algorithm
%
% Ref: Longstaff, F. A., and E.S. Schwartz, "Valuing American Options by
% Simulation: A Simple Least-Squares," The Review of Financial Studies, 14
% no. 1 (Spring 2001), pp. 113-147
% Approach,
%
% Inputs:
%
%   S0 - Current Asset price
%   K  - Strike price of option
%   r  - risk free interest rate
%   sigma - Volatility of underltying stock
%   T  - Maturity of option
%   N  - Number of timesteps to take
%   M  - Number of paths to simulate.

% Place holders for Monte Carlo simulation variables
dt = []; % Time step
R = []; % Returns
S = []; % Monte carlo paths
t = []; % Timebase

% Create place holders for variables we use in the nested functions
X = []; % Current asset values of in the money options
Y = []; % Discounted cashflows from future pay off of in the money options
x = []; % x-data for graphs
C = []; % Predicted cashflow from regression model
a = []; % Coefficients of model
handles = []; % Handles to graphical objects
CurrentTime = []; % Current time step.
Idx = []; % Indices of in the money asset paths
CF = [];
LegH = []; % Handles to objects to show in legend in third axes
LegendHandle = [-1 -1 -1]; % Handle to Legend

% Exercise time for each path - record when we exercise an option early
ExTime = (N+1)*ones(M,1);

fs = 12; % Fontsize to use in labels

% Function handles for Laguerre Polynomials
L_0 = @(x) ones(size(x));
L_1 = @(x) (1-x);
L_2 = @(x) 1/2*(2-4*x-x.^2);
% handle to function which will generate regression matrix
L = @(x) [L_0(x) exp(-x/2).*L_1(x) exp(-x/2).*L_2(x)];

% Plotting flag - if the run button is pressed we want to bypass all
% plotting steps
pFlag = true;
% Current state of the algorithm - used when stepping through the algorithm
status = 1;

% Call the creation function
nCreate;
nSetup;

% -------------------------------------------------------------------------
    function nHelp(src,evt) %#ok
        % Open a web browser pointing to the original Longstaff-Schwartz

    end

% -------------------------------------------------------------------------
    function nIdentify

        % Identify "in the money" options
        Idx = S(CurrentTime,:) < K;
        if pFlag
            % Clear Axis
            nClearAxes(2);
            for jj = 1:M
                if S(CurrentTime,jj) > K
                    % Plot out of the money asset paths as dotted blue
                    % lines
                    oom = plot(t(1:ExTime(jj)),S(1:ExTime(jj),jj),'g:',...
                        'parent',handles.Axes(2));
                else
                    % Plot in the money asset paths as black lines up to
                    % the current time, and then red lines for the rest of
                    % their life
                    tIdx = 1:CurrentTime; tJdx = CurrentTime:ExTime(jj);
                    itm = plot(t(tIdx),S(tIdx,jj),'k','parent',handles.Axes(2));
                    plot(t(tJdx),S(tJdx,jj),'r','parent',handles.Axes(2));
                end
            end
            LegendHandle(2) = legend([oom(1),itm(1)],...
                'Out of the money paths','In the money paths',...
                'location','NW');
            % Update the message box
            nMessage(['Identify asset paths which are in the money at ',...
                'current time step']);

        end
    end % end of nIdentify
% -------------------------------------------------------------------------
    function nPlot
        % Find the values of the current asset prices from in the money
        % paths and then find the discounted cash flow from holding on the
        % the option for another time step
        X = S(CurrentTime,Idx)';
        Y = CF(CurrentTime+1,Idx)'*exp(-r*dt); % Discounted cashflow
        if pFlag
            nClearAxes(3);
            LegH = [];
            % Plot the points we are going to try and fit with a regression
            % model
            LegH(1) = plot(X,Y,'rx','markersize',12,'parent',handles.Axes(3));
            LegendHandle(3) = legend(LegH,'Discounted payoff');
            nMessage(['Find the discounted Payoff from the option if we choose ',...
                'to hold it for another time step.']);
        end

    end % end of nPlot
% -------------------------------------------------------------------------
    function nRegress
        % Regress the discounted cash flow on to the current price and then
        % plot the resulting line
        R = L(X/S0); % Regression matrix
        a = R\Y; % Linear regression step
        C = R*a; % Cash flows as predicted by the model
        if pFlag
            xlim = get(handles.Axes(3),'xlim');
            x = linspace(xlim(1),xlim(2),51)';
            xM = L(x/S0);
            LegH(2) = plot(x,xM*a,'parent',handles.Axes(3));
            nMessage(['Fit a model to the discounted payoff']);
            LegendHandle(3) = legend(LegH,'Discounted payoff','Payoff model');
        end
    end % end of nRegress
% -------------------------------------------------------------------------
    function nExercise

        % Plot the pay off from early exercise of the option
        if pFlag
            LegH(3) = plot(x,max(K-x,0),'k','parent',handles.Axes(3));
            nMessage(['Plot the payoff achieved through early ',...
                'exercise of the option']);
            LegendHandle(3) = legend(LegH,'Discounted payoff','Payoff model','Early exercise payoff');
        end

    end % end of nExercise
% -------------------------------------------------------------------------
    function nModel
        % Plot the values of the predicted cash flow and the values
        % obtained by early exercise of the option.
        if pFlag
            nClearAxes(3);
            LegH = [];
            plot(x,max(K-x,0),'k','parent',handles.Axes(3));
            LegH(1) = plot(X,max(K-X,0),'ko','parent',handles.Axes(3));
            plot(x,L(x/S0)*a,'b','parent',handles.Axes(3));
            LegH(2) = plot(X,L(X/S0)*a,'bx','markersize',18,'parent',handles.Axes(3));
            nMessage(['Plot the expected payoff from holding the option and',...
                ' the pay off from early exercise for in ',...
                'the money asset paths']);
            LegendHandle(3) = legend(LegH,'Early Exercise values','Estimated payoff values');
        end

    end % End of nModel
% -------------------------------------------------------------------------
    function nDecision
        % Plot the maximum of Expected cashflow and holding on to the
        % option.
        if pFlag
            nClearAxes(3);
            LegH = [];
            plot(x,max(K-x,0),'k','parent',handles.Axes(3));
            plot(x,L(x/S0)*a,'b','parent',handles.Axes(3));
            Y1 = L(X/S0)*a;
            Y2 = max(K-X,0);
            ind = Y1 > Y2;
            LegH = plot(X(ind),Y1(ind),'g.',X(~ind),Y2(~ind),'m.',...
                'markersize',18,'parent',handles.Axes(3));
            nMessage(['Select maximum of payoff through immediate ',...
                'exercise and payoff through holding the option']);
            if numel(LegH) == 2
                strs = {'Hold Option','Exercise Option'};
            else
                strs = {'Hold Option'};
            end
            LegendHandle(3) = legend(LegH,strs);
        end
    end % end of nDecision
% -------------------------------------------------------------------------
    function nUpdate

        % Work out which options to exercise early
        Jdx = max(K-X,0) > C; % Immediate exercise better than predicted cashflow
        Idx = find(Idx); % Convert from logical to double indices
        nIdx = setdiff((1:M),Idx(Jdx)); % Find the options we are not exercising
        CF(CurrentTime,Idx(Jdx)) = max(K-X(Jdx),0); % Cash flows from exercised options
        ExTime(Idx(Jdx)) = CurrentTime; % Exercise time for the exercised options
        % Now compute cash flows from non-exercised option - take cash flow
        % at next step and discount it according to interest rate
        CF(CurrentTime,nIdx) = exp(-r*dt)*CF(CurrentTime+1,nIdx);
        oom = [];
        itm = [];
        ep = [];
        nClearAxes(2);

        for jj = 1:M
            if S(CurrentTime,jj) > K
                % Plot options which are not in the money as blue dotted lines
                oom = plot(t(1:ExTime(jj)),S(1:ExTime(jj),jj),'g:','parent',handles.Axes(2));
            else
                % Plot in the money options as black lines
                tIdx = 1:ExTime(jj);
                itm = plot(t(tIdx),S(tIdx,jj),'k','parent',handles.Axes(2));
                if ExTime(jj) == CurrentTime
                    % Exercise the option now - demonstrate these paths
                    % by plotting their end points in red
                    ep = plot(t(ExTime(jj)),S(ExTime(jj),jj),'r.',...
                        'markersize',12,'parent',handles.Axes(2));
                end
            end
        end
        LegH = [];
        strs = {};
        lcount = 1;
        if ~isempty(oom)
            LegH(lcount) = oom(1);
            strs{lcount} = 'Out of the money paths';
            lcount = lcount+1;
        end
        if ~isempty(itm)
            LegH(lcount) = itm(1);
            strs{lcount} = 'In the money paths';
            lcount = lcount+1;
        end
        if ~isempty(ep)
            LegH(lcount) = ep(1);
            strs{lcount} = 'Exercised options';
        end
        LegendHandle(2) = legend(LegH,strs,...
            'location','NW');
        nMessage(['Terminate the asset paths where ',...
            'we exercise the option.']);

        drawnow;
        pause(0.05);
        % Now move current time to the previous time step (i.e. step back
        % in time).
        CurrentTime = CurrentTime-1;
        if CurrentTime == 1
            % We are at the first time time - now more stepping possible,
            % change the controls to reflect this
            nFinish;
            set(handles.Controls(1),'enable','off');
            set(handles.Controls(2),'string','Restart','callback',@nRestart);
        end
    end % end of nUpdate
% -------------------------------------------------------------------------
    function nStep(src,evt) %#ok
        % Step through the various stages of the algorithm
        switch status
            case 1
                nIdentify;
            case 2
                nPlot;
            case 3
                nRegress;
            case 4
                nExercise;
            case 5
                nModel;
            case 6
                nDecision;
            case 7
                nUpdate;
        end
        if CurrentTime == N+1
            % If at the final time step the algorithm won't work - need
            % future cash flow, so just knock us back to the previous time
            % step and keep status as 1. We don't start at CurrentTime = N
            % since the Identify step with CurrentTime = N+1 plots all the
            % asset paths and shows which are in the money and which aren't
            CurrentTime = N;
            status = 1;
        else
            status = rem(status,7)+1;
        end
    end % end of nStep
% -------------------------------------------------------------------------
    function nRun(src,evt) %#ok
        % Run through the algorithm to conclusion
        pFlag = false; % Turn off all graphing
        nClearAxes([2,3]);
        for ii = CurrentTime:-1:2
            nIdentify;
            if CurrentTime == N+1
                CurrentTime = N;
            else
                nPlot;
                nRegress;
                nExercise;
                if CurrentTime == 2
                    pFlag = true; % We want final display
                end
                nUpdate;
            end
        end
    end % end of nRun
% -------------------------------------------------------------------------
    function nFinish(src,evt) %#ok
        % Calculate the price of the option and place it in the base
        % workspace
        P_LSM = mean(CF(2,:))*exp(-r*dt);
        assignin('base','P_LSM',P_LSM);
        str = ['Price = ',num2str(P_LSM)];
        disp(str);
        nMessage(str);

    end % end of nFinish
% -------------------------------------------------------------------------
    function nRestart(src,evt) %#ok
        % Reset the process
        nSetup;
        set(handles.Controls(1),'enable','on');
        set(handles.Controls(2),'string','Run','callback',@nRun);
    end % end of nRestart
% -------------------------------------------------------------------------
    function nCreate
        % Create the GUI
        % Get the screensize so that we can draw this GUI in the centre
        Spos = get(0,'screensize');
        col = get(0,'DefaultUicontrolBackgroundcolor');
        fHgt = 600;
        fWdh = 800;

        f = figure('numbertitle','off',...
            'name','Longstaff Schwarz process',...
            'position',[(Spos(3)-fWdh)/2, (Spos(4)-fHgt)/2 fWdh fHgt],...
            'resize','off',...
            'menubar','none',...
            'color',col);

        wdh = fWdh-10; hgt = fHgt-10;

        background = axes('units','pixels',...
            'position',[5 5 wdh hgt],...
            'xtick',[],'ytick',[],'box','off',...
            'xlim',[0 1],'ylim',[0,1],...
            'vis','off',...
            'color',col);

        nPanel([0 0 80 30]);
        nPanel([85 0 705 30]);
        nPanel([0 35 470 275]);
        nPanel([0 315 470 275]);
        nPanel([475 35 315 555]);

        %------------------------------------------------------------------
        function nPanel(pos)
            % Build etched in panels
            x = pos(1); y =pos(2); w = pos(3); h = pos(4);
            dc = .1;
            xdata = x/wdh+w/wdh*[0 0 1];
            ydata = y/hgt+h/hgt*[0 1 1];
            line(xdata,ydata,'parent',background,'color',col-dc);

            xdata = x/wdh+w/wdh*[0 1 1];
            ydata = y/hgt+h/hgt*[0 0 1];
            line(xdata,ydata,'parent',background,'color',col+dc);

        end % end of nPanel
        %------------------------------------------------------------------
        % Create axis to display the asset paths
        ax(1) = axes('units','pixels',...
            'position',[60 360 390 200]);

        % Create second axis - this is used for displaying asset paths up
        % to the point when we exercise them
        ax(2) = axes('units','pixels',...
            'position',[60 90 390 200],'nextplot','add');

        % Create a grid to demonstrate the modelling step of the algorithm
        ax(3) = axes('units','pixels',...
            'position',[525 150 245 410]);

        NumControls = 2;
        ControlW = 60;
        ControlH = 20;
        xgap = (315-NumControls*ControlW)/(NumControls+1);

        % UIcontrols
        u = uicontrol(f,'style','push',...
            'string','Step',...
            'callback',@nStep,...
            'position',[475+xgap, 45, ControlW, ControlH]);

        u(2) = uicontrol(f,'style','push',...
            'string','Run',...
            'callback',@nRun,...
            'position',[475+2*xgap+ControlW, 45, ControlW, ControlH]);
        %{
        u(3) = uicontrol(f,'style','push',...
            'string','Help',...
            'callback',@nHelp,...
            'position',[475+3*xgap+2*ControlW, 45, ControlW, ControlH]);
            %}
            u(4) = text(10/wdh,15/hgt,'','parent',background);
            u(5) = text(90/wdh,15/hgt,'','parent',background);

            handles.Controls = u;
            handles.Axes = ax;
            handles.Figure = f;

    end % end of nCreate

% -------------------------------------------------------------------------
    function nSetup

        nClearAxes([1 2 3]);
        CurrentTime = N+1; % Current time step.
        % Calculate the time step
        dt = T/N;

        % Generate the asset paths
        R = exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*randn(N,M)); % Matrix of returns
        S = cumprod([S0*ones(1,M); R]); % Asset price paths
        t = linspace(0,T,N+1); % Time base
        CF = zeros(size(S)); % Cash flow matrix
        CF(end,:) = max(K-S(end,:),0); % Option only pays off if it is in the money

        % Plot asset paths
        plot(t,S,'b','parent',handles.Axes(1)); grid on;

        ylim = get(handles.Axes(1),'ylim'); % Want second axis to have these limits
        axes(handles.Axes(1));
        title('Monte-Carlo Simulation of Asset paths','fontsize',fs);
        xl = get(handles.Axes(1),'xlabel');
        yl = get(handles.Axes(1),'ylabel');
        set(xl,'string','Time','fontsize',fs);
        set(yl,'string','Price','fontsize',fs);

        set(handles.Axes(2),'ylim',ylim);
        axes(handles.Axes(2));
        title('Asset paths up to point of exercise','fontsize',fs);
        xl = get(handles.Axes(2),'xlabel');
        yl = get(handles.Axes(2),'ylabel');
        set(xl,'string','Time','fontsize',fs);
        set(yl,'string','Price','fontsize',fs);
        
        xl = get(handles.Axes(3),'xlabel');
        yl = get(handles.Axes(3),'ylabel');
        set(xl,'string','Asset Price','fontsize',fs);
        set(yl,'string','Payoff','fontsize',fs);
        axes(handles.Axes(3));
    end % end of nSetup

% -------------------------------------------------------------------------
    function nClearAxes(I)

        for k = 1:numel(I)
            if ishandle(LegendHandle(I(k)))
                delete(LegendHandle(I(k)));
            end
            cla(handles.Axes(I(k)));
            set(handles.Axes(I(k)),'NextPlot','Add');
        end

    end % end of nClear Axes
% -------------------------------------------------------------------------
    function nMessage(str)

        % Set the string in the message box
        set(handles.Controls(4),'string',['t =  ',num2str(t(CurrentTime))]);
        set(handles.Controls(5),'string',str);

    end % end of nMessage

end % End of file