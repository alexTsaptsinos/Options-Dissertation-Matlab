function PlotCRRTree(S0,K,r,T,sigma,N,flag)

% PlotCRRTree the values in the CRR tree give by P.
%   This function plots a tree that contains all the points in the Cox-Ross
%   Rubenstein tree given by by the input P. The nodes in the tree are
%   coloured according to the value of the derivative at that point in the
%   tree.
%
% Inputs:
%
%   S0      Initial asset price
%   K       Strike Price
%   r       Interest rate
%   T       Time to maturity of option
%   sigma   Volatility of underlying asset
%   N       Number of timesteps to take
%
%   flag - true (default) or false, governs whether we add text markers to
%   the tree to indicate the value of the options at each node of the tree

if nargin < 7, flag = true; end % Set up default value for flag

LB = []; UB = []; Values = [];
Controls = [];
ResultsText = [];

nCreate;
Price = [];
P = [];
S = [];
Time = [];
dt = [];
u = [];
d = [];
p = [];
a = [];
N = 8;

count = 1;

L = []; T = []; A = [];

nPlotTree;

% -------------------------------------------------------------------------
    function nPlotTree
        
        S0 = Values(1); K = Values(2); r = Values(3); T = Values(4); sigma = Values(5);
        [Price,P,S,Time,dt,u,d,p,a] = iCalcTree(S0,K,r,T,sigma,N);
        count = 1;
        L = zeros(N*(N+1)/2,1); 
        T = L; 
        A = zeros(N*(N-1)/2,1);% Place holders for the handles we create.

        for ii = 1:N
            for jj = 1:numel(P{ii})
                % create the node
                L(count) = line(ii/(N+1),(ii-2*(jj-1)-1)/(N+1),'color','k','marker','o');
                if flag % Create a text label if requested
                    status = 'on';
                else
                    status = 'off';
                end
                T(count) = text((ii+.1)/(N+1),(ii-2*(jj-1)-1.1)/(N+1),...
                    ['£',sprintf('%5.3f',P{ii}(jj))],'visible',status);
                if ii < N % Now draw the lines in the tree
                    set(L(count),'buttondownfcn',{@nButtonDown,ii,jj,count});
                    xdata = [ii+1;ii;ii+1]/(N+1);
                    ydata = (ii-2*(jj-1)-1+[1;0;-1])/(N+1);
                    ExpectedValue = 1/a*(p*P{ii+1}(jj)+(1-p)*P{ii+1}(jj+1));
                    EarlyExercise = max(K-S{ii}(jj),0);
                    if EarlyExercise <= ExpectedValue
                        A(count) = line(xdata,ydata,...
                            'color','b',...
                            'linestyle','-',...
                            'hittest','off');
                    else
                        A(count) = line(xdata,ydata,...
                            'color','r',...
                            'linestyle',':',...
                            'hittest','off');
                    end

                end
                count = count+1;
            end
        end


        str = ['P(up): ',sprintf('%5.3f',p),char(10),...
            'Up factor: ', sprintf('%5.3f',u),char(10),...
            'Down factor: ', sprintf('%5.3f',d)];

        ResultsText = text(0.01, .99,'','fontsize',14,'VerticalAlignment','Top');
        ParamText = text(0.01,-0.99,'','fontsize',14,'VerticalAlignment','Bottom');
        set(ParamText,'string',str);
        Pricestr = ['Option Price: £',sprintf('%5.3f',Price)];
        ttl = get(gca,'title');
        set(ttl,'string',Pricestr,'fontsize',16);

    end
% -------------------------------------------------------------------------
    function nButtonDown(src,evt,ii,jj,Idx) %#ok

        set(A,'linewidth',.5);
        set(A(Idx),'linewidth',5);

        ExpectedValue = 1/a*(p*P{ii+1}(jj)+(1-p)*P{ii+1}(jj+1));
        EarlyExercise = max(K-S{ii}(jj),0);
        if EarlyExercise <= ExpectedValue
            Decision = 'Hold Option';
            dCol = 'b';
        else
            Decision = 'Exercise Option';
            dCol = 'r';
        end

        str = [...
            'Share price: ',sprintf('%5.3f',S{ii}(jj)),char(10),...
            'Hold Value: ', sprintf('%5.3f',ExpectedValue),char(10),...
            'Exercise Value: ',sprintf('%5.3f',EarlyExercise),char(10),...
            Decision];
      
        set(ResultsText,'string',str,'color',dCol);

    end
% -------------------------------------------------------------------------
    function nShowPrice(src,evt) %#ok
        set(DisplayPanel,'str',Pricestr);
    end
% -------------------------------------------------------------------------
    function nEarlyExercise(src,evt) %#ok

        val = get(src,'checked');
        if strcmp(val,'off')
            set(src,'checked','on');
            count = 1;
            for kk = 1:N-1
                for ll = 1:numel(P{kk})
                    ExpectedValue = 1/a*(p*P{kk+1}(ll)+(1-p)*P{kk+1}(ll+1));
                    EarlyExercise = max(K-S{kk}(ll),0);
                    if EarlyExercise <= ExpectedValue
                        set(A(count),'linestyle','-');
                    else
                        set(A(count),'linestyle',':');
                    end
                    count = count+1;
                end
            end
        else
            set(src,'checked','off');
            set(A,'linestyle','-');
        end
    end

% -------------------------------------------------------------------------
    function nShowText(src,evt) %#ok
        
        val = get(src,'checked');
        if strcmp(val,'off')
            set(src,'checked','on');
            set(T,'visible','on');
            flag = true;
        else
            set(src,'checked','off');
            set(T,'visible','off');
            flag = false;
        end
    end
% -------------------------------------------------------------------------
    function nCreate
        % Create a figure and hard code the axis settings.
        spos = get(0,'screensize');
        f = figure('position',[spos(3)/2-300 spos(4)/2-300 600 600],...
            'color',get(0,'DefaultUicontrolBackgroundcolor'),...
            'name','Cox-Ross-Rubenstein Tree',...
            'numbertitle','off');

        uic = uicontextmenu('parent',f);

        if flag
            uimenu(uic,'Label','Show Text',...
                'callback',@nShowText,...
                'checked','on');
        end

        axes('xlim',[0,1],'ylim',[-1 1],'xtick',[],'ytick',[],...
            'box','on',...
            'units','pixels',...
            'position',[50 200 500 350],'uicontextmenu',uic);

        xlabel('Time','fontsize',16);
        ylabel('Share price','fontsize',16);

        Names = {'Share price','Strike price','Interest rate',...
            'Maturity','Volatility'};
        LB = [0,0,0,0,0];
        UB = [100,100,0.1,5,.4];
        Values = [36,40,.06,1,.2];
        xpos = [50 240 430 50 240];
        ypos = [90 90 90 20 20];
        
        for ii = 1:numel(Names)
            nWidget(ii);
        end
        
        uicontrol(f,'style','text',...
            'string','Timesteps',...
            'fontsize',14,...
            'position',[430 60 120 20]);
        uicontrol(f,'style','edit',...
            'string','8',...
            'backgroundcolor',[1 1 1],...
            'callback',{@nTimeSteps},...
            'position',[475 40 30 18]);
        
        % -----------------------------------------------------------------
        function nWidget(ii)
            % Create a widget to control a paramter in the model
            Controls.Slider(ii) = uicontrol(f,'style','slider',...
                'min',LB(ii),'max',UB(ii),...
                'value',Values(ii),...
                'position',[xpos(ii),ypos(ii),120 18],...
                'callback',{@nSlider,ii});
                
            Controls.Min(ii) = uicontrol(f,'style','edit',...
                'backgroundcolor',[1 1 1],...
                'string',num2str(LB(ii)),...
                'position',[xpos(ii),ypos(ii)+20,30,18],...
                'callback',{@nMinimum,ii});
            
            Controls.Value(ii) = uicontrol(f,'style','edit',...
                'backgroundcolor',[1 1 1],...
                'string',num2str(Values(ii)),...
                'position',[xpos(ii)+45,ypos(ii)+20,30,18],...
                'callback',{@nValue,ii});
            
            Controls.Max(ii) = uicontrol(f,'style','edit',...
                'backgroundcolor',[1 1 1],...
                'string',num2str(UB(ii)),...
                'position',[xpos(ii)+90,ypos(ii)+20,30,18],...
                'callback',{@nMaximum,ii});
            
            uicontrol(f,'style','text',...
                'string',Names{ii},...
                'fontsize',14,...
                'position',[xpos(ii),ypos(ii)+40,120,20]);            
            
        end
    end
% -------------------------------------------------------------------------
    function nMinimum(src,evt,Idx)  %#ok
        
        value = str2num(get(src,'string'));
        if value < 0 || isempty(value) || value > Values(Idx)
            set(src,'string',num2str(LB(Idx)));
            return
        end
        set(Controls.Slider(Idx),'min',value);
        LB(Idx) = value;
        
    end
% -------------------------------------------------------------------------
    function nValue(src,evt,Idx)    %#ok
        
        value = str2num(get(src,'string'));
        if value < LB(Idx) || isempty(value) || value > UB(Idx)
            set(src,'string',num2str(Values(Idx)));
            return
        end
        Values(Idx) = value;
        set(Controls.Slider(Idx),'value',value);
        cla;
        nPlotTree;
        
    end
% -------------------------------------------------------------------------
    function nMaximum(src,evt,Idx)  %#ok
        
        value = str2num(get(src,'string'));
        if value < Values(Idx) || isempty(value) || isinf(value)
            set(src,'string',num2str(UB(Idx)));
            return
        end
        UB(Idx) = value;
        set(Controls.Slider(Idx),'max',value);
        
    end
% -------------------------------------------------------------------------
    function nSlider(src,evt,Idx)   %#ok
        
        value = get(src,'value');
        set(Controls.Value(Idx),'string',num2str(value));
        Values(Idx) = value;
        cla;
        nPlotTree;
        
    end

% -------------------------------------------------------------------------
    function nTimeSteps(src,evt) %#ok
        
        Num = str2num(get(src,'string'));
        if isempty(Num) || Num ~= round(Num)
            set(src,'string',num2str(N));
            return
        end
        N = Num;
        cla;
        nPlotTree;
        
    end
        
end
% -------------------------------------------------------------------------
function [Price,P,S,Time,dt,u,d,p,a] = iCalcTree(S0,K,r,T,sigma,N)

dt = T/N;

u = exp(sigma*sqrt(dt)); d = 1/u;
a = exp(r*dt); p = (a-d)/(u-d);

% Create final Returns on the tree
S{N+1} = S0*u^N*d.^(0:2:2*N);
P{N+1} = max(K-S{N+1},0);
Time{N+1} = T*ones(1,N+1);
% Now move back through time and calculate the expected return at previous
% nodes on the tree. Compare this with the immediate return. Exercise the
% option if the immediate return is greater than the expected return

for ii = N:-1:1
    Q = zeros(1,ii);
    V = zeros(1,ii);
    for jj = 1:ii
        % Share price at current node
        V(jj) = S0*u^(ii-1)*d^(2*(jj-1));
        % Expected value of option due if we continue to hold
        E = p*P{ii+1}(jj)/a+(1-p)*P{ii+1}(jj+1)/a;
        % Value of early exercise
        I = max(K-V(jj),0);
        % Value of option at this Node
        Q(jj) = max(E,I);
    end
    S{ii} = V;
    P{ii} = Q;
    Time{ii} = ii*dt*ones(size(S{ii}));
end

Price = P{1}

end

% -------------------------------------------------------------------------
function [S0,K,r,T,sigma,N] = iDefaults

S0 = 36; K = 40; r = 0.06; T = 1; sigma = 0.2; N = 8;

end