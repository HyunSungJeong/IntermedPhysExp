function varargout = ParseViscosData(N,R,T,varargin)
    % <Description>
    % Parses raw Tracker data from Brownian motion experiment and gives 
    % processed data such as x, y displacement, average of distance sqaured, etc.
    %
    % <Input>
    % N : [integer] Number of tracked particles data
    % R : [numeric] Radius of microscopic particle used in experiment, in units of micrometer
    % T : [numeric] Temperature in Kelvins
    % 
    % <Options>
    % '-v' : If used, the <r^2> vs time is plotted together wit a linear fit
    %       (Default: not used)
    % 'CorrectDrift' ... : If used, the average drift of the particles are subtracted from the original data
    % 
    %                      When used together with a input of 1xN cell array of two dimensional vectors,
    %                      These vectors are regarded drift velocities of each particle,
    %                      and are used to correct drift of individual particles
    %                      (Default: not used)
    % 'MaxNumData', ... : [integer] maximum number of sequential data for each particle
    %       (Default: 1500)
    %
    % <Output>
    % There are three possible outputs:
    % 1. [time,X,Y,Viscosity]
    % 2. [time,X,Y,Viscosity,DistSqAvg,DistSqVar]
    % 3. [time,X,Y,Viscosity,DistSqAvg,DistSqVar,DriftVel] : can only be used when 'CorrectDrift' option is on
    %
    % Time : [cell array of numeric vectors] Each cell element is a vector of frame times, for each particle. 
    % X : [cell array of numeric vectors] 1xN cell array. 
    %           Each cell element is a numeric vector of the x-displacement of each particle 
    %           at times given by output variable 'time'
    % Y : [cell array of numeric vectors] 1xN cell array. 
    %           Each cell element is a numeric vector of the y-displacement of each particle 
    %           at times given by output variable 'time'
    % Viscosity : [numeric] Calculated viscosity of the medium
    % DistSqAvg : [numeric vector] Numeric vector of the averages of displacemet squared at times give by output variable 'time'
    % DistSqVar : [numeric vector] Numeric vector of the variance of displacement squared at times given by output variable 'time'

    %% Parse inputs
    if ~isnumeric(N)
        error('ERR: ''N'' must be a positive integer');
    elseif N <= 0
        error('ERR: ''N'' must be a positive integer');
    end

    if ~isnumeric(R)
        error('ERR: ''R'' must be a positive real number');
    elseif R <= 0
        error('ERR: ''R'' must be a positive real number');
    end

    if ~isnumeric(T)
        error('ERR: ''T'' must be a positive real number');
    elseif T <= 0
        error('ERR: ''T'' must be a positive real number');
    end

    %% Parse options

    % Default values of options
    plotFig = false;
    CorrectDrift_default = false;
    CorrectDrift_custom = false;
    MaxNumData = 1500;

    while ~isempty(varargin)
        switch varargin{1}
            case '-v'
                plotFig = true;
                varargin(1) = [];

            case 'CorrectDrift'
                if iscell(varargin{2})
                    CorrectDrift_custom = true;
                    DriftVel = varargin{2};
                    varargin(1:2) = [];
                else
                    CorrectDrift_default = true;
                    varargin(1) = [];
                end

            case 'MaxNumData'
                if isnumeric(varargin{2})
                    if varargin{2} > 0
                        MaxNumData = varargin{2};
                        varargin(1:2) = [];
                    else
                        error('ERR: ''MaxNumData'' must be a positive integer');
                    end
                else
                    error('ERR: ''MaxNumData'' must be a positive integer');
                end

            otherwise
                if ischar(varargin{1})
                    error(['ERR: Unknown option ''',varargin{1},'''']);
                else
                    error('ERR: Unknown input');
                end
        end
    end


    %% Parse raw data

    time = linspace(0,(MaxNumData-1)/10,MaxNumData);
    kb = 1.38*1e-23;
    
    X = cell(1,N);
    Y = cell(1,N);
    DistSq = cell(1,N);
    DataIdx = cell(1,N);
    
    for it = 1:N
        data = importdata(['C:\Users\82104\Documents\서울대학교\2025\중급물리실험 1\Optical Tweezers\',sprintf('%d',R),'um\Viscosdata',num2str(it),'.txt']);   % load data
        if isstruct(data)
            data = data.data;
        end
        DataIdx{it} = nan(1,size(data,1));
        X{it} = nan(1,MaxNumData);
        Y{it} = nan(1,MaxNumData);
        DistSq{it} = nan(1,MaxNumData);
    
        skipIdx = [];
        for itD = 1:size(data,1)
            Time_now = data(itD,1);
            X_now = data(itD,2);
            Y_now = data(itD,3);
            
            if itD > 1
                X_prev = data(itD-1,2);
                Y_prev = data(itD-1,3);
                Time_prev = data(itD-1,1);
    
                if sqrt((X_now-X_prev)^2 + (Y_now-Y_prev)^2) > 1e-6     % if there is a sudden jump of coordinates
        
                    % calculate the amout of sudden jump
                    TimeJump = Time_now - Time_prev;
                    Xjump = X_now - X_prev;
                    Yjump = Y_now - Y_prev;
        
                    skipIdx = [skipIdx, itD];   % indices to be skipped when processing data into displacements

                    % subtract sudden jump from following data
                    data(itD:end,1) = data(itD:end,1)-  repmat(TimeJump,[size(data,1)-itD+1, 1]);
                    data(itD:end,2) = data(itD:end,2) - repmat(Xjump,[size(data,1)-itD+1, 1]);
                    data(itD:end,3) = data(itD:end,3) - repmat(Yjump,[size(data,1)-itD+1, 1]);
                end
            end
        end % itD
    
        for itD = 1:size(data,1)        % process data into displacements
            if ~ismember(itD,skipIdx)   % if there was no sudden jump in this frame
                DataIdx{it}(itD) = round(10*data(itD,1)) + 1;
                X_tmp = 1e6*(data(itD,2) - data(1,2));
                Y_tmp = 1e6*(data(itD,3) - data(1,3));
                X{it}(DataIdx{it}(itD)) = X_tmp;
                Y{it}(DataIdx{it}(itD)) = Y_tmp;
            end
        end
    
    end % it
    

    %% Correct overall drift of particles

    DriftAvg_X = 0;
    DriftAvg_Y = 0;
    Drift_num = 0;
    for it = 1:N
        if ~isnan(X{it}(900))
            DriftAvg_X = DriftAvg_X + X{it}(900);
            DriftAvg_Y = DriftAvg_Y + Y{it}(900);
            Drift_num = Drift_num + 1;
        else
            error('Choose another time');
        end
    end
    DriftAvg_X = DriftAvg_X/Drift_num;
    DriftAvg_Y = DriftAvg_Y/Drift_num;
    if CorrectDrift_default
        DriftVel = [DriftAvg_X/90, DriftAvg_Y/90];
    elseif ~CorrectDrift_custom
        DriftVel = [0,0];
    end
    
    
    for it = 1:N
        if CorrectDrift_custom
            X{it} = X{it} - DriftVel{it}(1)*time;
            Y{it} = Y{it} - DriftVel{it}(2)*time;
        else
            X{it} = X{it} - DriftVel(1)*time;
            Y{it} = Y{it} - DriftVel(2)*time;
        end
        DistSq{it} = X{it}.^2 + Y{it}.^2;
    end
    
    DistSqAvg = nan(1,MaxNumData);
    DistSqVar = nan(1,MaxNumData);
    for itD = 1:MaxNumData
        cnt = 0;
        for it = 1:N
            if ~isnan(X{it}(itD))
                if isnan(DistSqAvg(itD))
                    DistSqAvg(itD) = DistSq{it}(itD);
                    DistSqVar(itD) = DistSq{it}(itD)^2;
                end
                DistSqAvg(itD) = DistSqAvg(itD) + DistSq{it}(itD);
                DistSqVar(itD) = DistSqVar(itD) + DistSq{it}(itD)^2;
                cnt = cnt + 1;
            end
        end
        if cnt > 0
            DistSqAvg(it) = DistSqAvg(it)/cnt;
            DistSqVar(it) = DistSqVar(it)/cnt - DistSqAvg(it);
        end
    end
    
    %% Calculate viscosity from processed data

    % discard frames without data
    
    Time = cell(1,N);
    for it = 1:N
        Time{it} = time;
        EmptyIdx = isnan(X{it});
        Time{it}(EmptyIdx) = [];
        X{it}(EmptyIdx) = [];
        Y{it}(EmptyIdx) = [];
        DistSq{it}(EmptyIdx) = [];
    end
    
    EmptyIdx = isnan(DistSqAvg);
    DistSqAvg(EmptyIdx) = [];
    DistSqVar(EmptyIdx) = [];
    time(EmptyIdx) = [];    

    %{
    for it = 1:N
        figure;
        hold on;
        plot(Time{it},DistSq{it});
    end
    %}
    
    % linear fit <r^2> - time data and calculate viscosity
    [coeff, Rsq] = linFit(time(500:end), DistSqAvg(500:end));
    Viscosity = (1e18)*2*kb*T/(3*pi*R*coeff(1));
    Stdev = sqrt(DistSqVar);
    

    %% Plot figure if '-v' option was used

    if plotFig
        figure;
        hold on;
        errorbar(time(1:10:end),DistSqAvg(1:10:end),Stdev(1:10:end),'-.','CapSize',10,'LineWidth',0.8);
        plot(time,DistSqAvg,'-','color','red','LineWidth',1);
        
        
        TimeFit = linspace(time(1),time(end),500);
        DistFit = polyval(coeff, TimeFit);
        plot(TimeFit,DistFit,'--','Color','black','LineWidth',1.5);
        linFitEq = ['$r^{2} = ', sprintf('%.3g',coeff(1)),' t ',sprintf('%+.3g',coeff(2)),'$'];
        RsqText = ['$ R^{2} = ',sprintf('%.4g',Rsq),'$'];
        text(50, 70, linFitEq,'Interpreter','latex','FontSize',15);
        text(50, 90, RsqText,'Interpreter','latex','FontSize',15);
        text(50, 110,['$\eta = ',sprintf('%.5g',Viscosity),'$'],'Interpreter','latex','FontSize',15);
        
        set(gca,'XScale','linear','YScale','linear','fontsize',20);
        xlabel('$\mathrm{time} \left( s \right)$','Interpreter','latex','FontSize',25);
        ylabel('$\mathrm{distance^{2}} \left( \mu \mathrm{m}^{2} \right)$','Interpreter','latex','FontSize',25);
        title('$\mathrm{Brownian \ motion}$','Interpreter','latex','FontSize',30);
        hold off;
    end
    
    switch nargout
        case 4
            varargout = {Time, X, Y, Viscosity};
        case 6
            varargout = {Time, X, Y, Viscosity, DistSqAvg, DistSqVar};
        case 7
            if CorrectDrift_default
                varargout = {Time, X, Y, Viscosity, DistSqAvg, DistSqVar, DriftVel};
            else
                error('ERR: ''DriftVel'' output can only be output when ''CorrectDrift'' option is used without input drift velocity');
            end
        otherwise
            error('ERR: number of output variables must be either 4, 6, or 7');
    end
    

    function [coeff, Rsq] = linFit(time, dist)
    
        coeff = polyfit(time,dist,1);
        yfit = polyval(coeff, time);            % Estimated Regression Line
        SStot = sum((dist-mean(dist)).^2);      % Total Sum-Of-Squares
        SSres = sum((dist-yfit).^2);            % Residual Sum-Of-Squares
        Rsq = 1-SSres/SStot;                    % R^2
    end
end