function [DriftVel,GaussianFit] = GaussianScatter(Time,X,Y,varargin)
    % <Description>
    % Obtains average drift velocity for each particle, by obtaining the
    % average velocity of all successive brownain motion
    %
    % <Input>
    % Time : [cell array of numeric vectors] Each cell element is a vector of frame times, for each particle. 
    % X : [cell array of numeric vectors] 1xN cell array. 
    %           Each cell element is a numeric vector of the x-displacement of each particle 
    %           at times given by output variable 'time'
    % Y : [cell array of numeric vectors] 1xN cell array. 
    %           Each cell element is a numeric vector of the y-displacement of each particle 
    %           at times given by output variable 'time'
    %
    % <Option>
    % '-v' : If used, the scatter plot of velocities of all successive
    %        frames are displayed for each particle
    %       (Default: not used)
    %
    % <Output>
    % DriftVel : [cell array of numeric vectors] Cell array of 1x2 numeric vectors.
    %               Each cell element is the drift velocity of each particle (Vx,Vy).
    % GuassianFit : [Cell array of numeric vectors] Cell array of 1x2 numeric vectors.
    %               Each cell element if the R^2 value of the Gaussian fit
    %               for the x- and y- components of velocities for each particle

    %% Parse input data
    if ~iscell(Time)
        error('ERR: ''Time'' must be a cell array');
    end

    if ~iscell(X)
        error('ERR: ''X'' must be a cell array');
    end

    if ~iscell(Y)
        error('ERR: ''Y'' must be a cell array');
    end

    if isequal(numel(Time), numel(X), numel(Y))
        N = numel(Time);
    else
        error('ERR: ''Time'', ''X'', and ''Y'' must be cell arrays of same lengths');
    end

    for it = 1:N
        if ~isequal(numel(Time{it}), numel(X{it}), numel(Y{it}))
            error('ERR: each elements of ''Time'', ''X'', and ''Y'' must be numeric arrays of same lengths');
        end
    end

    %% Parse options
    
    % Default values of options
    PlotScatter = false;

    while ~isempty(varargin)
        switch varargin{1}
            case '-v'
                PlotScatter = true;
                varargin(1) = [];

            otherwise
                if ischar(varargin{1})
                    error(['ERR: Unknown option ''',varargin{1},'''']);
                else
                    error('ERR: Unknown input');
                end
        end
    end

    %% Parse displacement data

    Vx = cell(1,N);
    Vy = cell(1,N);
    GaussianFit = cell(1,N);
    for it = 1:N
        Vx{it} = zeros(1,numel(Time{it})-1);
        Vy{it} = zeros(1,numel(Time{it})-1);
        for itD = 1:numel(Time{it})-1
            Vx{it}(itD) = (X{it}(itD+1) - X{it}(itD))/(Time{it}(itD+1) - Time{it}(itD));
            Vy{it}(itD) = (Y{it}(itD+1) - Y{it}(itD))/(Time{it}(itD+1) - Time{it}(itD));
        end

        %[GaussCoeff_X,Rsq] = GaussFit(Vx{it});
        %[GaussCoeff_Y,Rsq] = GaussFit(Vy{it});
        DriftVel{it} = [mean(Vx{it}), mean(Vy{it})];
    end
    
    %% Draw scatter plot if '-v' option was used
    for it = 1:N
        figure;
        hold on;
        WinSize = 10;
        xlim([-WinSize,WinSize]);
        ylim([-WinSize,WinSize]);
        for itD = 1:numel(Vx)
            plot(Vx{it},Vy{it},'.','Color','blue','Markersize',10);
        end
        plot(DriftVel{it}(1),DriftVel{it}(2),'.','color','red','MarkerSize',20);
        plot(DriftVel{it}(1),DriftVel{it}(2),'x','color','red','MarkerSize',20);
        plot([-WinSize,WinSize],[0,0],'--','Color','black','LineWidth',1);
        plot([0,0],[-WinSize,WinSize],'--','Color','black','LineWidth',1);
        hold off;
    end



    function [GaussCoeff,Rsq] = GaussFit(ProbData)
        Grid = linspace(-10,10,101);
        Weight = zeros(1,numel(Grid));
        for itData = 1:numel(ProbData)
            IntervIdx = floor(5*(ProbData(itData)+10)) + 1;
            Weight(IntervIdx) = Weight(IntervIdx) + 1;
        end
        logWeight = log(Weight);
        NonzeroIdx = Weight>0;
        Grid = Grid(NonzeroIdx);
        logWeight = logWeight(NonzeroIdx);

        coeff = polyfit(Grid,logWeight,2);
        a = coeff(1);
        b = coeff(2);
        c = coeff(3);
        GaussCoeff = zeros(1,3);
        GaussCoeff(1) = exp(c-(b^2/(4*a)));
        GaussCoeff(2) = -b/(2*a);
        GaussCoeff(3) = 1/sqrt(-a);

        yfit = polyval(coeff, Grid);            % Estimated Regression Line
        SStot = sum((logWeight-mean(logWeight)).^2);      % Total Sum-Of-Squares
        SSres = sum((logWeight-yfit).^2);            % Residual Sum-Of-Squares
        Rsq = 1-SSres/SStot;                    % R^2

    end


    function circle(x,y,r)
        th = 0:pi/50:2*pi;
        xunit = r * cos(th) + x;
        yunit = r * sin(th) + y;
        plot(xunit, yunit);
    end
end