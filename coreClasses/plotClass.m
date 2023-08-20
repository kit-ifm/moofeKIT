classdef plotClass < matlab.mixin.Copyable
    % settings for plotting 
    properties
        % settings if and when the system is plotted
        flag                    = true;
        savePlot                = struct('flag',false, ...
                                         'name','saveName',...
                                         'type','-depsc',...
                                         'flagPlotColorbar',false) % '-depsc'

        steps                   = 1; % number of steps after which solution is plotted
        everyNewtonStep         = true;
        docked                  = true;
        keepFormerPlots         = false;
        plotInitialConfig       = false;
        % parts of drawing that should be plotted
        boundary            = true;
        load                = true;
        contact             = true;
        axis                = true;
        colorData           = 'none';
        postPlotType = 'stress'; % 'temp'; 'phi'; 'D1'; 'D2'; 'D3'; 'none'
        stress = struct('type','Cauchy',... %'FirstPK'
                        'component',-1); % -1 = von Mises, 11, 22, 33, 12, 21 etc.
        % variables for formatting
        lineWidth = 1;
        lineWeight              = 0.6; %FIXME: ???
        border                  = 0.3;
        boundarySize            = 0.03;
        colorBarLimits          = [];
        view = [-217.5, 30];
        makeMovie = false;
        time = 'N1';%'R';
    end
    properties (SetAccess = private) %these properties are calculated
        xmin = zeros(3,1);
        xmax = zeros(3,1);
        xminxmaxset = false;
        deltaXY = 0;
        maxF = 0; 
    end
    methods
        %% set methods
        function set.flag(obj,val)
            assert(islogical(val),'plot must be logical')
            obj.flag = val;
        end
        function set.steps(obj,val)
            assert(floor(val)==ceil(val) && val>0,'plotSteps must be an Integer>=1')
            obj.steps = val;
        end
        function set.boundary(obj,val)
            assert(islogical(val),'plotBoundary must be logical')
            obj.boundary = val;
        end
        function set.load(obj,val)
            assert(islogical(val),'plotLoad must be logical')
            obj.load = val;
        end
        function set.contact(obj,val)
            assert(islogical(val),'plotContact must be logical')
            obj.contact = val;
        end
        function set.everyNewtonStep(obj,val)
            assert(islogical(val),'plotEveryNewton must be logical')
            obj.everyNewtonStep = val;
        end
        function set.axis(obj,val)
%             assert(islogical(val),'plotAxis must be logical')
            obj.axis = val;
        end
        function set.colorData(obj,val)
            assert(strcmpi(val,'none')||strcmpi(val,'stress')||strcmpi(val,'plasticstrain'), ...
                'plotColorData must be: none, stress or plasticstrain')
            obj.colorData = val;
        end
        function set.colorBarLimits(obj,val)
            assert((size(val,1) == 1 && size(val,2) == 2) || isempty(val),'colorBarLimits must have size 1x2 or be empty')
            obj.colorBarLimits = val;
        end
        function set.makeMovie(obj,val)
            assert(islogical(val),'makeMovie must be logical')
            obj.makeMovie = val;
        end
        %% setting min and max values for plotting
        function setXminXmax(obj,xmin,xmax)
            % check input
            if size(xmin,1)~=3 || size(xmin,2)~=1 || size(xmax,1)~=3 || size(xmax,2)~=1
                error('xmin and xmax must have dimension 3x1')
            end
            if any(xmin>xmax)
                error('xmin must be > xmax')
            end
            obj.xmin = xmin;
            obj.xmax = xmax;
            obj.deltaXY = max(xmax-xmin);
            obj.xminxmaxset = true;
        end
        %% finding min and max values for plotting
        function findXminXmax(obj,solidObject)
            if ~obj.xminxmaxset
                for ii=1:numel(solidObject)

                    dimension = solidObject(ii).dimension;
                    if isa(solidObject,'stringClass') 
                        dimension = size(solidObject.qR,2);
                    end

                    for jj=1:dimension
                        % minimum and maximum
                        obj.xmin(jj) = min([obj.xmin(jj); solidObject(ii).meshObject.nodes(:,jj)]);
                        obj.xmax(jj) = max([obj.xmax(jj); solidObject(ii).meshObject.nodes(:,jj)]);
                        obj.deltaXY = max(obj.deltaXY,obj.xmax(jj)-obj.xmin(jj)); 
                    end
                end 
            end
        end
        function findFmax(obj,neumannObject)
            for ii=1:numel(neumannObject)
                % maximum force
                obj.maxF = max([obj.maxF; sqrt(sum(neumannObject(ii).nodalForce.^2,2))]);
            end 
        end
    end
end