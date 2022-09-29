classdef figureHandler < handle
    %FIGUREHANDLER handles the properties of a figure object. 
    %   Provides strucutres and functions to store and restore the
    %   properites of a given figure. 
    
    %% PROPERTIES
    properties (SetAccess=private)
        fig = []; %figure
        axes = []; %(first) axes
        cb = []; %colorbar
        legen = []; %legend
        
        %stored properties
        oldFigProps = struct();
        oldAxesProps = struct();
        oldColorbarProps = struct();
        oldLegendProps = struct();
    end
    

    %#######################################################################################################
    % METHODS - constructor
    %#######################################################################################################
    methods
        function obj = figureHandler(fig)
            %FIGUREHANDLER Construct an instance of this class
            %   fig has to be a figure and must be provided
            
            if nargin~=1 || ~isgraphics(fig,'figure')
                error('figureHandler must be initialized with a figure as argument')
            end
            
            %store the figure handles
            obj.fig = fig;
            
            %find the first axes and store it
            for i=1:numel(fig.Children)
                if isa(fig.Children(i),'matlab.graphics.axis.Axes')
                    obj.axes = fig.Children(i);
                    break
                end
            end            
            if isempty(obj.axes)
                error('no axes found')
            end
    
            %find a colorbar if available
            obj.cb = [];
            for i=1:numel(fig.Children)
                if isa(fig.Children(i),'matlab.graphics.illustration.ColorBar')
                    obj.cb = fig.Children(i);
                    break
                end
            end
            
            %find a legend if available
            obj.legen = [];
            for i=1:numel(fig.Children)
                if isa(fig.Children(i),'matlab.graphics.illustration.Legend')
                    obj.legen = fig.Children(i);
                    break
                end
            end
        end
    end
        
    %#######################################################################################################
    % METHODS - store properties
    %#######################################################################################################
    methods
        function storeProperties(obj)
            obj.storeFigProps;
            obj.storeAxesProps;
            obj.storeColorbarProps;
            obj.storeLegendProps;
        end
        function obj = storeFigProps(obj)
            %stores figure propterties that are modified during the code
            zwFig = obj.fig;
            obj.oldFigProps = struct(...
                'Units',zwFig.Units,...
                'PaperUnits',zwFig.PaperUnits,...
                'PaperSize',zwFig.PaperSize,...
                'PaperPosition',zwFig.PaperPosition,...
                'PaperPositionMode',zwFig.PaperPositionMode,...
                'PaperOrientation',zwFig.PaperOrientation,...
                'Position',zwFig.Position,...
                'InnerPosition',zwFig.InnerPosition,...
                'Color',zwFig.Color,...
                'numChilds',numel(zwFig.Children),...
                'childVisibility',vertcat(zwFig.Children.Visible));
        end
        function obj = storeAxesProps(obj)
            %stores axes propterties that are modified during the code
            zwAxes = obj.axes;
            obj.oldAxesProps = struct(...
                'Units',zwAxes.Units,...
                'Position',zwAxes.Position,...
                'XLim',zwAxes.XLim,...
                'YLim',zwAxes.YLim,...
                'ZLim',zwAxes.ZLim,...
                'CLim',zwAxes.CLim,...
                'XLimMode',zwAxes.XLimMode,...
                'YLimMode',zwAxes.YLimMode,...
                'ZLimMode',zwAxes.ZLimMode,...
                'CLimMode',zwAxes.CLimMode,...
                'XTickMode',zwAxes.XTickMode,...
                'YTickMode',zwAxes.YTickMode,...
                'ZTickMode',zwAxes.ZTickMode,...
                'XTickLabelMode',zwAxes.XTickLabelMode,...
                'YTickLabelMode',zwAxes.YTickLabelMode,...
                'ZTickLabelMode',zwAxes.ZTickLabelMode,...
                'numChilds',numel(zwAxes.Children),...
                'childVisibility',vertcat(zwAxes.Children.Visible),...
                'DataAspectRatioMode',zwAxes.DataAspectRatioMode,...
                'DataAspectRatio',zwAxes.DataAspectRatio,...
                'childDisplayNames',{cell(numel(zwAxes.Children),1)});
            %store the displaynames
            for i=1:numel(zwAxes.Children)
                if figureHandler.isproperty(zwAxes.Children(i),'DisplayName')
                    obj.oldAxesProps.childDisplayNames{i} = zwAxes.Children(i).DisplayName;
                else
                    obj.oldAxesProps.childDisplayNames{i} = '';
                end
            end
        end
        function obj = storeColorbarProps(obj)
            %stores colorbar propterties that are modified during the code
            zwCb = obj.cb;
            if ~isempty(zwCb)
                obj.oldColorbarProps = struct(...
                    'Units',zwCb.Units,...
                    'Position',zwCb.Position,...
                    'AxisLocationMode',zwCb.AxisLocationMode,...
                    'Limits',zwCb.Limits,...
                    'LimitsMode',zwCb.LimitsMode);
            end
        end
        function obj = storeLegendProps(obj)
            %stores colorbar propterties that are modified during the code
            zwLegen = obj.legen;
            if ~isempty(zwLegen)
                obj.oldLegendProps = struct(...
                    'AutoUpdate',zwLegen.AutoUpdate,...
                    'Color',zwLegen.Color,...
                    'Box',zwLegen.Box,...
                    'TextColor',zwLegen.TextColor,...
                    'LineWidth',zwLegen.LineWidth,...
                    'EdgeColor',zwLegen.EdgeColor,...
                    'FontSize',zwLegen.FontSize,...
                    'FontWeight',zwLegen.FontWeight,...
                    'FontName',zwLegen.FontName,...
                    'FontAngle',zwLegen.FontAngle,...
                    'TitleStruct',zwLegen.Title,...
                    'NumColumns',zwLegen.NumColumns,...
                    'NumColumnsMode',zwLegen.NumColumnsMode,...
                    'Visible',zwLegen.Visible,...
                    'Interpreter',zwLegen.Interpreter,...
                    'Position',zwLegen.Position,...
                    'Units',zwLegen.Units,...
                    'Orientation',zwLegen.Orientation,...
                    'Location',zwLegen.Location);
            end
        end
    end
    
    %#######################################################################################################
    % METHODS - restore properties
    %#######################################################################################################
    methods
        function restoreGraphicsState(obj)
            %restore the properties
            obj.restoreFigureProps;
            obj.restoreAxesProps;
            obj.restoreColorbarProps;
            obj.restoreLegendProps;           
            
            %restore the visibility
            for i=1:numel(obj.fig.Children)
                obj.fig.Children(i).Visible = obj.oldFigProps.childVisibility(i,:);
            end
            for i=1:numel(obj.axes.Children)
                obj.axes.Children(i).Visible = obj.oldAxesProps.childVisibility(i,:);
            end
        end
        function restoreFigureProps(obj)
            %restore all properties of the saved figure
            if ~isempty(obj.fig)
                obj.restoreProps(obj.fig,obj.oldFigProps);
                for i=1:numel(obj.fig.Children)
                    obj.fig.Children(i).Visible = obj.oldFigProps.childVisibility(i,:);
                end
            end
        end
        function restoreAxesProps(obj)
            %restore all properties of the saved axes
            if ~isempty(obj.axes)                
                %standard properties
                obj.restoreProps(obj.axes,obj.oldAxesProps);
                
                %properties of children
                for i=1:numel(obj.axes.Children)
                    obj.axes.Children(i).Visible = obj.oldAxesProps.childVisibility(i,:);
                end
            end
        end
        function restoreColorbarProps(obj)
            %restore all properties of the saved colorbar
            if ~isempty(obj.cb)
                obj.restoreProps(obj.cb,obj.oldColorbarProps);
            end
        end
        function restoreLegendProps(obj)        
            %restore all properties of the saved legend
            if ~isempty(obj.legen)
                obj.restoreProps(obj.legen,obj.oldLegendProps);
                if ~isempty(obj.oldLegendProps.TitleStruct.String)
                    proptit = obj.oldLegendProps.TitleStruct;
                    title(obj.legen,proptit.String,'FontSize',proptit.FontSize,'FontSizeMode',proptit.FontSizeMode,...
                        'FontAngle',proptit.FontAngle,'FontAngleMode',proptit.FontAngleMode,...
                        'FontName',proptit.FontAngle,'FontNameMode',proptit.FontAngleMode,...
                        'FontWeight',proptit.FontWeight,'FontWeightMode',proptit.FontWeightMode,...
                        'Color',proptit.Color,'ColorMode',proptit.ColorMode,...
                        'Interpreter',proptit.Interpreter,'InterpreterMode',proptit.InterpreterMode,...
                        'Visible',proptit.Visible);
                end
            end
        end
        function restoreLegendNames(obj)
            %restore the names in the legend
            if ~isempty(obj.axes)
                if ~isempty(obj.legen)
                    zwUpdate = obj.legen.AutoUpdate;
                    obj.legen.AutoUpdate = 'on';
                end
                for i=1:numel(obj.axes.Children)
                    if isporperty(obj.axes.Children(i),'DisplayName')
                        obj.axes.Children(i).DisplayName = obj.oldAxesProps.childDisplayNames{i};
                    end
                end
                %reset legend updating 
                if ~isempty(obj.legen)
                    obj.legen.AutoUpdate = zwUpdate;
                end
            end
        end
    end
    methods(Static,Access=private)
        function restoreProps(obj,oldProps)
            props = fields(oldProps);
            for i=1:numel(oldProps)
                if figureHandler.isproperty(obj,props{i})
                    obj.(props{i}) = oldProps.(props{i});
                end
            end
        end
        
        function out = isproperty(obj,field)
            % isproperty: Analogues to "isfield"
            %
            % 03.01.2012 C.Hesch

            % Expand list of properties
            mc = metaclass(obj);
            try
                % Matlab 2011a
                pl = mc.PropertyList;
                nProps = numel(pl);
                [props{1:nProps}] = deal(pl.Name);
            catch
                % Matlab 2010b
                pl = mc.Properties;
                nProps = numel(pl);
                props = {};
                for ii = nProps:-1:1
                    props{ii} = pl{ii}.Name;
                end
            end

            if isa(field,'char')
                out = any(strcmpi(props,field));
            elseif isa(field,'cell')
                out = false(size(field));
                for ii = 1:numel(field)
                    out(ii) = any(strcmpi(props,field{ii}));
                end
            else
                error('input type not supported');
            end
        end

    end
    
    %#######################################################################################################
    % METHODS - axis limits
    %#######################################################################################################
    methods
        function setAxisLimits(obj)
            %std axes
            zwAxes = obj.axes;
            [xlim,ylim,zlim] = obj.getAxisLimits();
            zwAxes.XLim = xlim;
            zwAxes.YLim = ylim;
            zwAxes.ZLim = zlim;
            %colorbar
            if ~isempty(obj.cb)
                obj.cb.Limits = obj.cb.Limits;
                zwAxes.CLim = zwAxes.CLim;
            elseif strcmpi(zwAxes.CLimMode,'manual')
                zwAxes.CLim = zwAxes.CLim;
            end
        end
        function [xlim,ylim,zlim] = getAxisLimits(obj)
            xlim = obj.axes.XLim;
            ylim = obj.axes.YLim;
            zlim = obj.axes.ZLim;

            [xlimData,ylimData,zlimData] = getDataLimits(obj);

            %handle inf in axes.Lim
            xlim(isinf(xlim)) = xlimData(isinf(xlim));
            ylim(isinf(ylim)) = ylimData(isinf(ylim));
            zlim(isinf(zlim)) = zlimData(isinf(zlim));
        end
        function [xlim,ylim,zlim] = getDataLimits(obj)
            zwAxes = obj.axes;
            xlim = [Inf,-Inf];
            ylim = [Inf,-Inf];
            zlim = [Inf,-Inf];

            %find minimum of all data objects of given axes
            for ii=1:numel(zwAxes.Children)
                if isprop(zwAxes.Children(ii),'XData') && ~isempty(zwAxes.Children(ii).XData)
                    xlim(1) = min([xlim(1),min(zwAxes.Children(ii).XData(:))]);
                    xlim(2) = max([xlim(2),max(zwAxes.Children(ii).XData(:))]);

                    ylim(1) = min([ylim(1),min(zwAxes.Children(ii).YData(:))]);
                    ylim(2) = max([ylim(2),max(zwAxes.Children(ii).YData(:))]);

                    zlim(1) = min([zlim(1),min(zwAxes.Children(ii).ZData(:))]);
                    zlim(2) = max([zlim(2),max(zwAxes.Children(ii).ZData(:))]);
                end
            end

            %handle empty data
            xlim = obj.fixInfLimit(xlim);
            ylim = obj.fixInfLimit(ylim);
            zlim = obj.fixInfLimit(zlim);
        end
    end
    methods(Static,Access=private)
        function out = fixInfLimit(lim)
            out = lim;
            if isinf(lim(1))
                out(1) = 0;
            end
            if isinf(lim(2))
                out(2) = 1;
            end
        end
    end
end

