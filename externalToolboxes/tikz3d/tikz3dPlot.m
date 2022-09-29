function tikz3dPlot(varargin)
    %TIKZ3DPLOT creates a tikz 3d plot with external graphics fitted into a pgf
    %axis system. It uses the matlab tools "matlab2tikz" and "export_fig" 
    %https://de.mathworks.com/matlabcentral/fileexchange/22022-matlab2tikz-matlab2tikz
    %https://de.mathworks.com/matlabcentral/fileexchange/23629-export_fig
    %
    %TIKZ3DPLOT options
    %All options of matlab2tikz are supproted and can be applied completely 
    %analogously (type "help matlab2tikz" into the command window for a list of 
    %all options).
    %
    %TIKZ3DPLOT(FILENAME,...) or TIKZ3DPLOT('filename',FILENAME,...) stores the
    %LaTeX code and additional pictures under FILENAME. If 'standalone' option
    %is activated in matlab2tikz suffix '_pic' is added to additional picture.
    %Default: 'tikz3dPlot'.
    %
    %TIKZ3DPLOT('picformat',PICFORMAT,...) or TIKZ3DPLOT(PICFORMAT,...) sets
    %the fileformat for the external picture. Possibilities pdf {'pdf',
    %'-dpdf'}, png {'png', '-dpng'}, eps {'eps', '-depsc', '-depsc2'} or tikz
    %{'tikz', '-tikz'}. If tikz is chosen the 3d plot is created using 
    %matlab2tikz. Default: 'pdf' 
    %
    %TIKZ3DPLOT('figurehandle',FIGUREHANDLE,...) or TIKZ3DPLOT(FIGUREHANDLE,...) 
    %explicitly specifies the handle of the figure that is to be stored. 
    %Default: gcf
    %
    %TIKZ3DPLOT('-painters',...) or TIKZ3DPLOT('-opengl',...) sets the renderer
    %for creation of the external figure. Use '-painters' for real vector
    %graphics. Default '-painters'
    %
    %TIKZ3DPLOT('extraPlotCode',CHAR or CELLCHAR,...) explicitly adds extra code
    %before any other plot commands of the output file. Default: []
    %
    %TIKZ3DPLOT('extraPlotCodeEnd',CHAR or CELLCHAR,...) explicitly adds extra code
    %after any other plot commands of the output file. Default: []
    %
    %TIKZ3DPLOT('use_export_fig',...) use export_fig for export of the
    %png/eps/pdf image. Has no effect if 'picformat' is 'tikz'.
    %
    %TIKZ3DPLOT('debug',...) activates debug mode and shows points used for
    %tikz in 3d space and on canvas. They should all lie above of each other.
    %Works only for 3d images.
    %
    %TIKZ3DPLOT('magnify',FACTOR,...) implements magnification of bitmap
    %pictures. Works only if default exportFig is used in stead of 
    %export_fig. Default: 1

    %read the input
    opts = inputParser(varargin{:});
    
    %check if necessary programs exist
    if exist('matlab2tikz','file')~=2
        error('tikz3dPlot:missingExternalProgram','could not find essential function "matlab2tikz" for tikz3dPlot')
    end
    if exist('export_fig','file')~=2
        error('tikz3dPlot:missingExternalProgram','could not find essential function "export_fig" for tikz3dPlot')
    end

    %call the required routine to create files
    switch lower(opts.picformat)
        case 'tikz'
            matlab2tikz([fullfile(opts.filefolder,opts.filename),'.tex'],opts.matlab2tikzopts{:}) 
            tikzPoints = struct();
        case {'eps','png','pdf'}
            tikzPoints = tikz3dPlotWithGraphicFile(opts);           
        otherwise
            error('not implemented')
    end
    
    %amend matlab2tikz file
    amendMatlab2tikzFile(opts,tikzPoints);
end

%% main function to export graphics file
function out = tikz3dPlotWithGraphicFile(opts)        
    %prepare handles and store old settings
    figHandler = figureHandler(opts.figurehandle);
    figHandler.storeProperties;
    fig = figHandler.fig;
    axes = figHandler.axes;
    
    %fix axis limits
    figHandler.setAxisLimits;
    
    %determine if 3d axes
    axis3d = isAxis3D(axes);
    
    %step 1 - export tikz axes with matlab2tikz
    %--------------------------------------------------------------
    % create simple objects for legend (since they are otherwise not
    % available to matlab2tikz)
    makeTikzLegend(figHandler,opts.tikz3dPlotText);
    
    %switch all graphics objects off
    switchOFFgraphicsObjects(axes)
    switchONTikzLegendObjects(axes,opts.tikz3dPlotText)
    try 
        matlab2tikz([fullfile(opts.filefolder,opts.filename),'.tex'],opts.matlab2tikzopts{:})
    catch e
        %try to restore the graphics before exiting
        try 
            deleteTikzLegendObjects(figHandler,opts.tikz3dPlotText,legendTemp)
            figHandler.restoreGraphicsState(); 
        catch 
        end
        error('tikz3dPlot:externalProgramFailed',['tikz3dPlot - error in matlab2tikz\n',e.message])
    end
    deleteTikzLegendObjects(figHandler,opts.tikz3dPlotText)
    switchONgraphicsObjects(axes,figHandler.oldAxesProps.childVisibility)
    
    %step 2 - export the content
    %--------------------------------------------------------------
    %switch everything of in figure 
    switchOFFgraphicsObjects(fig);
        
    %export the graphic 
    %modifications according to pgfplotsmanual
    if strcmpi(opts.picformat,'pdf')
        zwUnit = fig.Units;
        fig.Units = fig.PaperUnits;
        fig.PaperSize = fig.Position(3:4);
        fig.Units = zwUnit;
    end
    
    %export using export fig
    if axis3d
        crop = [];
    else
        %determine the cropping range
        axes.Units = 'pixels';
        fig.Units = 'pixels';
        posAxes = axes.Position;
        posFig = fig.InnerPosition;
        posFig(1:2) = 1;
        
        %fix for given aspect ratio (position of axes modified)
        %https://de.mathworks.com/help/matlab/creating_plots/automatic-axes-resize.html
        if strcmpi(axes.DataAspectRatioMode,'manual')
            [xlim,ylim,~] = figHandler.getAxisLimits;
            scalFactors = posAxes(3:4)./[xlim(2)-xlim(1), ylim(2)-ylim(1)];
            scalFactors = min(scalFactors)./scalFactors;
        else
            scalFactors = [1,1];
        end
        deltaCrop = (posAxes(3:4)).*(1-scalFactors)/2;
        
        %correction of crop margin
        if strcmpi(opts.picformat,'png')
            margin = [0,0,0,0];
        else
            margin = [4,4,4,3];
        end           
        
        %compute the crop margin - export_fig requires
        %[top,right,bottom,left]
        crop = [...
            posFig(4)-posAxes(4)-posAxes(2)+deltaCrop(2)+margin(1);
            posFig(3)-posAxes(3)-posAxes(1)+deltaCrop(1)+margin(2);
            posAxes(2)-posFig(2)+deltaCrop(2)+margin(3);
            posAxes(1)-posFig(1)+deltaCrop(1)+margin(4)];
        crop = max(crop*opts.magnify,0);%do not allow negative crop range
    end
    
    %actually create the file
    [err,posCroppedRange] = createGraphicOfContent(fig,opts,crop);
    
    %handle errors of export_fig
    if ~isempty(err)
        try figHandler.restoreGraphicsState(); catch, end
        if opts.useexportfig
            error('tikz3dPlot:externalProgramFailed',['tikz3dPlot - error in export_fig\n',err.message])
        else
            error('tikz3dPlot:exportFigFailed',['tikz3dPlot - error in exportFig\n',err.message])
        end
    end
    
    %cannot switch on graphics objects now since this would possible
    %distort the graphic leading to faulty tikzpoints
    %switchONgraphicsObjects(fig,oldFigProps.childVisibility)
    
    %step 3 - compute the tikz points (3d) and corresponding canvas points
    %--------------------------------------------------------------
    if axis3d
        [tikzPointsReal,tikzPointsCanv] = supportPointsTikz(axes,posCroppedRange/opts.magnify);
    else
        tikzPointsReal = [axes.XLim',axes.YLim'];
        tikzPointsCanv = [];
    end
    
    %restore the old state of the figure (completely!)
    figHandler.restoreGraphicsState();
    
    %draw points on figure
    if opts.debug==true 
        if axis3d
            hold on
            projectionMatlab(axes.CameraTarget',axes,true);
            for i=1:size(tikzPointsReal,2)
                projectionMatlab(tikzPointsReal(:,i),axes,true);
            end
        else
            hold on
            drawAxesPosAnnotation(fig,axes,{'color','r','linewidth',2})
        end
    end
    
    %create output structure
    out = struct('Real',tikzPointsReal,'Canv',tikzPointsCanv,'is3dAxis',axis3d);
end 

function [err,rngCropped] = createGraphicOfContent(fig,opts,crop)
    %actually creates the graphic file, which is then included via tikz
    %(png/eps/pdf)
    err = [];
    
    fileNameExtPic = fullfile(opts.filefolder,opts.filenamepic);
    rngCropped = [];
    if opts.useexportfig
        %use export_fig to create the graphic
        %----------------------------------------------------
        %cropping
        if numel(crop)==4
            cropString = ['-c',sprintf('%d ',round(crop,0))];
        else
            cropString = '-nocrop';
            crop = zeros(1,4);
        end
        
        %call export_fig with error handling
        try
            if strcmpi(opts.renderer,'-opengl')
                %opengl often fails to prodcue transparent background
                export_fig(fileNameExtPic,['-',opts.picformat],cropString,fig,opts.renderer)
            else
                export_fig(fileNameExtPic,['-',opts.picformat],'-transparent',cropString,fig,opts.renderer)
            end
        catch e
            err = e;
        end
        
        %determine cropped range
        fig = opts.figurehandle;
        fig.Units = 'pixels';
        pos = [1,1,fig.InnerPosition(3:4)];
        rngCropped = [pos(1)+crop(4), pos(2)+crop(3), pos(3)-crop(4)-crop(2), pos(4)-crop(3)-crop(1)];
    else
        %use custom routine
        %----------------------------------------------------
        if isempty(crop)
            crop = nan(1,4);
        end
        try
            resEF = exportFig(fileNameExtPic,opts.picformat,'transparent',true,opts.renderer,fig,...
                'crop',crop,'magnify',opts.magnify);
            rngCropped = resEF.posCroppedRange;
        catch e
            err = e;
        end
    end
end
%% 
%#######################################################################################################
% SUBFUNCTIONS - amend matlab2tikz file
%#######################################################################################################
function amendMatlab2tikzFile(opts,tikzPoints)    
    %get fileName and fullFileName
    fileName = opts.filename;
    fullFileName = fullfile(opts.filefolder,fileName);
    
    %add 3d graphics plot
    %--------------------------------------------------------------
    fileNamePic = fullfile(opts.relativeDataPath,opts.filenamepic);
    fileNamePic = replace(fileNamePic,filesep,'/');
    if ~strcmpi(opts.picformat,'tikz')
        if tikzPoints.is3dAxis
            lines = cell(3+size(tikzPoints.Real,2)+2,1);

            %add information to include graphic
            indx = 0;
            lines{indx+1} = '%include 3d plot as picture by tikz3dPlot';
            lines{indx+2} = '\addplot3 graphics[';%double backslash for fprintf
            lines{indx+3} = sprintf('\tpoints={%% important');
            indx = indx+3;
            for i=1:size(tikzPoints.Real,2)
                lines{indx+i} = sprintf('\t\t(%7.5g,%7.5g,%7.5g) => (%7.5g,%7.5g)',...
                    tikzPoints.Real(1,i),tikzPoints.Real(2,i),tikzPoints.Real(3,i),...
                    tikzPoints.Canv(1,i),tikzPoints.Canv(2,i));
            end
            indx = indx+size(tikzPoints.Real,2);
            lines{indx+1} = sprintf('\t}]');
            lines{indx+2} = sprintf('\t{%s};',fileNamePic);
        else
            lines = cell(4,1);
            
            %add information to include graphic
            lines{1} = '%include 2d plot as picture by tikz3dPlot';
            lines{2} = '\addplot graphics[';%double backslash for fprintf
            lines{3} = sprintf('\txmin=%.5g,xmax=%.5g,ymin=%.5g,ymax=%.5g]',...
                tikzPoints.Real(1,1),tikzPoints.Real(2,1),tikzPoints.Real(1,2),tikzPoints.Real(2,2));
            lines{4} = sprintf('\t{%s};',fileNamePic);
        end
        
        %add to extraPlotCodeEnd
        if isempty(opts.extraPlotCodeEnd)
            opts.extraPlotCodeEnd = lines;
        else
            opts.extraPlotCodeEnd = vertcat(lines,{''},opts.extraPlotCodeEnd);
        end
    end
    
    %read the current matlab2tikz file
    %--------------------------------------------------------------
    filecontent = fileread([fullFileName,'.tex']);
    linesOld = strsplit(filecontent,'\n','collapsedelimiters',false)';
    
    %find line with begining of plot commands
    [~,endAxisOpt] = regexp(filecontent,'\\begin\{axis\}\[%.*?\n\]\n');
    linePlotStart = count(filecontent(1:endAxisOpt),newline)+1;
    
    %find line with end of plot commands
    linePlotEnd = find(strcmpi(linesOld,'\end{axis}'))-1;
    
    %construct new content
    lines = linesOld(1:linePlotStart-1);
    
    %add additional code before
    if ~isempty(opts.extraPlotCode)
        lines = vertcat(lines,{''},opts.extraPlotCode,{''});
    end
    
    %add original plot content
    if linePlotEnd>=linePlotStart
        lines = vertcat(lines,linesOld(linePlotStart:linePlotEnd));
    end
    
    %add extra code after
    if ~isempty(opts.extraPlotCodeEnd)
        lines = vertcat(lines,{''},opts.extraPlotCodeEnd,{''});
    end
    
    %finish content
    lines = vertcat(lines,linesOld(linePlotEnd+1:end));

    %rewrite file
    filecontentnew = strjoin(lines,'\n');
    [fileID,errmsg] = fopen([fullFileName,'.tex'],'w');
    if fileID==-1
        error('tikz3dPlot:fileIO',...
            ['could not open the tikz file %s for new output. Error msg of System: ',errmsg],[fullFileName,'.tex'])
    end
    
    fprintf(fileID,'%s',filecontentnew);
    fclose(fileID);
end

%% 
%#######################################################################################################
% SUBFUNCTIONS - input arguments
%#######################################################################################################
function opts = inputParser(varargin)
    %default options
    opts = struct(...
        'filename',         'tikz3dPlot',...
        'filenamepic',      '',...              %set automatically
        'filefolder',       '',...
        'matlab2tikzopts',  {{}},...
        'renderer',         '-painters',...
        'useexportfig',     true,...
        'debug',            false,...
        'picformat',        'pdf',...
        'figurehandle',     [],...
        'relativeDataPath', '',...
        'standalone',       false,...           %needed to amend filename of export_fig
        'extraPlotCode',    {{}},...
        'extraPlotCodeEnd', {{}},...
        'magnify',          1,...
        'tikz3dPlotText',   ['tikz3dPlot',num2str(randi(1e7))]); %needed to allow to create a legend using tikz
    
    %list of all matlab2 tikz options (without filename and filehandle
    %since not possible to set) - these options are passed to matlab2tikz
    m2topts = {'figurehandle','colormap','strict','strictFontSize','showInfo','showWarnings',...
        'imagesAsPng','externalData','dataPath','relativeDataPath','height','width','noSize','extraCode',...
        'extraCodeAtEnd','extraAxisOptions','extraColors','extraTikzPictureOptions','encoding','floatformat',...
        'maxChunkLength','parseStrings','parseStringsAsMath','showHiddenStrings','interpretTickLabelsAsTex',...
        'arrowHeadSize','tikzFileComment','addLabels','standalone','checkForUpdates','semanticLineWidths'};
    
    ii = 1;
    setPicformat = false;
    while ii<=numel(varargin)
        if numel(varargin{ii})==1 && isgraphics(varargin{ii},'figure')
            varargin = horzcat(varargin,{'figurehandle',varargin{ii}});%#ok<AGROW>
        elseif ischar(varargin{ii}) && isvector(varargin{ii}) && ~isempty(varargin{ii})
            m2toptfound = false;
            %check for matlab2tikz
            if any(strcmpi(varargin{ii},m2topts))
                opts.matlab2tikzopts = horzcat(opts.matlab2tikzopts,varargin(ii:ii+1));
                m2toptfound = true;
            end

            %all other options
            switch lower(varargin{ii})
                case 'filename'
                    if ischar(varargin{ii+1}) && isvector(varargin{ii+1})
                        [opts.filefolder,opts.filename,~] = fileparts(varargin{ii+1});%do not store optional file suffix
                        ii = ii+1;
                    else
                        errorBadOption('option "filename" must be followed by a character vector')
                    end
                case {'-painters','-opengl'}
                    opts.renderer = varargin{ii};
                case 'debug'
                    opts.debug = true;
                case 'use_export_fig'
                    opts.useexportfig = true;
                case 'magnify'
                    if isnumeric(varargin{ii+1}) && isscalar(varargin{ii+1}) && varargin{ii+1}>0
                        opts.magnify = varargin{ii+1};
                        ii = ii+1;
                    else
                        errorBadOption('option "magnify" must be followed by a scalar > 0')
                    end
                case {'picformat','pdf','eps','png','-dpng','-dpdf','-depsc','-depsc2','tikz','-tikz'}
                    %warning for multiple settings of picformat
                    if setPicformat
                        warning('tikz3dPlot:DuplicatePicformat',['tikz3dPlot: Multiple picformats provided.',...
                            'Will use last one. Use either tikz3dPlot(...,picformat,format)',...
                            'or one of: pdf,eps,png,-dpng,-dpdf,-depsc,-depsc2,tikz,-tikz'])
                    end
                    %handle the input
                    if strcmpi(varargin{ii},'picformat')
                        if ischar(varargin{ii+1}) && isvector(varargin{ii+1}) && ~isempty(varargin{ii+1})
                            opt = varargin{ii+1};
                            ii = ii+1;
                        else
                            errorBadOption('option "picformat" must be followed by a character vector')
                        end
                    else
                        opt = varargin{ii};
                    end
                    %save final options
                    switch lower(opt)
                        case {'pdf','-dpdf'}
                            opts.picformat = 'pdf';
                        case {'eps','-depsc','-depsc2'}
                            opts.picformat = 'eps';
                        case {'png','-dpng'}
                            opts.picformat = 'png';
                        case {'tikz','-tikz'}
                            opts.picformat = 'tikz';
                        otherwise 
                            errorBadOption('picformat not implemented')
                    end
                    setPicformat = true;
                case 'figurehandle'
                    if isgraphics(varargin{ii+1},'figure')
                        opts.figurehandle = varargin{ii+1};
                        ii = ii+1;
                    else
                        errorBadOption('option "figurehandle" must be followed by a non deleted figure handle. ')
                    end
                case 'relativedatapath'
                    if ischar(varargin{ii+1}) && isvector(varargin{ii+1})
                        opts.relativeDataPath = varargin{ii+1};
                        ii = ii+1;
                    else
                        errorBadOption('option "relativDataPath" must be followed by a character vector')
                    end
                case 'standalone'
                    if islogical(varargin{ii+1}) && isscalar(varargin{ii+1})
                        opts.standalone = varargin{ii+1};
                        ii = ii+1;
                    else
                        errorBadOption('option "standalone" must be followed by true/false')
                    end
                case {'extraplotcode','extraplotcodeend'}
                    if ischar(varargin{ii+1}) && isvector(varargin{ii+1})
                        zw = varargin(ii+1);
                    elseif iscell(varargin{ii+1}) 
                        zw = varargin{ii+1};
                        zw = zw(:);
                        if ~all(cellfun(@ischar,zw)) || ~all(cellfun(@isvector,zw))
                            errorBadOption(['options "extraPlotCode" and "extraPlotCodeEnd" must be either ',...
                                'followed by a character vector or a cell of character vectors'])
                        end
                    else
                        errorBadOption(['options "extraPlotCode" and "extraPlotCodeEnd" must be either ',...
                            'followed by a character vector or a cell of character vectors'])
                    end
                    opts.(varargin{ii}) = zw;
                    ii = ii+1;
                otherwise
                    if m2toptfound
                        ii = ii+1;
                    elseif strcmpi(opts.filename,'tikz3dPlot')
                        %set the filename to be found automatically 
                        opts.filename = varargin{ii};
                        varargin = horzcat(varargin,'filename',varargin{ii});%#ok<AGROW>
                    else
                        errorBadOption(['Unrecognized tikz3dPlot input option: ''' varargin{ii} '''']);
                    end
            end
        else
            %completely unknown option
            disp(varargin{ii})
            errorBadOption('Unrecognized tikz3dPlot input option (displayed above)');
        end
        ii = ii+1;
    end
    
    %try to get standard input for non set options
    if isempty(opts.figurehandle)
        opts.figurehandle = get(0, 'CurrentFigure');
        if isempty(opts.figurehandle)
            error('tikz3dPlot:BadOption','no open figure and no figurehandle provided')
        end
        %add also to matlab2tikz options
        opts.matlab2tikzopts = horzcat(opts.matlab2tikzopts,{'figurehandle',opts.figurehandle});
    end
    
    %check if the folder exists
    if ~isempty(opts.filefolder) && ~(exist(opts.filefolder,'dir')==7)
        error('tikz3dPlot:fileIO','tikz3dPlot: provided folder "%s" for file not found',opts.filefolder)
    elseif isempty(opts.filefolder)
        opts.filefolder=pwd;%set folder to current path
    end
    
    %check for number of axes on figure
    count = 0;
    for i=1:numel(opts.figurehandle.Children)
        if isa(opts.figurehandle.Children(i),'matlab.graphics.axis.Axes')
            count = count + 1;
        end
    end
    if count~=1
        warning('tikz3dPlot:numberOfAxes',['tikz3dPlot: figure should have exactly one axes object. ',...
            'unrespected results possible'])
    end
    
    %check if figure is docked
    if ~strcmpi(opts.figurehandle.WindowStyle,'normal')
        errorBadOption('Figures cannot be "docked" for tikz3dPlot to work')
    end
    
    %create the filename of the picture
    if isempty(opts.filenamepic)
        if opts.standalone
            opts.filenamepic = [opts.filename,'_pic'];
        else
            opts.filenamepic = opts.filename;
        end
        opts.filenamepic = fullfile(opts.filenamepic);
    end
end
function errorBadOption(message)
    error('tikz3dPlot:BadOption',['tikz3dPlot: ',message]);
end

%%
%#######################################################################################################
% SUBFUNCTIONS - tikz points
%#######################################################################################################
% defines points in 3d and finds corresponding canvas points

% find tikz points
function [pointsReal,pointsCanv]=supportPointsTikz(axes,croppedRange)
    %computes the support points in the physical space
    
    %get the size of the canvas
    fig = axes.Parent;
    zw = fig.Units;
    fig.Units = 'pixel';
    canvas = fig.InnerPosition;
    canvas(1:2) = 1;
    fig.Units = zw;
    
    %box points
    %-----------------------------------------------
    boxReal = plotBoxPointsOnCanvas(axes,canvas);
    boxCanv = zeros(2,size(boxReal,2));
    for i=1:size(boxReal,2)
        boxCanv(:,i) = projectionMatlab(boxReal(:,i),axes);
    end
    
    %exclude all points that (still - after plotBoxPointsOnCanvas()) lie
    %outside the canvas
    %-----------------------------------------------
    indx = checkPointsOnCanvas(boxCanv,canvas);
    boxReal = boxReal(:,indx);
    boxCanv = boxCanv(:,indx);
    
    if size(boxReal,2)<4 
        error('tikz3dPlot:findTikzPoints','could not find enough points for tikz 3d')
    end
    
    %exclude all points that lie above of each other
    %-----------------------------------------------
    [~,indx,~] = unique(boxCanv','rows');
    boxReal = boxReal(:,indx);
    boxCanv = boxCanv(:,indx);
    
    if size(boxReal,2)<4 
        error('tikz3dPlot:findTikzPoints','could not find enough unique points for tikz 3d')
    end
    
    %find the best suitable permutaion of points
    %-----------------------------------------------
    combination = nchoosek(1:numel(indx),4);
    indx = 0;
    condNr = 0;
    for i=1:size(combination,1)
        pReal = boxReal(:,combination(i,:));
        pCanv = boxCanv(:,combination(i,:));
        
        %linear system for tikz points (internal in tikz)
        [A,~] = linearSystemTikZpoints(pReal,pCanv);
        
        %find best defined system
        zw = rcond(A);
        if zw>condNr && det(A)~=0
            condNr = zw;
            indx = i;
        end
    end
    
    %select best set
    boxReal = boxReal(:,combination(indx,:));
    boxCanv = boxCanv(:,combination(indx,:));
    
    %take points that have sufficient differnce in canvas
    %coordinates as first points
    %-----------------------------------------------
    permutations = sortrows(perms(1:4));
    indx = 0;
    fact = 0.5;
    while indx==0
        fact = fact*0.8;
        for i=1:size(permutations,1)
            ind1 = permutations(i,1);
            ind2 = permutations(i,2);
            
            %check if points have sufficient distance
            if abs(boxCanv(1,ind1)-boxCanv(1,ind2))>canvas(3)*fact ...
            && abs(boxCanv(2,ind1)-boxCanv(2,ind2))>canvas(4)*fact
                indx = i;
                break
            end        
        end
    end
    
    if fact<0.1
        warning('tikz points might be poor')
    end
    
    pointsReal = boxReal(:,permutations(indx,:));
    pointsCanv = boxCanv(:,permutations(indx,:));
    
    %adjust canvas coordinates for cropping
    %-----------------------------------------------
    pointsCanv(1,:) = pointsCanv(1,:)-croppedRange(1)+1;
    pointsCanv(2,:) = pointsCanv(2,:)-croppedRange(2)+1;
    
    %rescale canvas points with tex points
    %-----------------------------------------------
    dpi = get(0,'ScreenPixelsPerInch');
    pointsCanv = pointsCanv * 72.2/dpi;
end
function [A,b] = linearSystemTikZpoints(pReal,pCanv)
    A = zeros(8,8);
    b = zeros(8,1);
    
    for i=1:4
        indx1 = (i-1)*2+1;
        indx2 = indx1+1;
        
        row = [1,pReal(:,i)'];
        A(indx1,1:2:end) = row;
        A(indx2,2:2:end) = row;
        
        b([indx1,indx2]) = pCanv(:,i);
    end
end
function points=plotBoxPointsOnCanvas(axes,canvas)
    %searches for a plot box that lies on the canvas

    %get the physical and canvas plot points of the axis limits
    boxRealStart = plotBoxPoints(axes);
    meanBoxReal = mean(boxRealStart,2);
    boxCanv = zeros(2,size(boxRealStart,2));
    
    %search a fitting plotbox
    done = false;
    scal = 1; %scaling factor
    boxReal = boxRealStart;
    while ~done
        %scale the current box
        for i=1:size(boxReal,2)
            boxReal(:,i) = meanBoxReal + scal*(boxRealStart(:,i)-meanBoxReal);
        end
        
        %get canvas points
        for i=1:size(boxReal,2)
            boxCanv(:,i) = round(projectionMatlab(boxReal(:,i),axes),0);
        end
        
        %check if canvas points 
        [~,violation] = checkPointsOnCanvas(boxCanv,canvas);
        if violation>1
            scal = scal/violation/1.02;%saftey margin of 2
        else
            done = true;
        end
    end
    
    points = boxReal;
end
function points=plotBoxPoints(axes)
    %returns corner points of the 3d box
    xlim = axes.XLim;
    ylim = axes.YLim;
    zlim = axes.ZLim;
    
    points = [...
        xlim([1,2,1,2,1,2,1,2]);
        ylim([1,1,2,2,1,1,2,2]);
        zlim([1,1,1,1,2,2,2,2])];
end
function [indx,violation] = checkPointsOnCanvas(points,canvas)
    %find points that are not on canvas
    indx = ~(points(1,:)<canvas(1) | points(2,:)<canvas(2)...
        | points(2,:)>canvas(3) | points(2,:)>canvas(4));
    
    %compute the violation
    if ~all(indx)
        violation = 0;
        violation = max(violation, max(abs(2*points(1,:)-canvas(3)+canvas(1)))/(canvas(3)-canvas(1)));
        violation = max(violation, max(abs(2*points(1,:)-canvas(4)+canvas(2)))/(canvas(4)-canvas(2)));
    else 
        violation = 0;
    end
end 

%% 
%#######################################################################################################
% SUBFUNCTIONS - projection
%#######################################################################################################
function canvasPoint = projectionMatlab(point,axes,debug)
    %conducts the projection of a 3D point as defined in matlab. 
    %Described in depth in:
    %https://de.mathworks.com/matlabcentral/answers/373545-how-can-i-convert-from-the-pixel-position-in-an-image-to-3d-world-coordinates
    
    %https://de.mathworks.com/matlabcentral/fileexchange/43896-3d-data-space-coordinates-to-normalized-figure-coordinates-conversion
    if nargin<=2
        debug=false;
    end
    
    if ~strcmpi(axes.Projection,'orthographic')
        error('tikz3dPlot:findTikzPoints',[...
        'currently only orthographic projections possible in pgfplots.',...
        'See pgfplots manual 1.16 page 73'])
    end
    
    %4d coordinates needed for transformation
    xHom = [point;ones(size(point,2))];
    
    % obtain data needed
    %--------------------------------------
    % camera
    viewAngle = axes.CameraViewAngle;
    viewTarget = axes.CameraTarget;
    viewPosition = axes.CameraPosition;
    viewUp = axes.CameraUpVector;
    % axes direction
    axesDirection = strcmpi({axes.XDir,axes.YDir,axes.ZDir}, 'normal');
    % data scale
    dataZLim = axes.ZLim;
    dataRatio = axes.DataAspectRatio;
    plotBoxRatio = axes.PlotBoxAspectRatio;
    % is perspective
    isPerspective = strcmp(axes.Projection, 'perspective');
    % axes position
    axesUnitsOriginal = axes.Units;
    axes.Units = 'normalized';
    positionNormal = axes.Position;
    axes.Units = 'pixel';
    positionPixel = axes.Position;
    axes.Units = axesUnitsOriginal;
    % stretch
    stretchMode = strcmpi({axes.CameraViewAngleMode, ...
        axes.DataAspectRatioMode, axes.PlotBoxAspectRatioMode}, 'auto');
    stretchToFill = all(stretchMode);
    stretchToFit = ~stretchToFill && stretchMode(1);
    stretchNone = ~stretchToFill && ~stretchToFit;
    
    % model view matrix
    %--------------------------------------
    % move data space center to viewTarget point
    matrixTranslate = eye(4);
    matrixTranslate(1:3, 4) = -viewTarget;
    % rescale data
    % note: matlab will rescale data space by dividing DataAspectRatio
    %       and normalize z direction to 1 to makeup the 'PlotBox'
    if stretchToFill
        %https://de.mathworks.com/matlabcentral/answers/326180-how-to-get-
        %accurate-dataaspectratio-when-stretch-to-fit-is-enabled
        dataRatioReal = dataRatio./plotBoxRatio;
        error('stretch to fill behavior does not work yet')
    else
        dataRatioReal = dataRatio;
    end
    scaleFactor = dataRatioReal / dataRatioReal(3) * (dataZLim(2) - dataZLim(1));
    scaleDirection = axesDirection * 2 - 1;
    matrixRescale = diag([scaleDirection ./ scaleFactor, 1]);
    % rotate the 'PlotBox'
    vecticesZUp = ...
        [matrixRescale*matrixTranslate * [viewPosition, 1]', matrixRescale*[viewUp, 1]'];
    zVector = vecticesZUp(1:3, 1);
    upVector = vecticesZUp(1:3, 2);
    viewDistance = norm(zVector);
    zDirection = zVector / viewDistance;
    yVector = upVector - zDirection * dot(zDirection, upVector);
    yDirection = yVector / norm(yVector);
    matrixRotate = blkdiag( ...
        [cross(yDirection, zDirection), yDirection, zDirection]', 1);
    
    % projection matrix
    %--------------------------------------
    % note: matlab will project the rotated 'PlotBox' to an area of 
    %       [-0.5, 0.5; -0.5, 0.5]
    matrixProjection = eye(4);
    matrixProjection(4, 3) = -isPerspective / viewDistance;
    projectionArea = 2 * tan(viewAngle * pi / 360) * viewDistance;
    matrixProjection = diag([ones(1, 3), projectionArea]) * matrixProjection;
    
    % stretch matrix
    %--------------------------------------
    % stretch the projective 'PlotBox' into the position retangle of the axes
    % note: stretch will first detect data region
    if stretchToFill || stretchToFit
        plotBox = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1]' - .5;
        plotBox = diag(plotBoxRatio / plotBoxRatio(3)) * plotBox;
        edgeVertices = matrixProjection * matrixRotate * [plotBox; ones(1, 8)];
        edgeVertices(1, :) = edgeVertices(1, :) ./ edgeVertices(4, :);
        edgeVertices(2, :) = edgeVertices(2, :) ./ edgeVertices(4, :);
        edgeVertices = edgeVertices(1:2, :)';
        % note: the low boundary and the high boundary of data region may be
        %       difference in perspective projection, so the figure should move
        %       to center, but here no need to do so, because matlab ignore it
        dataRegion = max(edgeVertices) - min(edgeVertices);
        % note: matlab have a strange addition stretch in stretch to fit mode.
        %       one side of the data region will hit the position rectangle,
        %       and matlab will assume data region of that side to be 1 keeping
        %       aspect ratio.
        strangeFactor = dataRegion ./ positionPixel(3:4);
        if stretchToFit %|| stretchToFill
            if strangeFactor(1) > strangeFactor(2)
                dataRegionMod = dataRegion / dataRegion(1);
            else
                dataRegionMod = dataRegion / dataRegion(2);
            end
        elseif stretchToFill
            dataRegionMod = dataRegion;
        end
    else
        % note: if no stretch, it will use projection area as data region
        dataRegionMod = [1,1];
    end
    % note: stretch than apply a stretchFactor to the data, such that it fit
    %       in the axes position retangle
    if stretchToFit || stretchNone
        stretchFactor = dataRegionMod ./ positionPixel(3:4);
        stretchFactor = stretchFactor ./ max(stretchFactor);
    else
        stretchFactor = [1,1];
    end
    matrixStretch = diag([stretchFactor ./ dataRegionMod, 1, 1]);
    
    % view port matrix
    %--------------------------------------
    % here change to pixel coordinates in comparison to ds2fig
%     matrixViewPortNorm = diag([positionNormal(3:4), 1, 1]);
%     matrixViewPortNorm(1:2, 4) = positionNormal(1:2) + positionNormal(3:4)/2;
    matrixViewPortPix = diag([positionPixel(3:4),1,1]);
    matrixViewPortPix(1:2, 4) = positionPixel(1:2) + positionPixel(3:4)/2;
    
    % return transformation matrix
    %--------------------------------------
    matrixTransform = matrixViewPortPix * matrixStretch * matrixProjection * ...
        matrixRotate * matrixRescale * matrixTranslate;
    
    % actual computation of the point on the canvas
    %--------------------------------------
    %normalized figure coordinates
    canvCoordNorm = matrixTransform*xHom;
    canvCoordNorm = canvCoordNorm(1:3)./canvCoordNorm(4);
    
    canvasPoint = canvCoordNorm(1:2);
    
    if debug==true
        fig = axes.Parent;
        drawAxesPosAnnotation(fig,axes,{'color','r','linewidth',2});
        
        %original points
        plot3(axes,point(1),point(2),point(3),...
            'Marker','x', 'MarkerFaceColor','none', 'MarkerEdgeColor','r', ...
            'LineStyle','none', 'MarkerSize',40);
        boxPoints = plotBoxPoints(axes);
        plot3(axes,boxPoints(1,:),boxPoints(2,:),boxPoints(3,:),...
            'Marker','o', 'MarkerFaceColor','none', 'MarkerEdgeColor','k', ...
            'LineStyle','none', 'MarkerSize',10);
        
        %camara line of sight
        plot3(axes,[viewPosition(1),viewTarget(1)],[viewPosition(2),viewTarget(2)],[viewPosition(3),viewTarget(3)],...
            'LineStyle','-','Color','c','MarkerFaceColor','c','MarkerSize',10,'Marker','o');
        
        %point after translation, rescale (model transform)
        
        
        %add axes with figure canvas
        axes3 = createHelpAxis(fig,axes);
        axes3.InnerPosition = [0,0,1,1];
        fig.Units = 'pixel';
        axes3.XLim = [1,fig.InnerPosition(3)];
        axes3.YLim = [1,fig.InnerPosition(4)];
        axes3.ZLim = [0,1];
        plotPointAndBorders(fig,axes3,canvCoordNorm(1:2),{'color','m','linestyle','-'},'*','m')
        
        %reset current axes
        fig.CurrentAxes = axes;
    end
end
% function canvasPoint = projectionMatlabOld(point,axes,debug)
%     if nargin<=2
%         debug=false;
%     end
%     %description of parallel projection
%     %https://en.wikipedia.org/wiki/Parallel_projection
%     if ~strcmpi(axes.Projection,'orthographic')
%         error('tikz3dPlot:findTikzPoints',[...
%         'currently only orthographic projections possible in pgfplots.',...
%         'See pgfplots manual 1.16 page 73'])
%     end
%         
%     %get parameters of projection
%     %--------------------------------------
%     %camera information
%     cp = axes.CameraPosition;
%     ct = axes.CameraTarget;
%     cuv = axes.CameraUpVector;
%     cva = axes.CameraViewAngle;
%     da = axes.DataAspectRatio;
%     
%     %direction of view (computed from camera not view!)
%     v = (cp-ct)'/norm(cp-ct);
%     %view = axes.View;
%     %v = [sind(view(1)),-cosd(view(1)),sind(view(2))]';
%     %v = v/norm(v);
%     
%     %projection plane (normal n and position d)
%     n = (cp-ct)'/norm(cp-ct);
%     d = n'*ct';
%     
%     %projection of point to canvas
%     %--------------------------------------
%     pointCv = point + ((d-point'*n)/(n'*v))*v;
%     if debug==true && abs(n'*pointCv-d)>100*eps
%         error('projected point not on canvas')
%     end
%     
%     %local coordinate system on canvas
%     %--------------------------------------
%     Eold = eye(3);
%     
%     %new coordinate system
%     Enew = zeros(3,3);
%     Enew(:,1) = cross(cuv,n)/norm(cross(cuv,n));
%     Enew(:,2) = cross(n,Enew(:,1));
%     Enew(:,3) = n;
%     
%     %rotation matrix
%     rotMat = Eold'*Enew;
%     if debug==true && abs(abs(det(rotMat))-1)>100*eps
%         error('no rotation matrix')
%     end
%     
%     %actual rotation
%     pointCvLoc = rotMat'*(pointCv-ct');
%     if debug==true && abs(pointCvLoc(3))>100*eps
%         error('transformation onto canvas faulty')
%     end
%     
%     %compute point in normalized coord sys
%     %--------------------------------------
%     %compute width of canvas with camera view angle and distance
%     %https://de.mathworks.com/help/matlab/creating_plots/defining-scenes-with-camera-graphics.html
%     
%     %transform all edge points of the 3d box to get the dimensions of the
%     %enclosing rectangle
%     %box points 
%     boxPoints = plotBoxPoints(axes);
%     boxPointsCv = boxPoints + repmat(((d-boxPoints'*n)/(n'*v))', 3,1).*repmat(v,1,size(boxPoints,2));
%     boxPointsCvLoc = rotMat'*(boxPointsCv-repmat(ct',1,size(boxPoints,2)));
%     sizeBox = [max(boxPointsCvLoc(1,:))-min(boxPointsCvLoc(1,:)),...
%                max(boxPointsCvLoc(2,:))-min(boxPointsCvLoc(2,:))];
%     
%     %get physical and canvas size of canvas plane
%     axes.Units = 'pixel';
%     widthReal = norm(cp-ct)*tand(cva/2)*2;
%     if abs(widthReal-sizeBox(1))/sizeBox(1)<1e-10
%         widthCanv = axes.Position(3);
%     elseif abs(widthReal-sizeBox(2))/sizeBox(2)<1e-10
%         widthCanv = axes.Position(4);
%     else
%         widthCanv = min(axes.Position(3:4));
%     end
%     
%     %perform transformation
%     scalfact = widthCanv/widthReal;
%     pointCvLocNorm = pointCvLoc(1:2)*scalfact./axes.Position(3:4)'+0.5;
%     
%     %compute point in pixels on figure coord sys
%     %--------------------------------------
%     axes.Units = 'pixel';
%     posPix = axes.InnerPosition;
%     axes.Units = 'normalized';
%     
%     pointCvPix = pointCvLocNorm.*(posPix(3:4)'-1) + posPix(1:2)';
%     
%     canvasPoint = pointCvPix;
%     
%     %plotting if requested
%     %--------------------------------------
%     if debug==true
%         fig = axes.Parent;
%         drawAxesPosAnnotation(fig,axes,{'color','r','linewidth',2});
%         
%         %original point
%         plot3(axes,point(1),point(2),point(3),...
%             'Marker','x', 'MarkerFaceColor','none', 'MarkerEdgeColor','r', ...
%             'LineStyle','none', 'MarkerSize',40);
%         plot3(axes,boxPoints(1,:),boxPoints(2,:),boxPoints(3,:),...
%             'Marker','o', 'MarkerFaceColor','none', 'MarkerEdgeColor','k', ...
%             'LineStyle','none', 'MarkerSize',10);
%         %on canvas
%         plot3(axes,pointCv(1),pointCv(2),pointCv(3),...
%             'Marker','square', 'MarkerFaceColor','none', 'MarkerEdgeColor','g', ...
%             'LineStyle','none', 'MarkerSize',40);
%         plot3(axes,boxPointsCv(1,:),boxPointsCv(2,:),boxPointsCv(3,:),...
%             'Marker','+', 'MarkerFaceColor','none', 'MarkerEdgeColor','k', ...
%             'LineStyle','none', 'Tag','BoxPoints', 'MarkerSize',10);
%         
%         %camara line of sight
%         plot3(axes,[cp(1),ct(1)],[cp(2),ct(2)],[cp(3),ct(3)],...
%             'LineStyle','-','Color','c','MarkerFaceColor','c','MarkerSize',10,'Marker','o');
%         
%         %axis cross on canvas
%         x = [ct'+Enew(:,1)*widthReal/2, ct'-Enew(:,1)*widthReal/2];
%         plot3(axes,x(1,:),x(2,:),x(3,:),...
%             'Marker','o', 'MarkerFaceColor','k', 'MarkerEdgeColor','k', ...
%             'LineStyle','-', 'color','k', 'MarkerSize',10);
%         x = [ct'+Enew(:,2)*widthReal/2, ct'-Enew(:,2)*widthReal/2];
%         plot3(axes,x(1,:),x(2,:),x(3,:),...
%             'Marker','o', 'MarkerFaceColor','k', 'MarkerEdgeColor','k', ...
%             'LineStyle','-', 'color','k', 'MarkerSize',10);
%         
%         %add axes on canvas without normalization 
%         axes1 = createHelpAxis(fig,axes);
%         axes1.XLim = [-1,1]*widthReal/2*axes.Position(3)/widthCanv;
%         axes1.YLim = [-1,1]*widthReal/2*axes.Position(4)/widthCanv;
%         axes1.ZLim = [0,1];
%         plotPointAndBorders(fig,axes1,pointCvLoc,{'color','c','linestyle','-'},'+','c')
%         
%         %add axes with normalized canvas
%         axes2 = createHelpAxis(fig,axes);
%         axes2.XLim = [0,1];
%         axes2.YLim = [0,1];
%         axes2.ZLim = [0,1];
%         plotPointAndBorders(fig,axes2,pointCvLocNorm,{'color','b','linestyle','-'},'o','b')
%         
%         %add axes with figure canvas
%         axes3 = createHelpAxis(fig,axes);
%         axes3.InnerPosition = [0,0,1,1];
%         fig.Units = 'pixel';
%         axes3.XLim = [1,fig.InnerPosition(3)];
%         axes3.YLim = [1,fig.InnerPosition(4)];
%         axes3.ZLim = [0,1];
%         plotPointAndBorders(fig,axes3,pointCvPix,{'color','m','linestyle','-'},'*','m')
%         
%         %reset current axes
%         fig.CurrentAxes = axes;
%     end
% end
%help functions for debbung of projection routine
function axes1 = createHelpAxis(fig,axes)
    %creates a new axes and prepares it for output
    axes1 = copyobj(axes,fig);
    reset(axes1)
    for i=1:numel(axes1.Children)
        delete(axes1.Children(1))
    end
    hold(axes1,'on')
end
function plotPointAndBorders(fig,axes,point,formatAnnotation,pointMarker,pointColor)
    %only for debugging
    %draws the point point in the axes and formats it, draws a box around
    %the limits of the axes and draws the annotation for the position of
    %the axes
    plot(axes,point(1),point(2),...
        'Marker',pointMarker, 'MarkerFaceColor','none', 'MarkerEdgeColor',pointColor, ...
        'LineStyle','none', 'Tag','BoxPoints', 'MarkerSize',40);
    drawAxesbox(axes,{'d-','color',pointColor})
    drawAxesPosAnnotation(fig,axes,formatAnnotation)
    axes.Visible = 'off';
end
function drawAxesPosAnnotation(fig,axes,format)
    %draws a rectangle at the position of the axes (only for debugging)
    %normalized units needed (without overwriting)
    zw = axes.Units;
    axes.Units = 'normalized';
    pos = axes.Position;
    axes.Units = zw;
    
    %check if annotation is completely on figure
    if pos(1)>=0 && pos(1)<=1 && pos(2)>=0 && pos(2)<=1 && pos(1)+pos(3)<=1 && pos(2)+pos(4)<=1
        %normalized units needed (without overwriting)
        zw = fig.Units;
        fig.Units = 'normalized';
        annotation(fig,'rectangle',pos,format{:})
        fig.Units = zw;
    end
end
function drawAxesbox(axes,style)
    plot(axes,...
        [axes.XLim(1),axes.XLim(2),axes.XLim(2),axes.XLim(1),axes.XLim(1)],...
        [axes.YLim(1),axes.YLim(1),axes.YLim(2),axes.YLim(2),axes.YLim(1)],style{:});
end


%%
%#######################################################################################################
% SUBFUNCTIONS - creating the legend using tikz
%#######################################################################################################
function legendTemp = makeTikzLegend(figHandler,randText)
    %function that adds lines/patches outside of the axis in the correct
    %format so that the legend can be created using matlab to tikz. 
    %(if everything is exported using a graphic, matlab2tikz cannot create
    %a legend)
    
    legen = figHandler.legen;
    axes = figHandler.axes;
    
    legendTemp = [];
    if ~isempty(legen) && strcmpi(legen.Visible,'on')
        numObjs = numel(axes.Children);

        %switch on update of the legend
        zwUpdate = legen.AutoUpdate;
        legen.AutoUpdate = 'off';

        %array for new graphics objects
        newObjs = cell(1,numel(legen.String));

        %copy all visible objects with a legend entry
        ii = 0;
        count = 0;
        while ii<numObjs
            legendObj = axes.Children(end-ii);%bottom to top to have the correct legend order
            if ~isempty(legendObj.DisplayName) && strcmpi(legendObj.Visible,'on') %has a legend entry
                graphObj = [];
                switch lower(class(legendObj))
                    case lower('matlab.graphics.chart.primitive.Line')
                        graphObj = addLinePlotFromHandleForLegend(axes,legendObj);
                    case lower('matlab.graphics.primitive.Patch')
                        graphObj = addPatchPlotFromHandleForLegend(axes,legendObj);
                end
                count = count+1;
                newObjs{count} = graphObj;
                
                % add the random string to identify objects added by tikz3dPlot
                if ~isempty(graphObj)
                    graphObj.UserData = randText;
                end
            end
            ii = ii+1;
        end

        %resotre update of the legend
        legen.AutoUpdate = zwUpdate;
        
        %create a new legend and switch off the old one
        legendTemp = createLegend(axes,[newObjs{:}],figHandler.oldLegendProps);
        legendTemp.AutoUpdate = 'off';
    end
end
function switchONTikzLegendObjects(axes,randText)
    %switches the visibility of the legend objects on
    for i=1:numel(axes.Children)
        if ischar(axes.Children(i).UserData) && isvector(axes.Children(i).UserData)...
        && strcmpi(axes.Children(i).UserData,randText)
            axes.Children(i).Visible = 'on';
        end
    end
end
function deleteTikzLegendObjects(figHandler,randText)
    axes = figHandler.axes;
    %deletes all temporary tikz legend objects
    
    i = 1;
    while i<=numel(axes.Children)
        if ischar(axes.Children(i).UserData) && isvector(axes.Children(i).UserData)...
        && strcmpi(axes.Children(i).UserData,randText)
            delete(axes.Children(i));
        else
            i = i+1;
        end
    end
    
    %resotre the legend
    if ~isempty(figHandler.legen)
        legendObjs = axes.Children(~cellfun(@isempty,figHandler.oldAxesProps.childDisplayNames));
        createLegend(axes,legendObjs,figHandler.oldLegendProps);
    end
end
function out = createLegend(axes,graphObjs,props)
    out = legend(axes,graphObjs,'AutoUpdate',props.AutoUpdate,'Color',props.Color,'Box',props.Box,...
        'TextColor',props.TextColor,'LineWidth',props.LineWidth,'EdgeColor',props.EdgeColor,...
        'FontSize',props.FontSize,'FontWeight',props.FontWeight,'FontName',props.FontName,...
        'FontAngle',props.FontAngle,'NumColumns',props.NumColumns,...
        'NumColumnsMode',props.NumColumnsMode,'Visible',props.Visible,'Interpreter',props.Interpreter,...
        'Position',props.Position,'Units',props.Units,'Orientation',props.Orientation,...
        'Location',props.Location);
    %add a title to the legend
    if ~isempty(props.TitleStruct.String)
        proptit = props.TitleStruct;
        title(out,proptit.String,'FontSize',proptit.FontSize,'FontSizeMode',proptit.FontSizeMode,...
            'FontAngle',proptit.FontAngle,'FontAngleMode',proptit.FontAngleMode,...
            'FontName',proptit.FontAngle,'FontNameMode',proptit.FontAngleMode,...
            'FontWeight',proptit.FontWeight,'FontWeightMode',proptit.FontWeightMode,...
            'Color',proptit.Color,'ColorMode',proptit.ColorMode,...
            'Interpreter',proptit.Interpreter,'InterpreterMode',proptit.InterpreterMode,...
            'Visible',proptit.Visible);
    end
end
function out = addLinePlotFromHandleForLegend(axes,lineObj)
    %adds a line plot outside of the axis to allow the legend with tikz
    out = plot(axes,axes.XLim(2)+abs(axes.XLim(2)),axes.YLim(2)+abs(axes.YLim(2)),...
        'color',lineObj.Color,'Marker',lineObj.Marker,...
        'MarkerEdgeColor',lineObj.MarkerEdgeColor,'MarkerFaceColor',lineObj.MarkerFaceColor,...
        'MarkerSize',lineObj.MarkerSize,'LineStyle',lineObj.LineStyle,'LineWidth',lineObj.LineWidth,...
        'DisplayName',lineObj.DisplayName);
end
function out = addPatchPlotFromHandleForLegend(axes,patchObj)
    %adds a patch object outside of the axes to allow the legend with tikz
    out = patch(axes,'Vertices',patchObj.Vertices,'faces',patchObj.Faces,...
        'EdgeColor',patchObj.EdgeColor,'EdgeAlpha',patchObj.EdgeAlpha,'facecolor',patchObj.FaceColor,...
        'FaceAlpha',patchObj.FaceAlpha,'Marker',patchObj.Marker,'MarkerEdgeColor',patchObj.MarkerEdgeColor,...
        'MarkerFaceColor',patchObj.MarkerFaceColor,'MarkerSize',patchObj.MarkerSize,'LineStyle',patchObj.LineStyle,...
        'LineWidth',patchObj.LineWidth,'DisplayName',patchObj.DisplayName);
    %move it outside of the axes
    if size(out.Vertices,2)>=1
        out.Vertices(:,1) = axes.XLim(2)+abs(axes.XLim(2));
    end
    if size(out.Vertices,2)>=2
        out.Vertices(:,2) = axes.YLim(2)+abs(axes.YLim(2));
    end
    if size(out.Vertices,2)>=3
        out.Vertices(:,3) = axes.ZLim(2)+abs(axes.ZLim(2));
    end
end

%%
%#######################################################################################################
% SUBFUNCTIONS - plot handels
%#######################################################################################################

function bool = isAxis3D(axisHandle)
    % this function is taken from matlab2tikz
    % Check if elevation is not orthogonal to xy plane
    axisView = get(axisHandle,'view');
    bool     = ~ismember(axisView(2),[90,-90]);
end
function switchOFFgraphicsObjects(obj)
    %delete all graphics objects from plots
    for i=1:numel(obj.Children)
        obj.Children(i).Visible = 'off';
    end
end
function switchONgraphicsObjects(obj,visibility)
    %restore visibility
    for i=1:numel(obj.Children)
        obj.Children(i).Visible = visibility(i,:);
    end
end

%% 
%#######################################################################################################
% SUBFUNCTIONS - external
%#######################################################################################################
% %https://de.mathworks.com/matlabcentral/answers/363670-is-there-a-way-to-determine-illegal-characters-for-file-names-based-on-the-computer-operating-system
% function [Valid, Msg] = CheckFileName(S)
%     Msg = '';
%     if ispc
%         BadChar = '<>:"/\|?*';
%         BadName = {'CON', 'PRN', 'AUX', 'NUL', 'CLOCK$', ...
%             'COM1', 'COM2', 'COM3', 'COM4', 'COM5', 'COM6', ...
%             'COM7', 'COM8', 'COM9', ...
%             'LPT1', 'LPT2', 'LPT3', 'LPT4', 'LPT5', 'LPT6', ...
%             'LPT7', 'LPT8', 'LPT9'};
%         bad = ismember(BadChar, S);
%         if any(bad)
%             Msg = ['Name contains bad characters: ', BadChar(bad)];
%         elseif any(S < 32)
%             Msg = ['Name contains non printable characters, ASCII:', sprintf(' %d', S(S < 32))];
%         elseif ~isempty(S) && (S(end) == ' ' || S(end) == '.')
%             Msg = 'A trailing space or dot is forbidden';
%         else
%             % "AUX.txt" fails also, so extract the file name only:
%             [~, name] = fileparts(S);
%             if any(strcmpi(name, BadName))
%                 Msg = ['Name not allowed: ', name];
%             end
%         end
%     else  % Mac and Linux:
%         if any(S == '/')
%             Msg = '/ is forbidden in a file name';
%         elseif any(S == 0)
%             Msg = '\0 is forbidden in a file name';
%         end
%     end
%     Valid = isempty(Msg);
% end