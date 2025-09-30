function out = exportFig(varargin)
    %EXPORTFIG creates a eps, pdf or png of a figure. 
    %
    %Based on export_fig
    %https://de.mathworks.com/matlabcentral/fileexchange/23629-export_fig
    %
    %EXPORTFIG options can be combined at will
    %
    %EXPORTFIG(FILENAME,...) or EXPORTFIG('filename',FILENAME,...) stores the
    %figure in the format provided under 'picformat'.
    %
    %EXPORTFIG('picformat',PICFORMAT,...) or EXPORTFIG(PICFORMAT,...) sets
    %the fileformat for the external picture. Possibilities pdf {'pdf',
    %'-dpdf'}, png {'png', '-dpng'} or eps {'eps', '-depsc', '-depsc2'}. 
    %Default: 'png' 
    %
    %EXPORTFIG('figurehandle',FIGUREHANDLE,...) or EXPORTFIG(FIGUREHANDLE,...) 
    %explicitly specifies the handle of the figure that is to be stored. 
    %Default: gcf
    %
    %EXPORTFIG('-painters',...) or EXPORTFIG('-opengl',...) sets the renderer
    %for creation of the external figure. Use '-painters' for real vector
    %graphics. Default '-painters'
    %
    %EXPORTFIG('transparent',true/false) transparent background. Default: true
    %
    %EXPORTFIG('crop',[top,right,bottom,left],...) crops the image with the
    %given amount of pixels. If 'NaN' is given, the crop amount is
    %maximized. Default: [NaN,NaN,NaN,NaN] (tight cropping)
    %
    %EXPORTFIG('nocrop',...) equivalent to EXPORTFIG('crop',[0,0,0,0],...)
    %
    %EXPORTFIG('magnify',FACTOR,...) implements magnification of bitmap
    %pictures. Default: 1
    
    %read the input
    opts = inputParser(varargin{:});
    
    %prepare output structure
    out = struct('width',0,'height',0,'posCroppedRange',[0,0,0,0],'boundingBoxRelative',[0,0,1,1]);
    
    %store the current figure state
    figHandler = figureHandler(opts.figurehandle);
    figHandler.storeProperties;
    fig = figHandler.fig;
    
    %modify the figure for correct printing
    bgColor = modifyFigure(figHandler,opts);
    
    %--------------------------------------------------
    %create arrays with bitmap data
    %--------------------------------------------------
    [imData,alpha] = print2array(fig,opts,bgColor);
    out.width = size(imData,2);
    out.height = size(imData,1);
    
    %crop the image if requested
    [imData,alpha,rngCrop,bbRel] = cropBorders(imData,alpha,opts,bgColor);
    
    %position of cropped range as defined by matlab e.g. figure.Position
    out.posCroppedRange = [rngCrop(3), out.height-rngCrop(2)+1, rngCrop(4)-rngCrop(3)+1, rngCrop(2)-rngCrop(1)+1];
    out.boundingBoxRelative = bbRel;
    
    %--------------------------------------------------
    %create the files
    %--------------------------------------------------
    if any(strcmpi(opts.picformat,'png'))
        err = createPNG(imData,alpha,opts);
        restoreOnError(err);
    end
    if any(strcmpi(opts.picformat,'pdf'))
        err = createPDF(fig,opts,out.posCroppedRange/opts.magnify);
        restoreOnError(err);
    end
    
    %restore the figure state
    try figHandler.restoreGraphicsState(); catch, end
end
function restoreOnError(err,figHandler)
    %restores the figure if an error occured
    if ~isempty(err)
        try figHandler.restoreGraphicsState; catch, end
        rethrow(err);
    end
end
%%
%#######################################################################################################
% SUBFUNCTIONS - create picture files
%#######################################################################################################
function err = createPNG(imData,alpha,opts)
    err = [];
    
    %screen resolution
    res = opts.magnify * get(0, 'ScreenPixelsPerInch') / 25.4e-3;
    fileName = fullfile(opts.filefolder,[opts.filename,'.png']);
    pngOptions = {fileName, 'ResolutionUnit','meter', 'XResolution',res, 'YResolution',res};
    if opts.transparent
        pngOptions = [pngOptions 'Alpha',double(alpha)];
    end
    
    %actual image creation
    try 
        imwrite(imData,pngOptions{:});
    catch e
        err = e;
    end
end
function err = createPDF(fig,opts,posCroppedRange)
    err = [];
    
    %options for exportgraphics
    fileName = fullfile(opts.filefolder,[opts.filename,'.eps']);
    if opts.transparent
        bgColor = 'none';
    else
        bgColor = 'current';
    end
    if strcmpi(opts.renderer,'-painters')
        type = 'vector';
    else
        type = 'image';
        fig.PaperUnits = 'points';
        fig.PaperPosition = posCroppedRange;
    end
    
    %actual image creation
    try 
        exportgraphics(fig,fileName,'ContentType',type,'BackgroundColor',bgColor);
    catch e
        err = e;
    end
end
%%
%#######################################################################################################
% SUBFUNCTIONS - print2array
%#######################################################################################################
function [imData,alpha] = print2array(fig,opts,bgColor)
    %prints the figure to a image array. 
    %Based on print2array by export_fig
    
    %check for large arrays
    if prod(fig.Position(3:4))*opts.magnify>30e6
        warning('exportFig:largePicture','Large figure with more tan 30M pixel. Might cause out of memory errors')
    end
    
    % Set the resolution parameter
    resStr = ['-r' num2str(ceil(get(0, 'ScreenPixelsPerInch')*opts.magnify))];
    
    % use matlabs print function
    try 
        imData = print(fig, opts.renderer, resStr, '-RGBImage');
    catch e
        if contains(e.message,'memory')
            error('exportFig:largePicture','exportFig: out of memory')
        else
            rethrow e
        end
    end
    
    % Ensure that the output size is correct
    if isequal(opts.magnify, round(opts.magnify))
        px = round([fig.Position([4 3])*opts.magnify 3]);  % round() to avoid an indexing warning below
        if any(size(imData) > px) %~isequal(size(A), px)
            imData = imData(1:min(end,px(1)),1:min(end,px(2)),:);
        end
    end
    
    %handle transparancy
    imgSize = size(imData);
    alpha =  255*ones(imgSize(1:2), 'uint8');
    
    if opts.transparent
        isBgColor = imData(:,:,1)==bgColor(1) & imData(:,:,2)==bgColor(2) & imData(:,:,3)==bgColor(3);
        % Set the bgcolor pixels to be fully-transparent
        imData(repmat(isBgColor,[1,1,3])) = 255; %=white % TODO: more memory efficient without repmat
        alpha(isBgColor) = 0;
    end
    alpha = single(alpha)/255; %switch to [0,1] value
end
function [imData,alpha,rngCrop,bbRel] = cropBorders(imData,alpha,opts,bgColor)
    %crops the image by the given amount or determines the cropamount 
    %based on crop_borders by export_fig
    
    %compute the range of the picture after cropping (rngCrop) and the relative
    %bounding box for EPS (bbRel)
    [rngCrop,bbRel] = determineCropRange(imData,bgColor,opts.crop);
    
    %actual cropping of image data
    padding = 0; %for now no padding (ommitted this feature from export_fig)
    vA = [rngCrop(1)-padding rngCrop(2)+padding rngCrop(3)-padding rngCrop(4)+padding];
    imData = imData(vA(1):vA(2), vA(3):vA(4), :);
    alpha = alpha(vA(1):vA(2), vA(3):vA(4));
end
function [rngCrop,bbRel] = determineCropRange(imData,bgColor,cropAmounts)
    %crops the image by the given amount or determines the cropamount 
    %based on crop_borders by export_fig
    
    %dimensions of imData
    [h, w, c] = size(imData);
    
    % Crop margin from left
    if ~isfinite(cropAmounts(4))
        bail = false;
        for l = 1:w
            for a = 1:c
                if ~all(colVec(imData(:,l,a,:)) == bgColor(a))
                    bail = true;
                    break;
                end
            end
            if bail
                break;
            end
        end
    else
        l = 1 + abs(cropAmounts(4));
    end

    % Crop margin from right
    if ~isfinite(cropAmounts(2))
        bgColor = imData(ceil(end/2),w,:,1);
        bail = false;
        for r = w:-1:l
            for a = 1:c
                if ~all(colVec(imData(:,r,a,:)) == bgColor(a))
                    bail = true;
                    break;
                end
            end
            if bail
                break;
            end
        end
    else
        r = w - abs(cropAmounts(2));
    end

    % Crop margin from top
    if ~isfinite(cropAmounts(1))
        bgColor = imData(1,ceil(end/2),:,1);
        bail = false;
        for t = 1:h
            for a = 1:c
                if ~all(colVec(imData(t,:,a,:)) == bgColor(a))
                    bail = true;
                    break;
                end
            end
            if bail
                break;
            end
        end
    else
        t = 1 + abs(cropAmounts(1));
    end

    % Crop margin from bottom
    bgColor = imData(h,ceil(end/2),:,1);
    if ~isfinite(cropAmounts(3))
        bail = false;
        for b = h:-1:t
            for a = 1:c
                if ~all(colVec(imData(b,:,a,:)) == bgColor(a))
                    bail = true;
                    break;
                end
            end
            if bail
                break;
            end
        end
    else
        b = h - abs(cropAmounts(3));
    end
    
    %range of cropped picture
    rngCrop = [t,b,l,r];
    
    % For EPS cropping, determine the relative BoundingBox - bbRel
    bbRel = [l-1 h-b-1 r+1 h-t+1]./[w h w h];
end
function A = colVec(A)
    A = A(:);
end
%% 
%#######################################################################################################
% SUBFUNCTIONS - input arguments
%#######################################################################################################
function opts = inputParser(varargin)
    %default options
    opts = struct(...
        'filename',         'exportFig',...
        'filefolder',       '',...
        'renderer',         '-painters',...
        'picformat',        'pdf',...
        'figurehandle',     [],...
        'transparent',      true,...
        'crop',             [0,0,0,0],...
        'magnify',          1);
    
    ii = 1;
    setPicformat = false;
    while ii<=numel(varargin)
        if numel(varargin{ii})==1 && isgraphics(varargin{ii},'figure')
            varargin = horzcat(varargin,{'figurehandle',varargin{ii}});%#ok<AGROW>
        elseif ischar(varargin{ii}) && isvector(varargin{ii}) && ~isempty(varargin{ii})
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
                case 'transparent'
                    if isscalar(varargin{ii+1}) && islogical(varargin{ii+1})
                        opts.transparent = varargin{ii+1};
                        ii = ii+1;
                    elseif ischar(varargin{ii+1}) && isvector(varargin{ii+1})
                        opts.transparent = true;
                    else
                        errorBadOption('option "transparent" must be followed by true/false')
                    end
                case 'magnify'
                    if isnumeric(varargin{ii+1}) && isscalar(varargin{ii+1}) && varargin{ii+1}>0
                        opts.magnify = varargin{ii+1};
                        ii = ii+1;
                    else
                        errorBadOption('option "magnify" must be followed by a scalar > 0')
                    end
                case 'crop'
                    if isvector(varargin{ii+1}) && numel(varargin{ii+1})==4 && isnumeric(varargin{ii+1})
                        opts.crop = varargin{ii+1};
                        opts.crop(isinf(opts.crop)) = NaN;
                        zwNan = isnan(opts.crop);
                        if any(floor(opts.crop(~zwNan))~=ceil(opts.crop(~zwNan))) || any(opts.crop(~zwNan)<0)...
                        || ~isreal(opts.crop(~zwNan))
                            errorBadOption('crop must include NaN/Inf or real positive integers')
                        end
                        ii = ii+1;
                    elseif isempty(varargin{ii+1})
                        opts.corp = zeros(1,4);
                        ii = ii+1;
                    else
                        errorBadOption('option "crop" must be followed by vector [top,right,left,bottom]')
                    end
                case {'nocrop','-nocrop'}
                    varargin = horzcat(varargin,'crop',[0,0,0,0]);%#ok<AGROW>
                case {'picformat','pdf','eps','png','-dpng','-dpdf','-depsc','-depsc2'}
                    %warning for multiple settings of picformat
                    if setPicformat
                        warning('tikz3dPlot:DuplicatePicformat',['tikz3dPlot: Multiple picformats provided.',...
                            'Will use last one. Use either tikz3dPlot(...,picformat,format)',...
                            'or one of: pdf,eps,png,-dpng,-dpdf,-depsc,-depsc2'])
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
                otherwise
                    if strcmpi(opts.filename,'exportFig')
                        %set the filename to be found automatically 
                        opts.filename = varargin{ii};
                        varargin = horzcat(varargin,'filename',varargin{ii});%#ok<AGROW>
                    else
                        errorBadOption(['Unrecognized exportFig input option: ''' varargin{ii} '''']);
                    end
            end
        else
            %completely unknown option
            disp(varargin{ii})
            errorBadOption('Unrecognized exportFig input option (displayed above)');
        end
        ii = ii+1;
    end
    
    %try to get standard input for non set options
    if isempty(opts.figurehandle)
        opts.figurehandle = get(0, 'CurrentFigure');
        if isempty(opts.figurehandle)
            error('exportFig:BadOption','no open figure and no figurehandle provided')
        end
    end
    
    %check for number of axes on figure
    count = 0;
    for i=1:numel(opts.figurehandle.Children)
        if isa(opts.figurehandle.Children(i),'matlab.graphics.axis.Axes')
            count = count + 1;
        end
    end
    if count~=1
        warning('exportFig:numberOfAxes',['exportFig: figure should have exactly one axes object. ',...
            'unrespected results possible'])
    end
    
    %check if the folder exists
    if ~isempty(opts.filefolder) && ~(exist(opts.filefolder,'dir')==7)
        error('exportFig:fileIO','exportFig: provided folder "%s" for file not found',opts.filefolder)
    elseif isempty(opts.filefolder)
        opts.filefolder=pwd;%set folder to current path
    end
    
    %apply transparency if color of figure is 'none'
    if ischar(opts.figurehandle.Color) && strcmpi(opts.figurehandle.Color,'none')
        opts.transparent = true;
    end
end
function errorBadOption(message)
    error('exportFig:BadOption',['exportFig: ',message]);
end

%%
%#######################################################################################################
% SUBFUNCTIONS - plot handels
%#######################################################################################################
function bgColor = modifyFigure(figHandler,opts)
    %applies modifications to the figure in order for it to be printed
    %appropriately
    fig = figHandler.fig;
    
    %set axis limits to manual
    figHandler.setAxisLimits();
    
    %needed for ghostscript to work (crashes if not pixels - see
    %export_fig)
    fig.Units = 'pixels'; 
    
    %ensure some paper properties (fixes the print(...) output)
    fig.PaperOrientation = 'portrait';
    fig.PaperPositionMode = 'auto';
    
    %deal with transparency
    if opts.transparent
        fig.Color = 'w'; %set to white
    end
    bgColor = fig.Color; 
    
    %modify bgColor
    if all(bgColor<=1)
        bgColor = uint8(bgColor*255); %convert to rgb 0-255
    end
end