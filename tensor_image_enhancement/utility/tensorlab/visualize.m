function varargout = visualize(data, varargin)
%VISUALIZE Visualize a higher-order tensor.
%   1D, 2D and 3D plots of a tensor or a decomposed tensor can be constructed.
%   The other dimensions or variables are considered to be free and can be
%   adjusted using sliders. The data can be any tensor. The modulus is taken
%   for complex data.
%    
%   Options can be added as key-value pairs or as a struct, e.g.:
%           
%       visualize(T, 'FixedAxes', true, 'ShowContours', false);
%    
%       options = struct();
%       options.FixedAxes = true;
%       options.ShowContours = false;    
%       visualize(T, options);
%    
%   The following options can be specified (let N be the order of the data,
%   options with a * have different defaults depending on the type of data):
%   - plotDimensions:       logical vector of length N indicating free
%                           dimensions (Default: [1 1 0 ... 0])
%   - plotSliders:          vector of length N with indices for each slider
%                           (Default: [1 1 ... 1])
%   - dimensionLabels*:     Cell of strings with labels for each dimension.
%                           If this is a single string, it will be replicated
%                           for each dimension. %d will be replaced by the
%                           dimension number. (Default: 'i_%d')
%   - DimensionTransform    Function to transform the axes of the plot.
%                           Either, a single function handle or array a cell
%                           with N function handles or arrays. Some
%                           higher-order plots may not work with array type
%                           dimension transform (e.g., voxel3). (Default:
%                           @(i) i).
%   - fixedAxes:            If true, do not reset the axes each time the
%                           figure is updated. (Default: true)
%   - fixedView:            If true, do not reset the view each time the
%                           figure is updated. (Default: true)
%   - fixedColors:          If true, do not reset the colors each time the
%                           figure is updated. (Default: true)    
%   - surfOptions:          Extra options for the surf command (in case of a
%                           2D plot). (Default: {})
%   - showContours:         Show contours (in case of a 2D plot). (Default:
%                           false)
%   - fixContoursAtBottom:  Plot the contours at the bottom of the figure.
%                           (Default: true)
%   - contourOptions:       Extra options for the contour command (in case of a
%                           2D plot). (Default: {})
%   - outputLimits:         Limits for the y-axis (1D) or z-axis (3D). If
%                           empty, the plot/surf command defaults are used.
%                           (Default: [])
%   - holdPreviousAxes:     If true, do not reset the axes when calling
%                           visualize. This can be useful the sliders are set
%                           externally. (Default: false)
%   - hoplot:               Plot command to be used for 3D visualization.
%                           (Default: @slice3).
%   - enableDimension*:     Vector of logicals of size N. If true, the
%                           dimension can be selected  as an axis for
%                           plotting. If one element, the value is
%                           replicated. (Default: [true])
%   - tensorType:           Indicates how the data should be interpreted.
%                           Possible values: ['auto', 'full', 'cpd', 'tt',
%                           'btd', 'lmlra', 'sra', 'custom']. If the type is
%                           'custom', the functions getndims, getsize and
%                           getdata have to be provided by the user. (Default:
%                           'auto')
%   - original:             Original data to be plotted as dots on top of the
%                           given data.
%   - originalFormat:       Plotting format (string) for original data. This
%                           options is given to plot(x,y,originalFormat).
%                           (Default: '.');
%   - originalOptions       Plotting options for original data as a cell of
%                           key-value pairs. (Default: {})
%   - support:              Function returning a true/false matrix indicating
%                           whether the supplied indices are valid points.
%                           The function should have signature:
%
%                              sup = support(dim, ind)
%
%                           where dim is a logical array indicating if a
%                           dimension is plotted, and ind the current value
%                           of the sliders.
%   - trace:                For the 1D plots, do not erace the previous plots
%                           when changing the sliders if trace is true.
%
%   Handling different tensor types. All types known by getstructure can be
%   handled automatically. The type will be determined if tensortype =
%   'auto'. Other available types are:
%   - 'sra':      data is a cell containing the factor matrices for a
%                 successive rank-1 approximation. The defaults for the
%                 dimensionLabels are changed to {'i_%d', ..., 'i_%d', 'r'}
%                 and for enableDimension to [true, ..., true, false].
%   - 'custom':   a special or non-implemented decomposition. See further.
%
%   Handling custom types:
%   Three functions have to be provided, using the 'getsize', 'getndims' and
%   'getdata' options:
%   - function N = getndims(data):
%     % returns the number of dimensions N of the data, i.e. the number of
%     % sliders that are required. For example for the SRA, this is the
%     % dimension of the tensor + 1.
%   - function sz = getsize(data, dim):
%     % returns the size of the data in the requested dimension.
%   - function res = getdata(data, dim, ind)   
%     % returns the 1D, 2D or 3D data array where dim is a logical vector
%     % with ones for the free dimensions, and ind a integer vector with the
%     % current values of each slider. The number of dimensions for the
%     % output data can be determined using dim.
    
% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2016/03/14   NV      Updated help tooltips
% - 2016/01/08   NV      Added structured tensors as data model
% - 2015/12/15   NV      Added support for dimension transform
% - 2015/07/15   NV      Added support, original data and trace
% - 2014/02/07   NV      Initial version
    
    p = inputParser;
    p.addRequired('data');
    p.addOptional('plotDimensions', true);
    p.addOptional('plotSliders', true);
    p.addOptional('dimensionLabels', 'i_%d');
    p.addOptional('dimensionTransform', @(i) i);
    p.addOptional('fixedAxes', true);
    p.addOptional('fixedView', true);
    p.addOptional('fixedColors', true);
    p.addOptional('surfOptions', {});
    p.addOptional('showContours', false);
    p.addOptional('fixContoursAtBottom', true);
    p.addOptional('contourOptions', {});
    p.addOptional('outputLimits', []);
    p.addOptional('datatype', 'auto');
    p.addOptional('holdPreviousAxes', false);
    p.addOptional('hoplot', @slice3);
    p.addOptional('getsize', '');
    p.addOptional('getndims', '');
    p.addOptional('getdata', '');
    p.addOptional('enableDimension', 1, @islogical);
    p.addOptional('original', []);
    p.addOptional('originalType', 'auto');
    p.addOptional('originalFormat', '.');
    p.addOptional('originalOptions', {});
    p.addOptional('support', []);
    p.addOptional('trace', false);
    p.parse(data, varargin{:});
    
    isfunc = @(f) isa(f, 'function_handle');
    
    %% data type
    data = p.Results.data;
    datatype = p.Results.datatype;
    if ischar(datatype) && strcmpi(datatype, 'auto')
        try 
            datatype = getstructure(data);
        catch
            error('visualize:unknownType', ['The data type of the tensor could ' ...
                                'not be determined.']);
        end 
    end 
    
    switch datatype
      case {'full', 'cpd', 'btd', 'lmlra', 'tt', 'incomplete', 'sparse','loewner','hankel'}
        getdatasize = @(dim) getsize_generic(data, dim);
        getndims = @(~) getndims_generic(data);
        getdata = @(dim,ind)getdata_generic(datatype, data, dim, ind);
      case 'sra'
        getdatasize = @(dim) getsize_sra(data, dim);
        getndims = @(~) getndims_sra(dim);
        getdata = @(dim,ind) getdata_sra(data,dim,ind);
        if any(strcmp('enableDimension', p.UsingDefaults))
            enableDimension = [true(1, getndims(data)-1), false];
        end
        if any(strcmp('dimensionLabels', p.UsingDefaults))
            dimensionLabels = [repmat({'i_%d'}, 1, getndims(data)-1), 'r'];
        end
      case 'custom'
        getdatasize = @(dim) p.Results.getsize(data,dim);
        getndims = @(~) p.Results.getndims(data);
        getdata = @(dim,ind) p.Results.getdata(data, dim, ind);
        if ~isfunc(getdatasize) || ~isfunc(getndims) || ~isfunc(getdata)
            error('visualizetensor:customtype:notallfunctions', ...
                  'Not all given functions are functions.')
        end
      otherwise
        error('visualize:unknownType', ['The data type of the tensor could ' ...
                            'not be determined.']);
    end
        
    %% Original data
    original = p.Results.original;
    originalType = p.Results.originalType;
    if ischar(originalType) && strcmpi(originalType, 'auto')
        if isempty(original)
            originalType = 'none';
        else 
            try 
                originalType = getstructure(original);
            catch
                error('visualize:unknownType', ['The original data type of the tensor could ' ...
                                    'not be determined.']);
            end 
        end
    end 
    switch originalType
      case {'full', 'cpd', 'btd', 'lmlra', 'tt', 'incomplete', 'sparse'}
        getdataoriginal = @(dim,ind)getdata_generic(originalType, original, dim, ind);
        originalPresent = true;
      case 'none'
        originalPresent = false;
      otherwise 
        error('visualize:unknownType', ['The original data type of the tensor ' ...
                            'could not be determined.']);
    end
    originalFormat = p.Results.originalFormat;
    originalOptions = p.Results.originalOptions;
    
    %% Other options
    supportPresent = false;
    if isfunc(p.Results.support)
        support = p.Results.support;
        supportPresent = true;
    end
    
    plotDimensions = p.Results.plotDimensions;
    if islogical(plotDimensions)
        plotDimensions = [1 1 zeros(1,getndims(data)-2)];
    end
    plotSliders = p.Results.plotSliders;
    if islogical(plotSliders)
        plotSliders = ones(1,getndims(data));
    end
    if ~exist('dimensionLabels', 'var')
        dimensionLabels = p.Results.dimensionLabels;
    end
    dimTrans = p.Results.dimensionTransform;
    if ~iscell(dimTrans)
        dimTrans = repmat({dimTrans}, 1, getndims(data));
    end 
    fixedAxes = p.Results.fixedAxes;
    fixedView = p.Results.fixedView;
    fixedColors = p.Results.fixedColors;
    surfOptions = p.Results.surfOptions;
    if isstruct(surfOptions)
        surfOptions = [fieldnames(surfOptions), struct2cell(surfOptions)]';
    end
    showContours = p.Results.showContours;
    fixContoursAtBottom = p.Results.fixContoursAtBottom;
    contourOptions = p.Results.contourOptions;
    if isstruct(contourOptions)
        contourOptions = [fieldnames(contourOptions), struct2cell(contourOptions)]';
    end
    outputLimits = p.Results.outputLimits;
    holdPreviousAxes = p.Results.holdPreviousAxes;
    hoplot = p.Results.hoplot;

    if ~exist('enableDimension', 'var')
        enableDimension = p.Results.enableDimension;
        if length(enableDimension) == 1
            enableDimension = repmat(enableDimension, 1, getndims(data));
        elseif length(enableDimension) ~= getndims(data)
            error('visualize:enableDimension', ['Invalid number of entries: ' ...
                                'length(enableDimension) should be 1 or getndims']);
        end
    end
    sliderHeight = 30;
    
    N = getndims(data);
    if ~iscell(dimensionLabels), 
        dimensionLabels = repmat({dimensionLabels}, 1, N);
    end
    
    fh = gcf;
    clf;
    set(fh, 'menubar', 'figure');
    set(fh, 'toolbar', 'figure');
    defaultBackground = get(0,'defaultUicontrolBackgroundColor');
    set(fh,'Color',defaultBackground)
    
    ah = axes('Parent',fh,'Position',[0.15 0.15 0.7 0.7]);
    
    lbl = cell(N, 1);
    sh = cell(N, 1);
    eth = cell(N, 1);
    chb = cell(N, 1);
    slidercontainer = cell(N, 1);
    
    slidercontainer = uipanel(fh, ...
                              'units', 'pixels', ...
                              'BorderType', 'none', ...
                              'Tag', 'Sliders');
    defaultBackground = get(0,'defaultUicontrolBackgroundColor');
    set(slidercontainer,'BackgroundColor',defaultBackground)
    set(fh, 'ResizeFcn', @repositionsliders);
    
    function repositionsliders(~,~)
        old_units = get(fh,'Units');
        set(fh,'Units','pixels');
        figpos = get(fh,'Position');
        upos = [0, max(figpos(4) - N*sliderHeight,0), figpos(3), N*sliderHeight];
        set(slidercontainer,'Position',upos);
        upos = [0.15, 0.15, 0.7, max((figpos(4)-N*sliderHeight)/figpos(4)*0.7,0)];
        set(ah, 'Position', upos);
        set(fh,'Units', old_units);
    end

    y = sliderHeight*(N-1);
    
    for d = 1:N
        if getdatasize(d) == 1
            sliderStep = [0 0];
        else 
            sliderStep = [1 1] / (getdatasize(d)-1);
        end
        
        cbDimension = @(obj, evt) changePlotDimension(d, obj, evt);
        cbSlider = @(obj, evt) changePlotSlider(obj, evt, d);
                
        chb{d} = uicontrol(slidercontainer,'Style','checkbox',...
                           'String','',...
                           'Value',plotDimensions(d),'Position',[10 y 20 25],...
                           'Callback', cbDimension, ...
                           'ToolTipString', sprintf('test\nbla'));
        lbl{d} = uicontrol(slidercontainer,'Style','text',...
                           'String',sprintf(dimensionLabels{d}, d),...
                           'Position',[30 y 60 20],...
                           'Callback', cbDimension, ...
                           'HorizontalAlignment', 'left');
        sh{d} = uicontrol(slidercontainer,'Style','slider',... 
                          'Max',getdatasize(d), 'Min',1,'Value',plotSliders(d),...
                          'SliderStep',sliderStep,...
                          'Position',[110 y 150 25],...
                          'Callback',cbSlider,...
                          'KeyPressFcn', cbSlider,...
                          'ToolTipString', sprintf('Moving slider changes data'));      
        eth{d} = uicontrol(slidercontainer,'Style','edit',...
                           'String',num2str(dimTrans{d}(plotSliders(d))),...
                           'Position',[270, y, 50, 25],...
                           'KeyPressFcn',cbSlider, ...
                           'Callback',cbSlider);

        if ~enableDimension(d) 
            set(chb{d}, 'Enable', 'off')
        end
        
        y = y - sliderHeight;
    end

    previous = '';
    repositionsliders();
    changePlotDimension(1,chb{1})
    updateGUI();
    
    if nargout == 1
        varargout{1} = handle;
    end
    
    function changePlotSlider(hObject, eventdata, dim)
        try 
            % HG 2 version
            needUpdate = false;
            if strcmpi(get(hObject,'Style'), 'edit')
                % convert value to nearest index
                idx = dimTrans{dim}(1:getdatasize(dim));
                v = str2double(get(hObject, 'String'));
                [~,plotSliders(dim)] = min(abs(idx-v));
                
                if plotSliders(dim) < 1, plotSliders(dim) = 1; end
                if plotSliders(dim) > getdatasize(dim), 
                    plotSliders(dim) = getdatasize(dim);
                end
                needUpdate = true;
            else 
                temp = round(get(hObject, 'Value'));
                if plotSliders(dim) ~= temp
                    plotSliders(dim) = temp;
                    needUpdate = true;
                end
            end
        catch 
            % HG 1 version
            needUpdate = false;
            if isempty(eventdata) 
                plotSliders(dim) = round(get(hObject, 'Value'));
                needUpdate = true;
            else
                if strcmpi(eventdata.Key, 'rightarrow') || ...
                        strcmpi(eventdata.Key, 'leftarrow') 
                    temp = round(get(hObject, 'Value'));
                    if plotSliders(dim) ~= temp
                        plotSliders(dim) = temp;
                        needUpdate = true;
                    end
                end
                if strcmpi(eventdata.Key, 'return')
                    % convert value to nearest index
                    idx = dimTrans{dim}(1:getdatasize(dim));
                    v = str2double(get(hObject, 'String'));
                    [~,plotSliders(dim)] = min(abs(idx-v));

                    if plotSliders(dim) < 1, plotSliders(dim) = 1; end
                    if plotSliders(dim) > getdatasize(dim), 
                        plotSliders(dim) = getdatasize(dim);
                    end
                    needUpdate = true;
                end
            end
        end
        
        if needUpdate,
            set(eth{dim}, 'String', num2str(dimTrans{dim}(plotSliders(dim))));
            set(sh{dim},  'Value', plotSliders(dim));
            updateGUI();
        end
    end
    
    function changePlotDimension(dim, hObject, ~)
        plotDimensions(dim) = get(hObject, 'Value');
        if sum(plotDimensions) <= 2 || (sum(plotDimensions) == 3 && ...
                                       strcmpi(func2str(hoplot),'voxel3'))
            for i = 1:length(chb)
                if plotDimensions(i)
                    set(eth{i}, 'Enable', 'off');
                    set(sh{i}, 'Enable', 'off');
                    set(eth{i}, 'String', '');
                else
                    set(eth{i}, 'Enable', 'on');
                    set(sh{i}, 'Enable', 'on');
                    set(eth{i}, 'String', num2str(dimTrans{i}(plotSliders(i))));
                end
                set(chb{i}, 'Enable', 'on');
            end
            if sum(plotDimensions) == 1
                d = find(plotDimensions);
                set(chb{d}, 'Enable', 'off');
            end
        else 
            for i = 1:length(chb)
                set(eth{i}, 'Enable', 'on');
                set(sh{i},  'Enable', 'on');
                set(chb{i}, 'Enable', 'on');
                set(eth{i}, 'String', num2str(dimTrans{i}(plotSliders(i))));
            end
        end
        for i = 1:length(chb)
            if plotDimensions(i)
                if any(strcmpi(func2str(hoplot), {'surf3'}))
                    txt1 = sprintf(['Move slider to change plotted surface.\n' ...
                                    'Some surfaces may be out of view.\nMove ' ...
                                    'sliders to see the surfaces']);
                    txt2 = sprintf(['Change value to move plotted surface.\n' ...
                                    'Some surfaces may be out of view.\nChange ' ...
                                    'values to see the surfaces']);
                else 
                    txt1 = 'Move slider to change plotted surface.';
                    txt2 = 'Change value to move plotted surface.';

                end
                set(sh{i},  'ToolTipString', txt1);
                set(eth{i}, 'ToolTipString', txt2);
            else                 
                txt = sprintf(['The data tensor is sliced at %s = %d.\n',...
                               'Move slider to change %s.'], ...
                              sprintf(dimensionLabels{i},i), ...
                              dimTrans{i}(plotSliders(i)), ...
                              sprintf(dimensionLabels{i},i));
                set(sh{i}, 'ToolTipString', txt);
                txt = sprintf(['The data tensor is sliced at %s = %d.\n',...
                               'Update value to change %s.'], ...
                              sprintf(dimensionLabels{i},i), ...
                              dimTrans{i}(plotSliders(i)), ...
                              sprintf(dimensionLabels{i},i));
                set(eth{i}, 'ToolTipString', txt);
            end
            if get(chb{i}, 'Value')
                txt = sprintf(['This dimension is an axis in the plot.\n', ...
                               'Uncheck to use this dimension to slice the ' ...
                               'data.']);
                set(chb{i}, 'ToolTipString', txt);
                set(lbl{i}, 'ToolTipString', txt);
            else 
                txt = sprintf(['This dimension is used to slice the data.\n', ...
                               'Check to plot this dimension.']);
                set(chb{i}, 'ToolTipString', txt);
                set(lbl{i}, 'ToolTipString', txt);
            end
        end
            
        updateGUI();
    end
    
    function updateGUI()
        if sum(plotDimensions) == 0
            previous = 'none';
        elseif sum(plotDimensions) == 1
            if ~p.Results.trace;
                if strcmpi(previous, 'plot') || holdPreviousAxes
                    set(ah,'NextPlot','replacechildren')
                else 
                    set(ah,'NextPlot','replace')                
                end          
            end
            lim = ylim;

            res = getdata(plotDimensions, plotSliders);
            d = find(plotDimensions);
            x = dimTrans{d(1)}(1:length(res(:)));
            if supportPresent
                ind = cellfun(@(d,i) d(i), dimTrans, num2cell(plotSliders), ...
                              'UniformOutput', false);
                ind{d(1)} = x;
                res(~support(plotDimensions, ind)) = nan;
            end
            
            handle = plot(x(:), res(:)); 
            fulldims = find(plotDimensions);
            xlabel(sprintf(dimensionLabels{fulldims(1)}, fulldims(1)));
            if ~isempty(outputLimits), ylim(outputLimits); end

            miny = min(res(:));
            maxy = max(res(:));

            if p.Results.trace, hold on; 
            else hold off; end
            
            if originalPresent
                if p.Results.trace
                    set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex') -1);
                end
                holding = ishold;
                hold on
                res = getdataoriginal(plotDimensions, plotSliders);
                plot(x, res(:), originalFormat, originalOptions{:});
                if ~holding, hold off; end
                
                miny = min(miny, min(res(:)));
                maxy = max(maxy, max(res(:)));
            end
            
            if strcmpi(previous, 'plot') || holdPreviousAxes
                if fixedAxes, 
                    ylim([min(miny, lim(1)) max(maxy,lim(2))]);
                end
            else 
                xlim([min(x), max(x)]);
            end
            
            previous = 'plot';            
        elseif sum(plotDimensions) == 2
            if strcmpi(previous, 'surf') || holdPreviousAxes
                set(ah,'NextPlot','replacechildren')
                [az,el] = view;
                lim = zlim;
            else 
                set(ah,'NextPlot','replace')
            end
            
            if showContours
                sopt = [surfOptions(:).', contourOptions(:).'];
                surfcmd = @surfc;
            else 
                sopt = surfOptions(:).';
                surfcmd = @surf;
            end
            res = getdata(plotDimensions, plotSliders);
            d = find(plotDimensions);
            x = dimTrans{d(2)}(1:size(res,2));
            y = dimTrans{d(1)}(1:size(res,1));
            [X,Y] = meshgrid(x,y);
            if supportPresent
                ind = cellfun(@(d,i) d(i), dimTrans, num2cell(plotSliders), ...
                              'UniformOutput', false);
                ind{d(1)} = y;
                ind{d(2)} = x;
                res(~support(plotDimensions, ind)) = nan;
            end
            if ~isempty(surfOptions)
                handle = surfcmd(X,Y, res, sopt{:});
            else
                handle = surfcmd(X,Y, res, sopt{:});
            end
            if ~isempty(outputLimits), zlim(outputLimits); end
            if fixedColors, caxis(zlim); end
            
            miny = min(res(:));
            maxy = max(res(:));
            
            if originalPresent
                hold on
                res = getdataoriginal(plotDimensions, plotSliders);
                plot3(X(:),Y(:),res(:),originalFormat,originalOptions{:});
                hold off
                
                miny = min(miny, min(res(:)));
                maxy = max(maxy, max(res(:)));
            end
            
            xlim([min(x) max(x)]);
            ylim([min(y) max(y)]);

            fulldims = find(plotDimensions);
            ylabel(sprintf(dimensionLabels{fulldims(1)}, fulldims(1)));
            xlabel(sprintf(dimensionLabels{fulldims(2)}, fulldims(2)));
            if strcmpi(previous, 'surf') || holdPreviousAxes
                if fixedView, view(az,el); end
                if fixedAxes, 
                    zlim([min(miny, lim(1)) max(maxy,lim(2))]);
                end
                if fixedColors, caxis(lim); end
            end
            
            previous = 'surf';            
        elseif sum(plotDimensions) == 3
            if strcmpi(previous, 'voxel3') || holdPreviousAxes
                set(ah,'NextPlot','replacechildren')
            else 
                set(ah,'NextPlot','replace')                
            end            
            lim = axis;

            % dim transform 
            d = find(plotDimensions);
            ind = plotSliders(d);
            if any(strcmpi(func2str(hoplot), {'slice3', 'voxel3'}))
                opts = {'DimensionTransform', dimTrans(d)};
                for n = 1:3
                    ind(n) = dimTrans{d(n)}(plotSliders(d(n)));
                end
            else 
                opts = {};
            end
            if any(strcmpi(func2str(hoplot), {'surf3', 'slice3', 'voxel3'}))
                opts = [opts(:).', {'Values', ind}];
            end
            if any(strcmpi(func2str(hoplot), {'voxel3'}))
                fast = true;
                for n = 1:3
                    tmp = dimTrans{n}(1:getdatasize(d(n)));
                    tmp = diff(tmp);
                    if any(tmp~=tmp(1)), fast = false; end
                end
                opts = [opts(:).', {'fast', fast}];
            end
            if strcmpi(func2str(hoplot), 'surf3')
                out = hoplot(getdata(plotDimensions, plotSliders), opts{:}); 
                xlim(out.scales(1:2));
                ylim(out.scales(3:4));
                zlim(out.scales(5:6));
            else 
                hoplot(getdata(plotDimensions, plotSliders), opts{:}); 
            end
            fulldims = find(plotDimensions);
            xlabel(sprintf(dimensionLabels{fulldims(1)}, fulldims(1)));
            if ~isempty(outputLimits), ylim(outputLimits); end

            if strcmpi(previous, 'voxel3') || holdPreviousAxes
                if fixedAxes, axis(lim); end
            end
            
            previous = 'voxel3';     
            fulldims = find(plotDimensions);
            zlabel(sprintf(dimensionLabels{fulldims(1)}, fulldims(1)));
            ylabel(sprintf(dimensionLabels{fulldims(2)}, fulldims(2)));
            xlabel(sprintf(dimensionLabels{fulldims(3)}, fulldims(3)));           
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DATA Models
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function res = getdata_generic(~, data, dim, ind)
        ind = num2cell(ind);
        ind(dim==1) = repmat({':'}, 1, sum(dim));
        res = ful(data,ind{:});
        res = squeeze(res);
        if ~isreal(res), res = abs(res); end
    end
    
    function sz = getsize_generic(data, dim)
        sz = getsize(data);
        sz = sz(dim);
    end 
    
    function N = getndims_generic(data)
        N = getorder(data);
    end
    
    %% Successive Rank-1 Approximation
    function res = getdata_sra(data, dim, ind) 
        ind = num2cell(ind);
        r = ind{getndims(data)};
        val = ones(1,r);
        for d = find(~dim(1:end-1))
            val = val .* data{d}(ind{d}, 1:r);
        end
        if sum(dim) == 1
            res = data{dim}(:,1:r)*val';
        elseif sum(dim) == 2
            dim = find(dim);
            res = data{dim(1)}(:,1:r)*diag(val)*data{dim(2)}(:,1:r)';
        else 
            dim = find(dim);
            res = cpdgen({data{dim(1)}(:,1:r),data{dim(2)}(:,1:r),data{dim(3)}(:,1:r), val});
        end
        if ~isreal(res), res = abs(res); end
    end
    
    function sz = getsize_sra(data, dim)
        if dim <= getndims(data) - 1 
            sz = size(data{dim}, 1);
        else
            sz = size(data{1}, 2);
        end
    end
    
    function N = getndims_sra(data)
        N = length(data)+1;
    end
    
end
