function slice3(T,varargin)
%SLICE3 Visualize a third-order tensor with slices.
%   SLICE3(T) visualizes the third-order tensor T by drawing its mode-1,
%   -2, and -3 slices using sliders to define their respective indices.
%   Press 'h' to show/hide the figure's controls.
%   
%   SLICE3(T,options) or SLICE3(T,'key',value) can be used to set the following
%   (optional) options:
%
%   - External              Array of length 3 containing the indices of the
%                           slices to plot. 
%   - DimensionTransform    Function to transform the axes of the plot.
%                           Either, a single function handle or array a cell
%                           with N function handles or arrays. (Default:
%                           @(i) i).
%
%   All other options are passed on as options to the plot.

%   Authors: Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   Version history:
%   - 2016/01/15   NV    Added External and DimensionTransform

% Check options
p = inputParser;
p.addOptional('DimensionTransform', @(i) i);
p.addOptional('Values', nan(1,3));
p.KeepUnmatched = true;
p.parse(varargin{:});
dimTransf = p.Results.DimensionTransform;
if ~iscell(dimTransf)
    dimTransf = repmat({dimTransf}, 1, 3);
elseif length(dimTransf) == 1;
    dimTransf = repmat(dimTransf, 1, 3);
end
if length(dimTransf) ~= 3
    error('slice3:dimensionTransform', ['DimensionTransform should be a single ' ...
                        'function handle or arrays, or a cell of exactly ' ...
                        'three function handles or arrays.']);
end
sliders = p.Results.Values;
if all(~isnan(sliders))
    for n = 1:ndims(T)
        sliders(n) = find(dimTransf{n}(1:size(T,n)) == sliders(n));
    end
end

% Check the dimensions.
T = ful(T);
if ndims(T) <= 2
    imagesc(T,varargin{:}); return;
elseif ndims(T) ~= 3
    error('slice3:T','ndims(T) should be 3.');
end

% If the data is complex, convert it to the modulus.
% Remove any NaNs and Infs.
if any(~isreal(T(:))), T = abs(T); end
T((isinf(T) & T < 0)) = min(T(~isinf(T(:))));
T(isinf(T) & T > 0) = max(T(~isinf(T(:))));
mn = min(T(:));
mx = max(T(:));

% Set up the figure.
idx = size(T); idx(2) = 1;
T = double(permute(T,[2 3 1]));
cax = newplot;
if all(~isnan(sliders))
    init = sliders;
else
    init = [size(T,3), 1, size(T,2)];
end 
set (gcf,'Toolbar','figure');
set(datacursormode(gcf),'UpdateFcn',@datacursor);
zoom off; pan off; rotate3d off; datacursormode off;
obj{1} = uicontrol('Style','text','Position',[5 45 15 15],'String','i');
obj{2} = uicontrol('Style','slider','Position',[25 45 120 15], ...
                   'Min',1,'Max',size(T,3),'Value',init(1), ...
                   'SliderStep',1/(size(T,3)-1)*[1 1], ...
                   'Callback',{@redraw,1}, ...
                   'KeyPressFcn',@(obj,evt)toggle(evt.Key));
obj{3} = uicontrol('Style','text','Position',[5 25 15 15],'String','j');
obj{4} = uicontrol('Style','slider','Position',[25 25 120 15], ...
                   'Min',1,'Max',size(T,1),'Value',init(2), ...
                   'SliderStep',1/(size(T,1)-1)*[1 1], ...
                   'Callback',{@redraw,2}, ...
                   'KeyPressFcn',@(obj,evt)toggle(evt.Key));
obj{5} = uicontrol('Style','text','Position',[5 5 15 15],'String','k');
obj{6} = uicontrol('Style','slider','Position',[25 5 120 15], ...
                   'Min',1,'Max',size(T,2),'Value',init(3), ...
                   'SliderStep',1/(size(T,2)-1)*[1 1], ...
                   'Callback',{@redraw,3}, ...
                   'KeyPressFcn',@(obj,evt)toggle(evt.Key));
if ~exist('sliders', 'var')
    isVisible = true;
    set(gcf,'KeyPressFcn',@(obj,evt)toggle(evt.Key));
else 
    isVisible = false;
end 
if any(~isnan(sliders))
    redraw(obj{2},0,1);
    redraw(obj{4},0,2);
    redraw(obj{6},0,3);
    for i = 1:length(obj), set(obj{i},'Visible','off'); end
else
    redraw();
end

function toggle(key)
    if strcmpi(key,'h')
        if isVisible, val = 'off'; else val = 'on'; end
        for i = 1:length(obj), set(obj{i},'Visible',val); end
        isVisible = ~isVisible;
    end
end

function txt = datacursor(~,event_obj)
    pos = get(event_obj,'Position');
    ind(3) = find(dimTransf{1}(1:size(T,3)) == pos(3));
    ind(2) = find(dimTransf{2}(1:size(T,1)) == pos(2));
    ind(1) = find(dimTransf{3}(1:size(T,2)) == pos(1));
    txt = {['I: ' int2str(pos(3))], ...
           ['J: ' int2str(pos(2))], ...
           ['K: ' int2str(pos(1))], ...
           ['Value: ' num2str(T(ind(2),ind(1),ind(3)))]};
end

function redraw(hobj,~,ax)
    
    % Draw slices.
    if nargin >= 1, idx(ax) = round(get(hobj,'Value')); end
    x = dimTransf{2}(1:size(T,1));
    y = dimTransf{3}(1:size(T,2));
    z = dimTransf{1}(1:size(T,3));
    [x,y,z] = meshgrid(y,x,z);
    if mn == mx, T(end) = T(end)*(1+eps); end
    slice(x,y,z,T,dimTransf{3}(idx(3)),dimTransf{2}(idx(2)),dimTransf{1}(idx(1)));
    shading flat;
    
    % Set axis properties.
    xlabel('k');
    ylabel('j');
    zlabel('i');
    xlim(dimTransf{3}([1 size(T,2)]));
    ylim(dimTransf{2}([1 size(T,1)]));
    zlim(dimTransf{1}([1 size(T,3)]));
    set(cax,'YDir','reverse');
    set(cax,'ZDir','reverse');
    
    % Display grid.
    step = max(1,round(size(T)/6));
    set(cax,'XTickMode','manual','YTickMode','manual','ZTickMode','manual');
    set(cax,'XTick',dimTransf{3}([1:step(2):size(T,2)-1 size(T,2)]));
    set(cax,'YTick',dimTransf{2}([1:step(1):size(T,1)-1 size(T,1)]));
    set(cax,'ZTick',dimTransf{1}([1:step(3):size(T,3)-1 size(T,3)]));
    grid on;
end

end
