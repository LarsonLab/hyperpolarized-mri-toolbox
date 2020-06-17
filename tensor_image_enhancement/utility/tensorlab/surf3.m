function varargout = surf3(T,varargin)
%SURF3 Visualize a third-order tensor with surfaces.
%   SURF3(T) visualizes the third-order tensor T by drawing its mode-1, -2,
%   and -3 slices using sliders to define their respective indices. Press
%   'h' to show/hide the figure's controls.
% 
%   SURF3(T,options) or SURF3(T,'key',value) can be used to set the following
%   (optional) options:
%
%   - Values                Array of length 3 containing the indices of the
%                           slices to plot. 
%
%   All other options are passed as options to the plot.

%   Authors: Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   Version history:
%   - 2016/01/15   NV    Added Values

% Check options
p = inputParser;
p.addOptional('Values', nan(1,3));
p.addOptional('DimensionTransform', @(i) i);
p.KeepUnmatched = true;
p.parse(varargin{:});
sliders = p.Results.Values;
transf = p.Results.DimensionTransform;
if ~iscell(transf)
    transf = repmat({transf}, 1, 3);
elseif length(transf) == 1;
    transf = repmat(transf, 1, 3);
end
if length(transf) ~= 3
    error('surf3:dimensionTransform', ['DimensionTransform should be a single ' ...
                        'function handle or arrays, or a cell of exactly ' ...
                        'three function handles or arrays.']);
end
if ~any(strcmpi(p.UsingDefaults, 'DimensionTransform'))
    for n = 1:length(transf)
        if isa(transf{n}, 'function_handle')
            transf{n} = transf{n}(1:getsize(T,n));
        end
        if length(transf{n}) ~= getsize(T,n)
            error('surf3:dimensionTransform', ['length(DimensionTransform{n}) ' ...
                                'should be getsize(T,n) if DimensionTransform ' ...
                                'contains arrays.'])
        end
        if any(transf{n}(:).' ~= 1:getsize(T,n))
            error('surf3:dimensionTransform', ['Surf3 does not support dimension ' ...
                                'transform at this moment.']);
        end
    end
end


% Check the dimensions.
T = ful(T);
if ndims(T) <= 2
    imagesc(T,p.Unmatched); return;
elseif ndims(T) ~= 3
    error('surf3:T','ndims(T) should be 3.');
end

% If the data is complex, convert it to the modulus.
% Remove any NaNs and Infs.
if any(~isreal(T(:))), T = abs(T); end
T(isnan(T) | (isinf(T) & T < 0)) = min(T(~isinf(T(:))));
T(isinf(T) & T > 0) = max(T(~isinf(T(:))));

% Set up the figure.
idx = size(T); idx(2) = 1;
T = double(permute(T,[2 3 1]));
if all(~isnan(sliders))
    init = sliders;
else
    init = [size(T,3), 1, size(T,2)];
end
mn = min(T(:));
mx = max(T(:));
cax = newplot;
scales = zeros(1, 6);
set(gcf,'Toolbar','figure');
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
areControlsVisible = true;
if any(isnan(sliders))
    set(gcf,'KeyPressFcn',@(obj,evt)toggle(evt.Key));
else 
    toggle('h')
    areControlsVisible = false;
end 
if all(~isnan(sliders))
    redraw(obj{2},0,1);
    redraw(obj{4},0,2);
    redraw(obj{6},0,3);
    for i = 1:length(obj), set(obj{i},'Visible','off'); end
    
    if nargout > 0
        varargout{1} = struct;
        limits = [zlim ylim xlim];
        tmp = -scale(T(:,:,idx(1)), 3)+idx(1);
        limits(1) = min(limits(1), min(tmp(:)));
        limits(2) = max(limits(2), max(tmp(:)));
        tmp = scale(T(idx(2),:,:), 1);
        limits(3) = min(limits(3), min(tmp(:)));
        limits(4) = max(limits(4), max(tmp(:)));
        tmp = scale(T(:,idx(3),:), 2);
        limits(5) = min(limits(5), min(tmp(:)));
        limits(6) = max(limits(6), max(tmp(:)));
        varargout{1}.scales = limits;
    end
else
    redraw();
end

if nargout == 0
    varargout = {};
else 
    varargout = {struct('scales', [xlim ylim zlim])};
end

function toggle(key)
    if strcmpi(key,'h')
        if areControlsVisible, val = 'off'; else val = 'on'; end
        for i = 1:length(obj), set(obj{i},'Visible',val); end
        areControlsVisible = ~areControlsVisible;
    end
end

function slice = scale(slice,n)
    if all(slice(:) > 0)
        slice = slice-min(slice(:));
    elseif all(slice(:) < 0)
        slice = slice+max(slice(:));
    end
    if mx-mn > 0, slice = (slice-mn)/(mx-mn)*size(T,n)/5; end
end

function redraw(hobj,~,ax)

    % Draw slices.
    if nargin >= 1, idx(ax) = round(get(hobj,'Value')); end
    sJK = T(:,:,idx(1));
    surf(cax, ...
         transf{3}(1:size(T,2)),...
         transf{2}(1:size(T,1)),...
         transf{1}(-scale(sJK,3)+idx(1)),...
         sJK); 
    hold on;
    sIK = reshape(T(idx(2),:,:),size(T,2),[]).';
    sIK = surf(cax,...
               transf{3}(1:size(T,2)),...
               transf{2}(idx(2)-1+(1:size(T,3))),...
               transf{1}(-scale(sIK,1)+1),...
               sIK);
    rotate(sIK,[1 0 0],90,[1 idx(2) 1]);
    sIJ = reshape(T(:,idx(3),:),size(T,1),[]);
    sIJ = surf(cax,...
               transf{3}(idx(3)-1+(1:size(T,3))),...
               transf{2}(1:size(T,1)),...
               transf{1}(scale(sIJ,2)+1),...
               sIJ);
    rotate(sIJ,[0 1 0],-90,[idx(3) 1 1]); hold off;
    shading flat;
    
    % Set axis properties.
    caxis([mn mx]);
    xlabel('k');
    ylabel('j');
    zlabel('i');
    xlim(transf{3}([1 size(T,2)]));
    ylim(transf{2}([1 size(T,1)]));
    zlim(transf{1}([1 size(T,3)]));
    set(cax,'YDir','reverse');
    set(cax,'ZDir','reverse');

    % Display grid.
    step = max(1,round(size(T)/6));
    set(cax,'XTickMode','manual','YTickMode','manual','ZTickMode','manual');
    set(cax,'XTick',transf{3}([1:step(2):size(T,2)-1 size(T,2)]));
    set(cax,'YTick',transf{2}([1:step(1):size(T,1)-1 size(T,1)]));
    set(cax,'ZTick',transf{1}([1:step(3):size(T,3)-1 size(T,3)]));
    grid on;

end

end
