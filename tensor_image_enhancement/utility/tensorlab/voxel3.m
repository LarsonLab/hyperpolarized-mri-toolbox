function voxel3(T,varargin)
%VOXEL3 Visualize a third-order tensor with voxels.
%   voxel3(T) visualizes the third-order tensor T by plotting its elements
%   as voxels whose color and opacity are proportional to their value. The
%   figure contains two sliders for setting the parameters thresh and
%   degree (press 'h' to hide/show them). Let alpha = (T(i,j,k)-min(T(:))/
%   (max(T(:))-min(T(:))), then the opaqueness of each voxel, where 0 is
%   transparent and 1 is opaque, is computed as
%      
%      alpha                            if alpha >= thresh
%      thresh^(1-degree)*alpha^degree   if alpha <  thresh
%
%   VOXEL3(T, options) or VOXEL3(T, 'key', value) can be used to set the
%   following options:
%
%   - Fast = [{true}|false]      - If false, a slower but more accurate
%                                  algorithm to draw the voxel plot is used
%   - Degree = 1                 - Degree parameter as defined above
%   - Tresh = 1                  - Threshold parameter as defined above
%   - DimensionTransform =       - Function to transform the axes of the
%     [@(i) i]                     plot. Either, a single function handle or
%                                  a cell with three function handles.
    
%   Authors: Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
% 
% Version History:
% - 2016/03/06   NV      Added dimension transform & option parsing

    
% Check options.
p = inputParser;
p.addOptional('Tresh', 0.5);
p.addOptional('Degree', 1);
p.addOptional('Fast', true);
p.addOptional('DimensionTransform', @(i) i);
p.addOptional('Values', nan); % ignored here, but needed for visualize
p.KeepUnmatched = true;
p.parse(varargin{:});
degree = p.Results.Degree;
tresh  = p.Results.Tresh;
fast   = p.Results.Fast;
transf = p.Results.DimensionTransform;
if ~iscell(transf)
    transf = repmat({transf}, 1, 3);
elseif length(transf) == 1;
    transf = repmat(transf, 1, 3);
end
if length(transf) ~= 3
    error('voxel3:dimensionTransform', ['DimensionTransform should be a single ' ...
                        'function handle or arrays, or a cell of exactly ' ...
                        'three function handles or arrays.']);
end

% If the data is complex, convert it to the modulus.
% Remove any NaNs and Infs.
T = ful(T);
if any(~isreal(T(:))), T = abs(T); end
T(isnan(T) | (isinf(T) & T < 0)) = min(T(~isinf(T(:))));
T(isinf(T) & T > 0) = max(T(~isinf(T(:))));
st = [size(T) 1 1];
mn = min(T(:));
mx = max(T(:));

% Compute dimension transforms for arrays (such that +- 0.5 works)
for i = 1:length(transf)
    if ~isnumeric(transf{i}), continue; end
    [dt,j] = sort(transf{i});
    dt = dt(:).';
    % check if it can be plotted using fast
    tmp = diff(dt);
    if fast && any(tmp ~= tmp(1)) 
        error('voxel3:incompatibleOptions', ...
              ['The option fast=true and non-equidistant dimension transforms ' ...
               'cannot be combined, yet. Set the option fast = false.']);
    end
    % sort data for increased efficiency
    sub = repmat({':'}, 1, 3);
    sub{i} = j;
    T(sub{:}) = T;
    % compute mid points
    midp = 0.5*(dt(1:end-1)+dt(2:end));
    midp = [2*dt(1)-midp(1), midp];
    % insert in dt
    dt = reshape([midp(:).'; dt(:).'], 1, []);
    dt = [dt 2*dt(end)-midp(end)];
    % convert to function
    transf{i} = @(j) dt(round(j*2));
end

% Compute the vertices and the minimal number of colors and faces.
if ~fast
    sv = st(1:3)+1;
    verts = zeros(prod(sv),3,'single');
    verts(:,3) = transf{1}(repmat((.5:1:st(1)+.5).',size(verts,1)/sv(1),1));
    verts(:,2) = transf{2}(repmat(kron((.5:1:st(2)+.5).',ones(sv(1),1)),sv(3),1));
    verts(:,1) = transf{3}(kron((.5:1:st(3)+.5).',ones(size(verts,1)/sv(3),1)));
    idx = {(1:prod(st(1:3)+[1 0 0]))',(1:prod(st(1:3)+[0 1 0]))', ...
        (1:prod(st(1:3)+[0 0 1]))'};
    off = [0 cumsum(cellfun(@numel,idx))];
    color = zeros(sum(cellfun(@numel,idx)),1,'single');
    faces = zeros(sum(cellfun(@numel,idx)),4,'single');
    Tint = cat(1,T(1,:,:),max(T(1:end-1,:,:),T(2:end,:,:)),T(end,:,:));
    [i,j,k] = ind2sub(size(Tint),idx{1});
    color(off(1)+idx{1}) = single(Tint(:));
    faces(off(1)+idx{1},:) = [sub2ind(sv,i,j  ,k)   sub2ind(sv,i,j+1,k) ...
                              sub2ind(sv,i,j+1,k+1) sub2ind(sv,i,j  ,k+1)];
    Tint = cat(2,T(:,1,:),max(T(:,1:end-1,:),T(:,2:end,:)),T(:,end,:));
    [i,j,k] = ind2sub(size(Tint),idx{2});
    color(off(2)+idx{2}) = single(Tint(:));
    faces(off(2)+idx{2},:) = [sub2ind(sv,i  ,j,k)   sub2ind(sv,i+1,j,k) ...
                              sub2ind(sv,i+1,j,k+1) sub2ind(sv,i  ,j,k+1)];
    Tint = cat(3,T(:,:,1),max(T(:,:,1:end-1),T(:,:,2:end)),T(:,:,end));
    [i,j,k] = ind2sub(size(Tint),idx{3});
    color(off(3)+idx{3}) = single(Tint(:));
    faces(off(3)+idx{3},:) = [sub2ind(sv,i  ,j  ,k) sub2ind(sv,i+1,j,k) ...
                              sub2ind(sv,i+1,j+1,k) sub2ind(sv,i  ,j+1,k)];
	alpha = (color-mn)/(mx-mn);
    cutoff = 1e-4;
end

% Set up the figure.
cax = newplot;
options = {'FaceColor','texturemap','CDataMapping','scaled', ...
           'FaceAlpha','texturemap','AlphaDataMapping','scaled', ...
           'EdgeColor','none','Parent',cax};
k = degree; t = tresh;
set(gcf,'Toolbar','figure');
zoom off; pan off; rotate3d off; datacursormode off;
lbl1 = uicontrol('Style','text','Position',[5 25 65 15], ...
                 'HorizontalAlignment','left');
sld1 = uicontrol('Style','slider','Position',[75 25 120 15], ...
                 'Min',0,'Max',1,'Value',0.5, ...
                 'SliderStep',[0.1 0.25], ...
                 'Callback',{@redraw,'t'}, ...
                 'KeyPressFcn',@(obj,evt)toggle(evt.Key));
lbl2 = uicontrol('Style','text','Position',[5 5 65 15], ...
                 'HorizontalAlignment','left');
sld2 = uicontrol('Style','slider','Position',[75 5 120 15], ...
                 'Min',1,'Max',5,'Value',1, ...
                 'SliderStep',[0.25 0.5], ...
                 'Callback',{@redraw,'k'}, ...
                 'KeyPressFcn',@(obj,evt)toggle(evt.Key));
isVisible = true;
set(gcf,'KeyPressFcn',@(obj,evt)toggle(evt.Key));

redraw();

function toggle(key)
    % Toggle show or hide of uicontrols.
    if strcmpi(key,'h')
        onoff = 'on'; if isVisible, onoff = 'off'; end
        cellfun(@(c)set(c,'Visible',onoff),{lbl1,sld1,lbl2,sld2});
        isVisible = ~isVisible;
    end
end

function adata = adata(cdata)
    % Compute opaqueness from color data.
    adata = (squeeze(cdata)-mn)/(mx-mn);
    adata(adata<t) = t^(1-k)*adata(adata<t).^k;
end

function redraw(hobj,~,param)
    
    % Update the parameter.
    if nargin > 1
        switch param
            case 't', t = get(hobj,'Value');
            case 'k', k = get(hobj,'Value');
        end
    end
    if exist('lbl1','var')
        set(lbl1,'String',sprintf('thresh = %g',t));
        set(lbl2,'String',sprintf('degree = %g',k));
    end

    % Save axis information
    if ~any(strcmpi(p.UsingDefaults, 'Values'))
        xlab = get(get(gca, 'XLabel'), 'String');
        ylab = get(get(gca, 'YLabel'), 'String');
        zlab = get(get(gca, 'ZLabel'), 'String');
    end
    
    % Speed up rendering a bit.
    plot3(cax,nan,nan,nan,p.Unmatched);
    set(gcf,'DoubleBuffer','off');
    set(cax,'XLimMode','manual','YLimMode','manual', ...
            'ZLimMode','manual','CLimMode','manual','ALimMode','manual');
    
    if ~any(strcmpi(p.UsingDefaults, 'Values'))
        xlabel(xlab);
        ylabel(ylab);
        zlabel(zlab);
    else 
        xlabel('k');
        ylabel('j');
        zlabel('i');
    end


    if all(cellfun(@isnumeric, transf))
        xlim(transf{3}([1 st(3)])+[-0.5 0.5]);
        ylim(transf{2}([1 st(2)])+[-0.5 0.5]);
        zlim(transf{1}([1 st(1)])+[-0.5 0.5]);
    else % all functions
        xlim(transf{3}([0.5 st(3)+0.5]));
	ylim(transf{2}([0.5 st(2)+0.5]));
	zlim(transf{1}([0.5 st(1)+0.5]));	      
    end
    
    set(cax,'YDir','reverse');
    set(cax,'ZDir','reverse');
    caxis([mn mx]);

    % Display grid.
    step = max(1,round(st/6));
    set(cax,'XTickMode','manual','YTickMode','manual','ZTickMode','manual');
    tk = transf{3}([1:step(3):st(3)-1 st(3)]);
    if length(tk) > 2 && tk(end)-tk(end-1) < .3*transf{3}(step(3))
        tk = transf{3}([tk(1:end-2) tk(end)]);
    end
    set(cax,'XTick',tk);
    tk = transf{2}([1:step(2):st(2)-1 st(2)]);
    if length(tk) > 2 && tk(end)-tk(end-1) < .3*transf{3}(step(2))
        tk = transf{2}([tk(1:end-2) tk(end)]);
    end
    set(cax,'YTick',tk);
    tk = transf{1}([1:step(1):st(1)-1 st(1)]);
    if length(tk) > 2 && tk(end)-tk(end-1) < .3*transf{3}(step(1))
        tk = transf{1}([tk(1:end-2) tk(end)]);
    end
    set(cax,'ZTick',tk);
    grid on;
    
    % Draw the voxels.
    if fast
        surface(transf{3}(st(3)*[0 1;0 1]+.5),...
                transf{2}(st(2)*[0 0;1 1]+.5),...
                transf{1}(.5*ones(2)), ...
                squeeze(T(1,:,:)),'AlphaData',adata(T(1,:,:)),options{:});
        surface(transf{3}(st(3)*[0 1;0 1]+.5), ...
                transf{2}(st(2)*[0 0;1 1]+.5), ...
                transf{1}((st(1)+.5)*ones(2)), ...
                squeeze(T(end,:,:)),'AlphaData', ...
            adata(T(end,:,:)),options{:});
        for n = 1:size(T,1)-1
            cdata = squeeze(max(T(n,:,:),T(n+1,:,:)));
            surface(transf{3}(st(3)*[0 1;0 1]+.5),...
                    transf{2}(st(2)*[0 0;1 1]+.5), ...
                    transf{1}((n+.5)*ones(2)),...
                    cdata,'AlphaData',adata(cdata),options{:});
        end
        surface(transf{3}(st(3)*[0 1;0 1]+.5), ...
                transf{2}(.5*ones(2)),...
                transf{1}(st(1)*[0 0;1 1]+.5), ...
                squeeze(T(:,1,:)),'AlphaData',adata(T(:,1,:)),options{:});
        surface(transf{3}(st(3)*[0 1;0 1]+.5), ...
                transf{2}((st(2)+.5)*ones(2)), ...
                transf{1}(st(1)*[0 0;1 1]+.5), ...
                squeeze(T(:,end,:)),'AlphaData', ...
                adata(T(:,end,:)),options{:});
        for n = 1:size(T,2)-1
            cdata = squeeze(max(T(:,n,:),T(:,n+1,:)));
            surface(transf{3}(st(3)*[0 1;0 1]+.5),...
                    transf{2}((n+.5)*ones(2)), ...
                    transf{1}(st(1)*[0 0;1 1]+.5), ...
                    cdata,'AlphaData',adata(cdata), options{:});
        end
        surface(transf{3}(.5*ones(2)),...
                transf{2}(st(2)*[0 1;0 1]+.5), ...
                transf{1}(st(1)*[0 0;1 1]+.5), ...
                squeeze(T(:,:,1)),'AlphaData',adata(T(:,:,1)),options{:});
        surface(transf{3}((st(3)+.5)*ones(2)), ...
                transf{2}(st(2)*[0 1;0 1]+.5), ... 
                transf{1}(st(1)*[0 0;1 1]+.5), ...
                squeeze(T(:,:,end)),'AlphaData', ...
                adata(T(:,:,end)),options{:});
        for n = 1:size(T,3)-1
            cdata = squeeze(max(T(:,:,n),T(:,:,n+1)));
            surface(transf{3}((n+.5)*ones(2)),...
                    transf{2}(st(2)*[0 1;0 1]+.5), ...
                    transf{1}(st(1)*[0 0;1 1]+.5), ...
                    cdata,'AlphaData',adata(cdata), ...
                options{:});
        end
    else
        opaq = alpha;
        opaq(alpha<t) = t^(1-k)*alpha(alpha<t).^k;
        opaq(opaq > 0.9) = 1; % Fix Matlab render bug.
        patch('Vertices',verts, ...
              'Faces',faces(opaq > cutoff,:), ...
              'EdgeColor','none','FaceAlpha','flat','FaceColor','flat', ...
              'FaceVertexAlphaData',opaq(opaq > cutoff), ...
              'FaceVertexCData',color(opaq > cutoff));
    end
    
end

end
