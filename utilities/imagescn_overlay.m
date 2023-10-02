function imagescn_overlay(base, basescale, I, scale, dims, threshold, transparency, overlay_colormap, FigureWidth, colorbar_on)
% imagescn_overlay(base, basescale, I, scale, dims, threshold, transparency, overlay_colormap, FigureWidth)
% 
% function to overlay I on base images
% base and I can be 2-d, i.e., single image
%            3-d,       array of images
%            4-d,       2-d array of images
%
% allow base 3-d while I 4-d, assume 3rd dim is slice, 4th dim is time
%
% basescale and scale are images scales for base and I
%
% dims: rows and cols of displayed images
%
% threshold, transparency: see color_overlay
%
% overlay_colormap: colormap for overlay images, colormap for base images
% is gray
%
% user specified scale [min max] or is applied to all images (i.e., clims)
% or scale=[min1 max1; min2 max2; ...] will be applied to images 1,2,...,etc.
% defaults to independent scaling of each sub-image to full range (for scale=[])
%
% sub-images are displayed in N1 rows by N2 columns
%   N1 and N2 may be specified in function input as "rows" and "cols"
%		if input values for rows & cols input are such that
%		rows*cols < full number of images, then the 1-st rows*cols images are displayed
%	or N1 and N2 may be calculated by default:
%		for 3-d array (N1 and N2 are calculated to make approx square image array)
%		for 4-d array (N1 and N2 are 3rd and 4th array dimensions, respectively)
%
% FigureWidth sets the width in inches (defaults to 6 inches). It also sets the paper width
% so that exported figures (to files) will have the specified width.
%
% usage: 
%
% imagescn written by: 	Peter Kellman  (kellman@nih.gov)
%				Laboratory for Cardiac Energetics
%				NIH NHBI
% imagescn  tools by: 	Dan Herzka  (herzkad@nhlbi.nih.gov)
%				Laboratory for Cardiac Energetics
%				NIH NHBI
%   
% remove movie feature and add color_overlay by: Shuyu Tang (shuyutang2012@gmail.com)
%
Nd=ndims(I);
if Nd==2 % case of single image
    N=1;
	N1=1; N2=1;
elseif Nd==3 % case of array of images
    N=size(I,3);
    N2=ceil(sqrt(N)); N1=ceil(N/N2);
elseif Nd==4 % case of 2-d array of images
    N1=size(I,3);
    N2=size(I,4);
    N=N1*N2;
end
if exist('dims','var')
	if length(dims)==2; rows=dims(1);cols=dims(2);
		N1=rows;N2=cols;
	else
		if ~isempty(dims);disp('Error: must enter [rows cols] for dimensions'); return;end
	end
end
if ~exist('scale','var'); scale=[];end
if 1<size(scale,1) && size(scale,1) < min(N,N1*N2)
    disp('scale vector must be either: empty, size 1x2, i.e. [min max],')
    disp('or a matrix [min1 max1;min2 max2;...] with # rows equal to number of plots');
    return
end

if ~exist('colorbar_on')
    colorbar_on = 0;
end
if ~exist('overlay_colormap')
    overlay_colormap = 'jet';
end

set(0,'Units','Inches');
scnsize=get(0,'ScreenSize'); % [left,bottom,width,height] % full screen size available
ScreenWidth=scnsize(3); ScreenHeight=scnsize(4); % width & height in inches
Xsize=size(I,2);Ysize=size(I,1); % size of pixels in image (Xsize x Ysize)
border_percent=.01;
deltaX =(border_percent*Xsize); deltaY=(border_percent*Ysize); % calculate pixel borders as % of size
X0size=N2*Xsize+(N2-1)*deltaX; Y0size=N1*Ysize+(N1-1)*deltaY; % full figure size in pixels (before display)
aspect_ratio=Y0size/X0size; % aspect ratio

% center figure on screen with specified figure width in inches
if ~exist('FigureWidth'); FigureWidth=6;end  % figure width in inches (default is 6")
if isempty(FigureWidth); FigureWidth=6;end
FigureHeight=FigureWidth*aspect_ratio; % figure height in inches
FigureBottom=(ScreenHeight-FigureHeight)/2;
FigureLeft=(ScreenWidth-FigureWidth)/2;
fig_handle=gcf;%figure;
set(fig_handle,'Units','Inches')
set(fig_handle,'Position',[FigureLeft FigureBottom FigureWidth FigureHeight])

% calculate sub-image dimensions in inches
SubImageWidth=FigureWidth*Xsize/X0size;
SubImageHeight=FigureHeight*Ysize/Y0size;
Xborder=FigureWidth*deltaX/X0size;
Yborder=FigureHeight*deltaY/Y0size;

% set background color to be white
set(fig_handle,'Color',[1 1 1]);

% temporary fix for v7 pre-release bug
% v=version;
% if v(1)=='7';
%     set(fig_handle,'Color','none','Renderer','Painters','InvertHardCopy','off');
% end

% calculate sub-image dimensions in normalized units
SubImageWidth=SubImageWidth/FigureWidth;
SubImageHeight=SubImageHeight/FigureHeight;
Xborder=Xborder/FigureWidth;
Yborder=Yborder/FigureHeight;

Nd_base = ndims(base);

for k=1:min(N,N1*N2)
	i=ceil(k/N2);
	j=k-(i-1)*N2;
    if Nd>3
    	i0=ceil(k/size(I,4));
		j0=k-(i0-1)*size(I,4);
    end        
	SubImageLeft=(j-1)*(SubImageWidth+Xborder);
	SubImageBottom=(N1-i)*(SubImageHeight+Yborder);
	subplot('Position',[SubImageLeft SubImageBottom SubImageWidth SubImageHeight])


    if Nd>=4   
        %imagesc(I(:,:,i0,j0)); axis image; axis off;
        overlay = I(:,:,i0,j0);
    else 		
        %imagesc(I(:,:,k)); axis image; axis off;
        overlay = I(:,:,k);
    end

    if Nd_base>=4   
        base_I = base(:,:,i0,j0);
    elseif Nd_base == 3 && Nd == 3
        base_I = base(:,:,k);
    elseif Nd_base == 3 && Nd == 4
        base_I = base(:,:,i0);
    else
        base_I = base;
    end    
    
    
	if isempty(basescale)
        base_low = [];
        base_high = [];
	elseif size(basescale,1)==1
        base_low = basescale(1);
        base_high = basescale(2);
	elseif size(basescale,1)==min(N,N1*N2);
        base_low = basescale(k,1);
        base_high = basescale(k,2);
    end    
    
	if isempty(scale)
        overlay_low = [];
        overlay_high = [];
	elseif size(scale,1)==1
        overlay_low = scale(1);
        overlay_high = scale(2);
	elseif size(scale,1)==min(N,N1*N2);
        overlay_low = scale(k,1);
        overlay_high = scale(k,2);
    end
    
    color_overlay(base_I.', overlay.', ...
        base_low, base_high, overlay_low, overlay_high, ...
        threshold, transparency, overlay_colormap, colorbar_on);
    %colorbar off; axis image; axis off;

end

% set(fig_handle, 'PaperPosition', [0 0 FigureWidth FigureHeight]); % old
set(fig_handle, 'PaperPosition', [1 1 FigureWidth+1 FigureHeight+1]);

%toolbox;
%colormap(gray)
end

function color_overlay(base_slice, overlay_slice, ...
		       base_low, base_high, ...
		       overlay_low, overlay_high, ...
		       threshold, transparency,overlay_cmap, colorbar_on);
% COLOR_OVERLAY Plots a transparent color overlay of one image on another
%
%   COLOR_OVERLAY(BASE_SLICE, OVERLAY_SLICE) is the simplest way to
%   use it, in which the overlay_slice is plotted as an overlay
%   over the base_slice.  
%
%   **  It is YOUR responsibility as the user to ensure that the   **
%   **  base image and overlay image are the same size!  Otherwise **
%   **  you will get a warning and some (probably) weird looking   **
%   **  image results.  So do any interpolation BEFOREHAND!        **
%
%
%   COLOR_OVERLAY(BASE_SLICE, OVERLAY_SLICE, BASE_LOW, BASE_HIGH,
%   OVERLAY_LOW, OVERLAY_HIGH, THRESHOLD, TRANSPARENCY, COLORBAR) 
%   is the full usage, with other parameters described here.  You can 
%   always leave them out for a default value, or enter a [] instead of 
%   a number.
%
%     BASE_LOW / BASE_HIGH  :      Allows you to window/level your
%                                  anatomical image by specifying
%                                  minimum and maximum values for the
%                                  base_slice.  Default is the full
%                                  range of pixel values.
%
%     OVERLAY_LOW / OVERLAY_HIGH : Specifies the upper and lower
%                                  ends of the overlaid image.
%                                  Anything higher than the upper
%                                  bound will be displayed as just
%                                  red, and everything lower will
%                                  be just blue.  Default is the
%                                  full range of pixel values.
%
%     THRESHOLD :                  The minimum absolute value in
%                                  the overlaid image that will
%                                  show up.  If you set this to
%                                  zero or leave it blank, your
%                                  background around the image will
%                                  show up as some sort of color.
%                                  Setting this to, say, 1, means
%                                  that any pixels with a value
%                                  between -1 and 1 will not be
%                                  displayed as an overlay at all.
%
%     TRANSPARENCY :               Alpha value between 0 and 1,
%                                  with 1 being solid and 0 being
%                                  invisible.  Default is 0.3.
%     
%
%        Michael C. Lee, Ph.D.
%        Department of Radiology
%        University of California, San Francisco
%
%        Last modified 30 March 2004
%

if (nargin < 2) 
  %help color_overlay;
  return; 
end;

if (nargin < 3 | isempty(base_low))
  base_low = min(base_slice(:));
end

if (nargin < 4 | isempty(base_high))
  base_high = max(base_slice(:));
end

if (nargin < 5 | isempty(overlay_low))
  overlay_low = min(overlay_slice(:));
end

if (nargin < 6 | isempty(overlay_high))
  overlay_high = max(overlay_slice(:));
end

if (nargin < 7 | isempty(threshold))
  threshold = min(abs(overlay_slice(:))) + 0.001;
end

if (nargin < 8 | isempty(transparency))
  transparency = 0.7;
end

if (nargin < 9 | isempty(overlay_cmap))
  overlay_cmap = 'jet';
end

if (nargin < 10 | isempty(colorbar_on))
  colorbar_on = 0;
end


if (size(base_slice) ~= size(overlay_slice))
  warning(sprintf('Base slice is %d x %d  ... overlay sice is %d x %d', ...
		  size(base_slice), size(overlay_slice)));
end
  
set(gcf,'Renderer','OpenGL','RendererMode','manual');

if (length(threshold) == 1)
  alphadata = transparency .* double(abs(overlay_slice) > threshold);
elseif (length(threshold) == 2)
  alphadata = transparency .* double(overlay_slice < threshold(1) | ...
                                     overlay_slice > threshold(2));
end

% ------------------------------------------------------------------------

% Matlab only allows one colormap per image, so what we end up
% having to do is to concatenate the grayscale and jet colormaps,
% and rescale the data to make sure that the overlay lives in the
% jet regime and the anatomical image lives in the grayscale regime

% The image will live in the range -full_range:0 and the
% overlay will be 0:fullrange

% Clamp down to the appropriate range

base_slice    = max(base_slice, base_low);
base_slice    = min(base_slice, base_high);

overlay_slice = max(overlay_slice, overlay_low);
overlay_slice = min(overlay_slice, overlay_high);

full_range    = max(base_slice(:)) - min(base_slice(:));

base_slice    = ...
    base_slice - max(base_slice(:)) - 0.1;
overlay_slice = ...
    full_range*(overlay_slice-overlay_low) / ...
    (overlay_high-overlay_low) + 0.1;

% Concatenate colormaps

% colormap jet(512);  cmapjet = colormap;
% colormap autumn(512); cmapautumn = colormap;
colormap gray(512); cmapgray = colormap;
% colormap cool(512); cmapcool = colormap;
% colormap hsv(512); cmaphsv = colormap;
% colormap hot(512); cmaphot = colormap;
colormap(overlay_cmap); cmapoverlay = colormap;
cmap = [cmapgray;cmapoverlay];
%cmap = [cmapgray;cmapautumn];

% ------------------------------------------------------------------------

% Plot the anatomical image

pcolor(base_slice);
shading flat; 
view([90 90]); 
axis equal; 
axis tight;

hold on;

% ------------------------------------------------------------------------

% Plot the overlay image

overlayh = pcolor(overlay_slice); shading flat;

% Set transparency data

set(overlayh, 'AlphaData', alphadata, ...
	      'AlphaDataMapping', 'none', ...
	      'FaceAlpha','flat');

colormap(cmap);
caxis([-full_range-0.1, full_range+0.1]);

% MCL : Added 30 March 2004 ... allows a correct color bar

if (colorbar_on == 1)
  h = colorbar;
  set(gcf,'Renderer','OpenGL','RendererMode','manual', ...
          'BackingStore','off');
  old_range = get(h,'YLim');
  new_range = [overlay_low overlay_high];
  label_pts = [0:(old_range(2)/5):old_range(2)];
  %disp(overlay_high);
  set(h, 'YLim', [0 old_range(2)]);
  set(h, 'YTickLabelMode', 'manual', 'YTickMode', 'manual');
  set(h, 'YTick', label_pts);
  for (i = 1:length(label_pts))
    label_names{i} = sprintf('%.1f', ...
   new_range(1)+(new_range(2)-new_range(1))/5*(i-1));
  end
  set(h, 'YTickLabel',label_names');
  set(h, 'Box', 'on');
end

axis off;

hold off;
end