function imagescn(I, scale, dims, FigureWidth)
% function imagescn(I,[min max],[rows cols],FigureWidth,timedimension)
%
% function to display multiple images
%	I can be 2-d, i.e., single image
%            3-d,       array of images
%            4-d,       2-d array of images
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
% usage:  imagescn(I)
%         imagescn(I,[],[],[],[])
%         imagescn(I,scale)
%         imagescn(I,[],[rows cols])
%         imagescn(I,scale,[rows cols])
%         imagescn(I,[],[],[],timedimension)
%         ...

% written by: 	Peter Kellman  (kellman@nih.gov)
%				Laboratory for Cardiac Energetics
%				NIH NHBI
%   tools by: 	Dan Herzka  (herzkad@nhlbi.nih.gov)
%				Laboratory for Cardiac Energetics
%				NIH NHBI
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
% 	if isempty(scale)
% 		imagesc(I(:,:,k)); axis image; axis off;
%     elseif size(scale,1)==1
% 		imagesc(I(:,:,k),scale); axis image;axis off
%     elseif size(scale,1)==min(N,N1*N2);
%         imagesc(I(:,:,k),scale(k,:)); axis image;axis off
% 	end
	if isempty(scale)
		if Nd>=4
            imagesc(I(:,:,i0,j0)); axis image; axis off;
		else 		
            imagesc(I(:,:,k)); axis image; axis off;
		end
	elseif size(scale,1)==1
		if Nd>=4
            imagesc(I(:,:,i0,j0),scale); axis image; axis off;
		else		
            imagesc(I(:,:,k),scale); axis image;axis off
		end
	elseif size(scale,1)==min(N,N1*N2);
		if Nd>=4
            imagesc(I(:,:,i0,j0),scale(k,:)); axis image; axis off;
		else		
            imagesc(I(:,:,k),scale(k,:)); axis image;axis off
		end
    end
end

% set(fig_handle, 'PaperPosition', [0 0 FigureWidth FigureHeight]); % old
set(fig_handle, 'PaperPosition', [1 1 FigureWidth+1 FigureHeight+1]);

%toolbox;
colormap(gray)