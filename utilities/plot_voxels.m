function plot_voxels(S, f, xsize, ysize, Sscale);
% plot_voxels(S, f, xsize, ysize, Smax)
% 
% Plots voxel data in individual plots within a figure.
%
% S - data, dimensions: 1=y, 2=x, 3=freq (or time)
% f - frequency or time axis
% xsize - plots voxels within the range of xsize(1) to xsize(2)
% ysize - plots voxels within the range of ysize(1) to ysize(2)
% S - amplitudes in plots
%
% Peder Larson, 7/2008

%S = squeeze(S);

if (nargin < 2) || isempty(f)
    f = 1:size(S,3);
end

if (nargin < 3) || isempty(xsize)
    xsize = [1 size(S,2)];
end

if (nargin < 4) || isempty(ysize)
    ysize = [1 size(S,1)];
end

if nargin < 5
    Sscale = [min(S(:)) max(S(:))];
end


Nx = xsize(2)-xsize(1) + 1;
Ny = ysize(2)-ysize(1) + 1;

sp_x = 1/Nx;
sp_y = 1/Ny;

clf;

for Ix = 1:Nx
    for Iy = 1:Ny
%        subplot( Ny,Nx, Ix + Nx*(Iy-1) ) 
        subplot('Position', [(Ix-1)*sp_x, (Ny-Iy)*sp_y, sp_x, sp_y] ) 
        plot(f, squeeze(S( Iy +ysize(1)-1, Ix + xsize(1)-1,:)))
        axis tight
        ylim(Sscale)
        set(gca,'xtick',[],'ytick',[])

    end
end
