% FERMI FILTER FUNCTION
function [filtered_image] = fermi_filter(image, padscale, dim)
% -------------------------------------------------------------------------
% this function performs a 2d fft on the image I
% then applies a 2D fermi filter using the function fermi
% finally it applies an inverse 2d fft and returns the
% filtered image as proj.
%
% Inputs:
%   padscale - size of filtered image / size of image
%   dim - filter dimension
%
% JHH 7/29/02
% Shuyu Tang added dim > 1 for fermi_filtering in any dimension 6/18/19

if nargin < 3
    dim = 2; % fermi 2d
end


if padscale > 1
    if dim == 2
    szim = size(image);
    if length(szim) > 3
        image = reshape_to_3d(image);
    end        
    % *****************************************************
    % perform fft
    % *****************************************************
   
    I_kspace = fftshift(fftshift( fft2(fftshift(fftshift(image,1),2)) ,1),2);


    % ********************************************************
    %zeropad the image following algorithm Sean Fain's code from radonfft2.m
    % ********************************************************
    [yi,xi,n_slice] = size(I_kspace);
     
    %Ipad = zeros(padscale*xi,padscale*yi);
    Ipad = zeros(padscale*yi,padscale*xi,n_slice);
    
    Nx = padscale*xi;
    Ny = padscale*yi;

    %Ipad((Nx-xi)/2+1:Nx-(Nx-xi)/2,(Ny-yi)/2+1:Ny-(Ny-yi)/2) = I_kspace;
     Ipad((Ny-yi)/2+1:Ny-(Ny-yi)/2,(Nx-xi)/2+1:Nx-(Nx-xi)/2,:) = I_kspace;
   
     
     % ****************************************************
    % get fermi filter
    % ****************************************************
    [filter] = fermi(Nx -1, Ny -1, 0.5, 2);
     %[filter] = fermi(Ny -1, Nx -1, 0.5, 2);
     for kk = 1:n_slice
        dummy(:,:,kk) = filter;
     end
     clear filter; filter = dummy; clear dummy;
    % ****************************************************
    % apply filter to k-space data
    % ****************************************************
    %
     filtered_transform = Ipad .* filter ;
%     for( i = 1 : Nx )    
%         for(j = 1: Nx)
%             filtered_transform(i,j) = Ipad(i,j) * filter(i,j) ;
%         end;
%     end;

    % **************************************************
    % apply 2D inverse fft
    % **************************************************
    filtered_image = fftshift(fftshift( ifft2(fftshift(fftshift(filtered_transform,1),2)) ,1),2);
    else % non 2d filter
        if dim > 1
        old_order = 1:ndims(image);
        new_order = [dim,old_order(old_order~=dim)];
        image = permute(image,new_order);
        [~,restore_order] = ismember(old_order, new_order);
        end
        I_kspace = fftc(image,1);
        siz_Ik = size(I_kspace);
        xi = siz_Ik(1);
        Nx = padscale*xi;
        Ipad = padarray(I_kspace,(Nx-xi)/2);
        [filter] = fermi(Nx -1, 0,0.5, 1);
        filter = repmat(filter(:),[1,siz_Ik(2:end)]);
        filtered_transform = Ipad .* filter ;
        filtered_image = ifftc(filtered_transform,1);
        if dim > 1
            filtered_image = permute(filtered_image,restore_order);
        end
    end
else
    filtered_image = image;
end

if dim == 2
    filtered_image = real(filtered_image*padscale^2);
    szim(1) = size(filtered_image,1);
    szim(2) = size(filtered_image,2);
    filtered_image = reshape(filtered_image,szim);    
else
    filtered_image = real(filtered_image*padscale);
end

end

function [f] = fermi(N1, N2, cutoff, dim)
% -------------------------------------------------------------------------
% Fermi Filter Window
%
% Sample call: f = fermi(N,cutoff, beta, dim)
% N1 = length in longest dimension
% N2 = length in shortest dimension
% xres = pixel width along N1
% yres = pixel width along N2
% cutoff = fraction of N/2 in passband
% beta = severity of stopband; smaller beta, more severe
% the fall off. beta = 0.07813 * N/2 is GE default.
% dim = 1D returns fermi vector
% dim =2D returns fermi radial symmmetric matrix
%
% JHH 7/29/02 (modified from SBF 2-6-99)
% RO 3-1-08, abbrieviated for volumeviewer call

% if 1; (N1 >= N2);
x = -N1/2:N1/2;
if dim == 2    
    y = -N2/2:N2/2;
end
    %size(x)
% else
%     x = -N1/2 + 1 : N2/2;
%     x(find(x == 0)) = [];
% end

%y = x;
a = cutoff*max(x);
beta= 0.07813 * max(x);
if dim == 1 %1d filter
    r = sqrt(x.^2);
else
    [X,Y] = meshgrid(x,y);
    r = sqrt(X.^2 + Y.^2);
end
f = 1./(1+exp((r-a)./beta));
%f = resample(f,N2,N1);
end
