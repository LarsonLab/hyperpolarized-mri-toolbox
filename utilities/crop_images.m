function I_cropped = crop_images(I,mask)
% crop images with a mask 
% Inputs:
%  mask - [xmin ymin width height]  or {ymin:ymax, xmin:xmax} 
% Outputs:
%  I_cropped - [height, width] or [ymax - ymin + 1, xmax - xmin + 1]
%
% Author: Shuyu Tang

I_cropped = [];

orig_sz = size(I);
if length(orig_sz) > 3
    I = reshape_to_3d(I);
end
if iscell(mask)
    % mask is {ymin:ymax, xmin:xmax} 
    crop_idx = mask;
    I_cropped = I(crop_idx{1},crop_idx{2},:);
else
    % mask is [xmin ymin width height]  
    for i = 1:size(I,3)
        I_cropped = cat(3,I_cropped, imcrop(I(:,:,i), mask));
    end
end

if length(orig_sz) > 3
    new_sz = orig_sz;
    new_sz(1) = size(I_cropped,1);
    new_sz(2) = size(I_cropped,2);    
    I_cropped = reshape(I_cropped,new_sz);
end