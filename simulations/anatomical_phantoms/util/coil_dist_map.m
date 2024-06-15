function [weights] = coil_dist_map(mask)
    
    maskSize = size(mask);
    weights = zeros(maskSize); 

    % create y gradient
    x = linspace(-1, 1, maskSize(1));
    y = linspace(-1, 1, maskSize(2));
    [~, Y] = meshgrid(x, y);
    lim = [0.6, 1.2];
    grad = 0.5*(lim(2) - lim(1))*Y + 0.5*(lim(1) + lim(2));
    
    %tic
    for z=1:maskSize(3)
        mask_sl = squeeze(mask(:,:,z));

        
        
        % get outline/perim of mask
        mask_sl = imfill(bwmorph(bwareaopen(mask_sl,300),"fill"),"holes");
        %figure, imagesc(mask); axis off square;
        bw2 = bwperim(mask_sl);
        %figure, imagesc(bw2)

        %reverse_mask
        mask_rev = imcomplement(mask_sl);
        
        % calculate weights based on distance from mask perim
        w = bwdist(bw2) .^0.5;
        weights(:,:,z) = ((1 - (w ./max(w(:)))) .* grad .* mask_sl) + mask_rev;
        %figure, imagesc(weights)
        
    end
    %toc

end