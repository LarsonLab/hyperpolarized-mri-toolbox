function [weights] = coil_dist_map(mask)
    
    maskSize = size(mask);
    weights = zeros(maskSize); 
    
    %tic
    for z=1:maskSize(3)
        mask_sl = squeeze(mask(:,:,z));
        
        % get outline/perim of mask
        mask_sl = imfill(bwmorph(bwareaopen(mask_sl,300),"fill"),"holes");
        %figure, imagesc(mask); axis off square;
        bw2 = bwperim(mask_sl);
        %figure, imagesc(bw2)
        
        % calculate weights based on distance from mask perim
        w = bwdist(bw2) .^0.5;
        weights(:,:,z) = (1 - (w ./max(w(:)))) .* mask_sl;
        %figure, imagesc(weights)
        
    end
    %toc

end