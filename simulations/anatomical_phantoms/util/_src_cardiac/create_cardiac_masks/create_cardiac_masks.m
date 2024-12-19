function [cardiac_masks, layered_masks] = create_cardiac_masks(vertices, tolerance, maskSize)
% CREATE_CARDIAC_MASKS creates cardiac masks from point cloud data for
% cardiac_metabolic_phantom
%   Parameters: 
%       vertices        = struct with fields mc, lv, rv, each of which are Nx3
%                         matrices of xyz coordinates
%       tolerance       = how "off" a point can be to be included in a
%                         slice. defaults to 1.
%       maskSize        = how large the masks should be. [xmin, xmax, xres;
%                         ymin, ymax, yres; zmin, zmax, zres;]. default =
%                         size chosen based on point cloud data
% 
%   Outputs:
%       cardiac_masks   = 4d array of masks (nx, ny, nz, tissue)
%       layered_masks   = 3d array of unseparated masks

    
    arguments 
        vertices struct 
        tolerance (1,1) double = 1
        maskSize (3,3) double = zeros(3)
    end
    
    % default maskSize
    if ~any(maskSize,'all')
        all_vertices = cat(1, vertices.mc, vertices.lv, vertices.rv);
        % get reasonable "window size"
        x_range = max(all_vertices(:,1),[],'all') - min(all_vertices(:,1),[],'all');
        y_range = max(all_vertices(:,2),[],'all') - min(all_vertices(:,2),[],'all');
        % get x y bounds
        ranges = [x_range, y_range];
        for dim = 1:numel(ranges)
            offset = 0.25 * ranges(1);
            maskSize(dim,1) = min(all_vertices(:,dim),[],'all') - offset;
            maskSize(dim,2) = max(all_vertices(:,dim),[],'all') + offset;
            maskSize(dim,3) = round(maskSize(dim,2) - maskSize(dim,1) + 1);
        end
        % get z bounds
        maskSize(3,1) = min(all_vertices(:,3),[],'all');
        maskSize(3,2) = max(all_vertices(:,3),[],'all');
        maskSize(3,3) = round(maskSize(3,2) - maskSize(3,1) + 1);
    end


    % read points into 3d array
    point_map = zeros(maskSize(2,3), maskSize(1,3), maskSize(3,3), 3);
    tissues = fieldnames(vertices);
    
    x = linspace(maskSize(1,1),maskSize(1,2),maskSize(1,3));
    y = linspace(maskSize(2,1),maskSize(2,2),maskSize(2,3));
    z = linspace(maskSize(3,1),maskSize(3,2),maskSize(3,3));
    [X,Y,Z] = meshgrid(x,y,z);

    for tissue = 1:numel(tissues)
        tissue_name = tissues{tissue}; % TODO: tissue vs tissue name is bad
        tissue_vertices = vertices.(tissue_name);
        tissue_point_map = zeros(size(X)); 
        
        % plot each point on the map
        for point = 1:size(tissue_vertices,1)
            % find distance from each point from the map
            xDiff = abs(X - tissue_vertices(point,1));
            yDiff = abs(Y - tissue_vertices(point,2));
            zDiff = abs(Z - tissue_vertices(point,3));
            totalDiff = xDiff + yDiff + zDiff;
            % add points within tolerance to map
            tissue_point_map(totalDiff <= tolerance) = 1;
        end
        point_map(:,:,:,tissue) = tissue_point_map;
    end


    % create convex hulls of each slice
    conv_hulls = zeros(size(point_map));
    slices = size(point_map, 3);
    for tissue = 1:numel(tissues)
        for slice = 1:slices
            conv_hulls(:,:,slice,tissue) = bwconvhull(point_map(:,:,slice,tissue));
        end
    end

    % layer hulls
    % mc = 1, lv = 2, rv = 3
    layered_masks = int8(conv_hulls(:,:,:,1)); % mc
    layered_masks(conv_hulls(:,:,:,2) == 1) = 2; % lv
    layered_masks(conv_hulls(:,:,:,3) == 1) = 3; % rv

    % separate rv and lv myocardium
    old_layered_masks = zeros(size(layered_masks));
    while ~isequal(old_layered_masks, layered_masks)
        old_layered_masks = layered_masks;
        for z = 1:size(layered_masks, 3)
            % get maps rv lv and mc maps
            slice = layered_masks(:,:,z);
            mc_map = slice == 1;
            lv_map = (slice == 2) | (slice == 4);
            rv_map = (slice == 3) | (slice == 5);
            kernel = [1,1,1 ; 1,0,1 ; 1,1,1];
            
            % get pixels around the lv/rv
            lv_conv_result = conv2(lv_map, kernel, "same");
            rv_conv_result = conv2(rv_map,kernel,"same");
            lv_mc = (lv_conv_result > 0) & mc_map;
            rv_mc = (rv_conv_result > 0) & mc_map;
            
            % update layered masks 
            new_slice = slice;
            new_slice(lv_mc == 1) = 4;
            new_slice(rv_mc == 1) = 5;
            layered_masks(:,:,z) = new_slice;
        end
    end

    % assume non-labeled myocardium is rv myocardium (TODO: maybe add more convolutions along the z axis)
    layered_masks(layered_masks == 1) = 4;
    
    layered_masks = layered_masks - 1;
    layered_masks(layered_masks < 0) = 0;
    
    % separate masks
    cardiac_masks = zeros(cat(2,size(layered_masks),4));
    cardiac_masks(:,:,:,1) = layered_masks == 1;
    cardiac_masks(:,:,:,2) = layered_masks == 2;
    cardiac_masks(:,:,:,3) = layered_masks == 3;
    cardiac_masks(:,:,:,4) = layered_masks == 4;

end