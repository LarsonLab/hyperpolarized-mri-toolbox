% script to create cardiac masks
function [masks, layered_masks] = create_cardiac_masks(zlevels)
    arguments
        zlevels (1,:) {mustBeNumeric}
    end
    
    [all, ~, ~, ~] = read_atlas("ED", false);
    connections = read_connectivity();
       
    LV_pts = 1:1500;
    RV_pts = 1501:3224;
    EPI_pts = 3225:5582;
    LV_faces = 1:3072;
    RV_faces = cat(2, 3073:4480, 4480:6752);
    EPI_faces = 6753:11616;

    dilation_radius = 1;
    slice_size = 150;

    layered_masks = zeros(slice_size, slice_size, numel(zlevels));
    for zlevel_idx = 1:numel(zlevels)
        zlevel = zlevels(zlevel_idx);
        LV_slice = create_single_slice(zlevel, slice_size, all, connections, LV_pts, LV_faces, dilation_radius);
        RV_slice = create_single_slice(zlevel, slice_size, all, connections, RV_pts, RV_faces, dilation_radius);
        EPI_slice = create_single_slice(zlevel, slice_size, all, connections, EPI_pts, EPI_faces, dilation_radius);

        % 1 = LV, 2 = RV, 3 = EPI
        combined_slice = EPI_slice * 3;
        combined_slice(LV_slice > 0) = 1;
        combined_slice(RV_slice > 0) = 2;
        layered_masks(:,:,zlevel_idx) = combined_slice;

    end
    % separate myocardium. 1 = lv, 2 = rv, 3 = lvmy, 4 = rvmy
    layered_masks = separate_myocardium(layered_masks);

    masks = zeros(slice_size, slice_size, numel(zlevels), 3);
    masks(:,:,:,1) = layered_masks == 1;
    masks(:,:,:,2) = layered_masks == 2;
    masks(:,:,:,3) = layered_masks == 3;
    masks(:,:,:,4) = layered_masks == 4;
end


function filled = create_single_slice(zlevel, slice_size, all, connections, pts, faces, dilation_radius)
    mesh_pts = trace_surface(zlevel, slice_size, all, connections, pts, faces);
    closed = close_small_holes(mesh_pts, dilation_radius);
    closed_fr = close_big_holes(closed);
    filled = fill_outline(closed_fr);
    fprintf("z slice at %i has been made\n", zlevel);
end



% ================
% == READ ATLAS ==
% ================
function [all, LV, RV, EPI] = read_atlas(stage, plot)  
    % ARGS:
    %   stage = either "ED" for end diastole or "ES" for end systole
    %   plot = logical flag to plot after reading
    arguments
        stage string {mustBeMember(stage, ["ED", "ES"])} = "ED"
        plot logical = false
    end

    pc = h5read('UKBRVLV.h5', '/COEFF'); % read the principal components
    ev = h5read('UKBRVLV.h5', '/LATENT'); % read the eigenvalues
    mu = h5read('UKBRVLV.h5', '/MU'); % read the mean shape
    
    % generate the first principal mode
    % with 1.5 times the standard deviation
    S = mu + (1.5 .* sqrt(ev(1)) .* pc(:,1)');
    
    % get ED & ES points, & convert to 3 columns matrix [x, y, z]
    N = length(S);
    if stage == "ED"
        all = reshape(S(1:N/2), 3, [])';
    else
        all = reshape(S((N/2+1):end), 3, [])';
    end
    
    % only take some of the points
    % these separations were from https://github.com/ComputationalPhysiology/ukb-atlas/blob/main/src/ukb/surface.py
    LV = all(1:1500, :);
    RV = all(cat(2, 1501:3224, 5730:5808), :);
    EPI = all(3225:5582, :);
    everything_else = all(5583:end, :);
    
    if plot
        % plot ED points in blue
        % plot ES points in red
        figure('Color', 'w');
        hold on;
        plot3(LV(:,1), LV(:,2), LV(:,3), 'r.');
        plot3(RV(:,1), RV(:,2), RV(:,3), 'g.');
        plot3(EPI(:,1), EPI(:,2), EPI(:,3), 'b.');
        plot3(everything_else(:,1), everything_else(:,2), everything_else(:,3), 'y.');
        
        axis vis3d
        axis equal
    end
end


% =======================
% == READ CONNECTIVITY ==
% =======================
function connectivity = read_connectivity(path)
    arguments
        path string = 'connectivity.txt'
    end
    connectivity = readmatrix(path) + 1; % +1 to account for the fact that this came from python, and python uses 0-indexing
end


% ===================
% == PROJECT EDGES ==
% ===================
function slice = trace_surface(zlevel, slice_size, pts, faces, pt_idxs, face_idxs)
    % ARGS:
    %   zlevel = z coordinate of slice
    %   slice_size = scalar of how big the slice should be
    %   pts = Nx3 array of points
    %   faces = Nx3 array of pt idxs that define triangles
    %   pt_idxs = 1xN array of idxs of all points to look at
    %   face_idxs = 1xN array of idxs of all faces to look at
    % OUT:
    %   slice = 2d matrix slice

    slice = zeros(slice_size);
    offset = slice_size / 2;
    
    for face_idx = face_idxs
        face = faces(face_idx, :);
        % if any point index in face isn't in pt_idxs, skip this face
        if any(~ismember(face, pt_idxs))
            continue
        end

        verts = pts(face, :);
        n_verts = size(verts, 2);

        % for each pair of vertices, find intersecting points
        for vert_idx = 1:n_verts
            vert1 = verts(vert_idx, :);
            vert2 = verts(mod(vert_idx, n_verts) + 1, :);
            if (vert1(3) >= zlevel && vert2(3) <= zlevel) || ...
                (vert1(3) <= zlevel && vert2(3) >= zlevel)
                intersection_coords = calc_intersection(vert1, vert2, zlevel);
                slice(round(intersection_coords(1)) + offset, round(intersection_coords(2)) + offset) = 1;
            end
        end
    end    
end


% ===============================================
% == CLOSE SMALL HOLES (morphological closing) ==
% ===============================================
function closed = close_small_holes(slice, dilation_radius)
    closed = imdilate(slice, strel('disk', dilation_radius));
    closed = bwskel(cast(closed, "logical"));
    
    % remove the weird bubble thing
    kernel = [0 1 0; 1 0 1; 0 1 0];
    n_adjacent = conv2(closed, kernel, 'same');
    bubbles = (n_adjacent == 4 & closed == 0);
    closed = closed + bubbles;

    % re-skeletonize
    closed = bwskel(cast(closed, 'logical'));
end


% =====================
% == close big holes ==
% =====================
function bresenhamified = close_big_holes(closed)
    bresenhamified = closed;
    % identify endpoints
    kernel = [1 1 1; 1 0 1; 1 1 1];
    surrounding_pixels = conv2(closed, kernel, 'same');
    
    endpoints = (surrounding_pixels == 1 & closed == 1);
    
    endpt_indices = find(endpoints);
    n_endpts = numel(endpt_indices);
    if n_endpts == 0 
        return
    end

    [endpt_x, endpt_y] = ind2sub(size(endpoints), endpt_indices);

    if mod(numel(endpt_x), 2) == 1
        error("there's a weird outcrop somewhere")
    end
    
    if numel(endpt_indices) > 2
        % multiple gaps case
        distances = zeros(n_endpts);
        for endpt_idx1 = 1:numel(endpt_indices)
            for endpt_idx2 = 1:numel(endpt_indices)
                % check if the points are connected
                pt1 = [endpt_x(endpt_idx1), endpt_y(endpt_idx1)];
                pt2 = [endpt_x(endpt_idx2), endpt_y(endpt_idx2)];
                is_connected = check_connection(closed, pt1, pt2);
                % if not connected, set it to distance between 2 points
                if ~is_connected
                    dx = pt1(1) - pt2(1);
                    dy = pt1(2) - pt2(2);
                    distances(endpt_idx1, endpt_idx2) = sqrt(dx.^2 + dy.^2);
                end
            end
        end

        remaining_endpts = 1:n_endpts;
        while any(distances > 0, 'all')
            % find minimum nonzero distance, fill that one in first
            nonzero_min = min(distances(distances > 0));
            [pt1_idx, pt2_idx] = ind2sub(size(distances), find(distances == nonzero_min, 1));
            bresenhamified = draw_line(bresenhamified, endpt_x(pt1_idx), endpt_y(pt1_idx), endpt_x(pt2_idx), endpt_y(pt2_idx));

            % set these guys to 0
            distances(pt1_idx,:) = 0;
            distances(:,pt1_idx) = 0;
            distances(pt2_idx,:) = 0;
            distances(:,pt2_idx) = 0;
            
            % recheck connectivity
            remaining_endpts(pt1_idx) = 0;
            remaining_endpts(pt2_idx) = 0;
            for endpt_idx1 = remaining_endpts(remaining_endpts > 0)
                for endpt_idx2 = remaining_endpts(remaining_endpts > 0)
                    pt1 = [endpt_x(endpt_idx1), endpt_y(endpt_idx1)];
                    pt2 = [endpt_x(endpt_idx2), endpt_y(endpt_idx2)];
                    is_connected = check_connection(bresenhamified, pt1, pt2); % todo: consider making a wrapper function for checking connectivity of a set of points
                    if is_connected
                        distances(endpt_idx1, endpt_idx2) = 0;
                    end
                end
            end
        end
    end
    
    % re-skeletonize (in case some corners were introduced and endpoints
    % can't be found)
    bresenhamified = bwskel(bresenhamified);
    surrounding_pixels = conv2(bresenhamified, kernel, 'same');
    endpoints = (surrounding_pixels == 1 & closed == 1);
    endpt_indices = find(endpoints);
    [endpt_x, endpt_y] = ind2sub(size(endpoints), endpt_indices);
    bresenhamified = draw_line(bresenhamified, endpt_x(1), endpt_y(1), endpt_x(2), endpt_y(2));
end


% ================
% == flood-fill ==
% ================
function filled = fill_outline(outline)
    inverse_fill = imfill(outline, [1,1]);
    filled = ~inverse_fill + outline;
end


% ======================
% == helper functions ==
% ======================
% quick wrapper function to draw a bresenham line
function i_new = draw_line(i, x1, y1, x2, y2)
    [x,y] = bresenham(x1, y1, x2, y2);
    line_indices = sub2ind(size(i), x, y);
    i_new = i;
    i_new(line_indices) = 1;
end

function xy = calc_intersection(a, b, z)
    % a and b are coordinates [x,y,z]
    % z is the plane to intersect
    xy = [...
        a(1) + (((z - a(3)) * (b(1) - a(1))) / (b(3) - a(3))), ...
        a(2) + (((z - a(3)) * (b(2) - a(2))) / (b(3) - a(3))) ...
    ];
end

% ============================
% == distinguish myocardium ==
% ============================
function new_mask = separate_myocardium(og_layered)
    % og_layered = row x col x slice, 1 = lv, 2 = rv, 3 = epi
    % outputs a new layered mask where 1 = lv, 2 = rv, 3 = lvmy, 4 = rvmy
    new_mask = og_layered;

    % find distances
    lv_dist = bwdist(new_mask == 1);
    imshow(lv_dist(:,:,20), [0, max(lv_dist, [], 'all')]);
    
    rv_dist = bwdist(new_mask == 2);
    imshow(rv_dist(:,:,20), [0, max(rv_dist, [], 'all')]);
    
    % for each voxel, see if a) it's epi, b) if if it's closer
    for row = 1:size(new_mask, 1)
        for col = 1:size(new_mask, 2)
            for slice = 1:size(new_mask, 3)
                if new_mask(row,col,slice) == 3 && rv_dist(row,col,slice) < lv_dist(row,col,slice)
                    new_mask(row,col,slice) = 4;
                end
            end
        end
    end
    
    imshow(new_mask(:,:,20), [0,4]);
end

%% run the script
clear;
zlevels = -40:2:50;
[masks, layered_masks] = create_cardiac_masks(zlevels);

matlab.io.savevariablestoscript('outputs/mask.m', {'masks'})

