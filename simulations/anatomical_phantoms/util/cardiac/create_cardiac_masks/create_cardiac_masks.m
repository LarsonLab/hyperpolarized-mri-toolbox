% ================
% == READ ATLAS ==
% ================
function [pts, LV, RV, EPI] = read_atlas(stage, plot)  
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

    deviations = rand([200,1]); % how many standard deviations to vary each pc
    scaled_ev = deviations .* sqrt(ev); % scaled eigenvalues
    variations = repmat(scaled_ev', [size(pc, 1), 1]) .* pc;
    summed_variations = sum(variations,2)';
    
    % generate the first principal mode
    % with 1.5 times the standard deviation
    S = mu + summed_variations;
    
    % get ED & ES points, & convert to 3 columns matrix [x, y, z]
    N = length(S);
    if stage == "ED"
        pts = reshape(S(1:N/2), 3, [])';
    else
        pts = reshape(S((N/2+1):end), 3, [])';
    end
    
    % only take some of the points
    % these separations were from https://github.com/ComputationalPhysiology/ukb-atlas/blob/main/src/ukb/surface.py
    % and i kinda trust them more than my own lol
    LV = pts(1:1500, :);
    RV = pts(cat(2, 1501:3224, 5730:5808), :);
    EPI = pts(3225:5582, :);
    everything_else = pts(5583:end, :);
    
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
    % tbh, that's the main reason why i wanted a wrapper function for this
    % lol
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
    %   slice = 2d matrix slice with points and (most) lines

    slice = zeros(slice_size);
    offset = slice_size / 2;
    
    for face_idx = face_idxs
        face = faces(face_idx, :);
        % if any point index in face isn't in pt_idxs, skip this face
        % because i guess there's points in connectivity.txt that don't exist/aren't relevant to the parts we're looking at
        if any(~ismember(face, pt_idxs))
            continue
        end

        verts = pts(face, :);
        n_verts = size(verts, 2); % i mean, this is gonna be 3 because everything is a triangle, but... NO MAGIC NUMBERS

        % for each pair of vertices, find intersecting points
        intersection_coords = {};
        for vert_idx = 1:n_verts
            vert1 = verts(vert_idx, :);
            vert2 = verts(mod(vert_idx, n_verts) + 1, :);
            % funny thing about this line: it's supposed to be something
            % like vert2 = verts(mod(vert_idx + 1, 3), :), but because
            % matlab uses 1-indexing, it turns into 
            % mod((vert_idx + 1) - 1, n_verts) + 1, so that's what that does
            if (vert1(3) >= zlevel && vert2(3) <= zlevel) || ...
                (vert1(3) <= zlevel && vert2(3) >= zlevel)
                xy = find_plane_intersection(vert1, vert2, zlevel);
                intersection_coords{numel(intersection_coords) + 1} = [round(xy(1)) + offset, round(xy(2)) + offset];
                slice(round(xy(1)) + offset, round(xy(2)) + offset) = 1; % yes, it's code duplication, but i'm lazy
            end

            % if there's 2 intersection coordinates, draw a line between them
            if numel(intersection_coords) == 2
                x1 = intersection_coords{1}(1);
                y1 = intersection_coords{1}(2);
                x2 = intersection_coords{2}(1);
                y2 = intersection_coords{2}(2);
                slice = draw_line(slice, x1, y1, x2, y2);
            end
        end
    end    
end

% ======================
% == HELPER FUNCTIONS ==
% ======================
function xy = find_plane_intersection(a, b, z)
    % a and b are coordinates [x,y,z]
    % z is the plane at z=z
    % xy are the xy coordinates where line ab intersects z
    xy = [...
        a(1) + (((z - a(3)) * (b(1) - a(1))) / (b(3) - a(3))), ...
        a(2) + (((z - a(3)) * (b(2) - a(2))) / (b(3) - a(3))) ...
    ];
end

% quick wrapper function to draw a bresenham line
function I_new = draw_line(I, x1, y1, x2, y2)
    [x,y] = bresenham(x1, y1, x2, y2);
    line_indices = sub2ind(size(I), x, y);
    I_new = I;
    I_new(line_indices) = 1;
end

function [new_endpt_x, new_endpt_y] = remove_erroneous_endpoints(I, endpt_x, endpt_y)
    % removes pixels that could be considered offshoots
    % this is a very "this just needs to work" kind of function
    kernel = ones(5);
    n_neighbors = conv2(I, kernel, 'same'); % includes the corner itself

    new_endpt_x = [];
    new_endpt_y = [];
    for endpt_idx = 1:numel(endpt_x)
        if n_neighbors(endpt_x(endpt_idx), endpt_y(endpt_idx)) <= 3
            new_endpt_x(numel(new_endpt_x) + 1, 1) = endpt_x(endpt_idx);
            new_endpt_y(numel(new_endpt_y) + 1, 1) = endpt_y(endpt_idx);
        else
            disp("removed erroneous endpoint")
        end
    end
end

function I = remove_blobs(I, blob_size)
    % removes small blobs of isolated pixels. intended to work with even
    kernel = ones(blob_size + 2) .* 100;
    kernel(2:(blob_size + 1), 2:(blob_size + 1)) = 1;
    conved = conv2(I, kernel, 'same');
    blob_locations = conved < 100 & conved ~= 0;
    if mod(blob_size, 2) == 0
        new_kernel = zeros(blob_size + 1);
        new_kernel(2:end, 2:end) = 1; % see, kernels with even dimensions are weird, so this manipulates the anchor of the kernel so it actually catches the right pixels
    else
        new_kernel = ones(blob_size);
    end
    to_remove = conv2(blob_locations, new_kernel, 'same');
    I(to_remove ~= 0) = 0;
end



% =====================
% == CLOSE BIG HOLES ==
% =====================
function bresenhamified = close_big_holes(I)
    I = bwskel(logical(I));
    I = remove_blobs(I, 5); % i think 5 is reasonable
    bresenhamified = I;
    % identify endpoints
    kernel = [1 1 1; 1 0 1; 1 1 1];
    surrounding_pixels = conv2(I, kernel, 'same');
    endpoints = (surrounding_pixels == 1 & I == 1);
    endpt_indices = find(endpoints);
    n_endpts = numel(endpt_indices);
    if n_endpts == 0 
        return
    end

    [endpt_x, endpt_y] = ind2sub(size(endpoints), endpt_indices);

    if mod(numel(endpt_x), 2) == 1
        figure; imshow(I);
        error("there's a weird outcrop somewhere idk help")
    end
    
    % multiple gaps case
    if numel(endpt_indices) > 2
        distances = zeros(n_endpts);
        for endpt_idx1 = 1:numel(endpt_indices)
            for endpt_idx2 = 1:numel(endpt_indices)
                % check if the points are connected
                pt1 = [endpt_x(endpt_idx1), endpt_y(endpt_idx1)];
                pt2 = [endpt_x(endpt_idx2), endpt_y(endpt_idx2)];
                is_connected = check_connection(I, pt1, pt2);
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
                    is_connected = check_connection(bresenhamified, pt1, pt2); % TODO: consider making a wrapper function for checking connectivity of a set of points
                    if is_connected
                        distances(endpt_idx1, endpt_idx2) = 0;
                    end
                end
            end
        end
    end

    
    % now there should only be two endpoints
    % re-skeletonize (in case some corners were introduced and endpoints
    % can't be found)
    bresenhamified = bwskel(logical(bresenhamified));
    surrounding_pixels = conv2(bresenhamified, kernel, 'same');
    endpoints = (surrounding_pixels == 1 & bresenhamified == 1);
    endpt_indices = find(endpoints);
    [endpt_x, endpt_y] = ind2sub(size(endpoints), endpt_indices);
    if numel(endpt_x) ~= 2
        warning('something looks awefully suspicious. double check this layer')
        [endpt_x, endpt_y] = remove_erroneous_endpoints(bresenhamified, endpt_x, endpt_y);
    end
    bresenhamified = draw_line(bresenhamified, endpt_x(1), endpt_y(1), endpt_x(2), endpt_y(2));
end


% ================
% == FLOOD-FILL ==
% ================
function filled = fill_outline(outline)
    inverse_fill = imfill(outline, [1,1]);
    filled = ~inverse_fill + outline;
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
    
    rv_dist = bwdist(new_mask == 2);
    
    % for each voxel, see if a) it's epi, b) if it's closer
    for row = 1:size(new_mask, 1)
        for col = 1:size(new_mask, 2)
            for slice = 1:size(new_mask, 3)
                if new_mask(row,col,slice) == 3 && rv_dist(row,col,slice) < lv_dist(row,col,slice)
                    new_mask(row,col,slice) = 4;
                end
            end
        end
    end
end


function [filled, closed, traced] = create_one_slice(zlevel, slice_size, pts, connections, pt_idxs_to_use, face_idxs_to_use)
    traced = trace_surface(zlevel, slice_size, pts, connections, pt_idxs_to_use, face_idxs_to_use);
    closed = close_big_holes(traced);
    filled = fill_outline(closed);
end


%%
% ==========
% == MAIN ==
% ==========
clear; close all;

% read atlas and connectivity, declare ranges
[pts, ~, ~, ~] = read_atlas("ED", false);
connections = read_connectivity();

LV_pts = 1:1500;
RV_pts = 1501:3224;
EPI_pts = 3225:5582;
LV_faces = 1:3072;
RV_faces = cat(2, 3073:4480, 4480:6752);
EPI_faces = 6753:11616;

% parameters for slices
z_values = -40:2:50;
slice_size = 150;

% creating masks
layered_masks = zeros(slice_size, slice_size, numel(z_values));
traced_masks = zeros(slice_size, slice_size, numel(z_values));
closed_masks = zeros(slice_size, slice_size, numel(z_values));

for z_idx = 1:numel(z_values)
    z = z_values(z_idx);
    % generate each slice
    [LV_slice, LV_closed, LV_traced] = create_one_slice(z, slice_size, pts, connections, LV_pts, LV_faces);
    [RV_slice, RV_closed, RV_traced] = create_one_slice(z, slice_size, pts, connections, RV_pts, RV_faces);
    [EPI_slice, EPI_closed, EPI_traced] = create_one_slice(z, slice_size, pts, connections, EPI_pts, EPI_faces);

    % layer everything
    layered_slice = EPI_slice .* 3;
    layered_slice(LV_slice > 0) = 1;
    layered_slice(RV_slice > 0) = 2;
    layered_masks(:,:,z_idx) = layered_slice;

    % also do the traced masks (for debugging purposes)
    closed_masks(:,:,z_idx) = EPI_closed;
    traced_masks(:,:,z_idx) = EPI_traced;
    fprintf("z slice at %i (slice #%i) has been made\n", z, z_idx);
end

layered_masks = separate_myocardium(layered_masks);

% separate out masks
masks = zeros(slice_size, slice_size, numel(z_values), 4);
masks(:,:,:,1) = layered_masks == 1;
masks(:,:,:,2) = layered_masks == 2;
masks(:,:,:,3) = layered_masks == 3;
masks(:,:,:,4) = layered_masks == 4;

%%
for slice = 1:46
    figure; imshow(layered_masks(:,:,slice), [0, max(layered_masks, [], 'all')]);
end

%%
matlab.io.saveVariablesToScript('outputs/cardiac_mask_10.m', {'masks'})
