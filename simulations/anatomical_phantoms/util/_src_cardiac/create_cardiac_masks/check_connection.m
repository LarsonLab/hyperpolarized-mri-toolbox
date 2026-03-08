function isconnected = check_connection(BW, pt1, pt2)
    % checks if 2 points are connected by 1s on a black and white image.
    % basically just a breadth-first search pathfinding algorithm that
    % returns true if a path exists and false if path doesn't exist.
    %
    % ARGS:
    %   BW = binary image. size [N x M]
    %   pt1 = coordinates of the first point. size [1 2]
    %   pt2 = coordinates of the second point. size [1 2]
    
    % if the 2 points are the same, just say they're connected
    if pt1 == pt2
        isconnected = true;
        return
    end

    isconnected = false;
    directions = [0 1; 0 -1; 1 0; -1 0; ... % 4 cardinal directions
              1 1; 1 -1; -1 1; -1 -1]; % 4 corners

    queue = {pt1};
    
    visited = zeros(size(BW));
    visited(pt1(1), pt1(2)) = 1; % mark initial point as visited
    
    while ~isempty(queue)
        current_point = queue{1};
        % check each direction
        for direction_n = 1:size(directions, 1)
            neighbor_idx = current_point + directions(direction_n, :);
            % ensure valid idx
            if any(neighbor_idx < 1) || any(neighbor_idx > size(BW))
                continue
            end
            % add to queue if unsearched point
            if BW(neighbor_idx(1), neighbor_idx(2)) == 1 & ...
                    visited(neighbor_idx(1), neighbor_idx(2)) == 0
                queue{end + 1} = neighbor_idx;
                visited(neighbor_idx(1), neighbor_idx(2)) = 1;
                if neighbor_idx == pt2
                    isconnected = true;
                    return
                end
            end
        end
        % pop point of queue off the queue
        queue = queue(2:end);
    end
end
