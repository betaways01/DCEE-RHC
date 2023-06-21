function point = getPointByIndex(currentPosition, index)
    % Define the 8 directions in which the agent can move
    directions = [0 1; 1 1; 1 0; 1 -1; 0 -1; -1 -1; -1 0; -1 1];
    
    % Get the direction corresponding to the index
    direction = directions(index, :);
    
    % Calculate the coordinates of the point
    point = currentPosition + direction;
end
