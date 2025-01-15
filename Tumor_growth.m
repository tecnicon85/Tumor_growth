     
clear ; close all;
% Read input data from CSV files
cells = csvread('body_elems.csv');    % Cell connectivity information
coordinates = csvread('body_point.csv');  % Coordinates of cell vertices
neighbors = csvread('body_esuel.csv'); % Neighbors of each cell

% Compute area for each cell
numCells = length(cells);
areaCells = zeros(numCells, 1); % Initialize array to store cell areas

for i = 1:numCells
    % Extract vertex coordinates for the current cell
    vertices = coordinates(cells(i,:), :);
    
    % Calculate two vectors forming the triangle
    vector1 = vertices(2,:) - vertices(1,:);
    vector2 = vertices(3,:) - vertices(1,:);
    
    % Compute cross product to get area vector
    crossproduct = cross(vector1, vector2);
    
    % Calculate and store the area of the triangle
    areaCells(i) = 0.5 * norm(crossproduct);
end

% Define probability constants for state transitions
p1 = 0.1;  % Probability for normal cell to become cancerous
p2 = 0.1;  % Probability for cancer cell to become complex
p3 = 0.1;  % Probability for complex cell to become necrotic

% Initialize cell states
% 0: normal, 1: cancer, 2: complex, 3: necrotic
cellStatus = zeros(numCells, 1);
cancerCell = 10; % Index of the first cancer cell
cellStatus(cancerCell) = 1;  % Set the first cancer cell

numIterations = 100;  % Number of time steps to simulate
rng(1); % Set random number generator seed for reproducibility

figure; % Create a new figure for visualization

for t = 1:numIterations
    newCellStatus = cellStatus; % Create a copy to update cell statuses 

    for i = 1:numCells
        neighborIndices = neighbors(i, :); % Get indices of neighboring cells
        
        % Rule 1: Normal cell becomes cancer cell
        if cellStatus(i) == 0
            cancerousNeighbors = cellStatus(neighborIndices) == 1; % Check for cancerous neighbors
            if any(cancerousNeighbors)
                areaCancerousNeighbors = sum(areaCells(neighborIndices(cancerousNeighbors))); % Total area of cancerous neighbors
                totalNeighborArea = sum(areaCells(neighborIndices)); % Total area of all neighbors
                p = p1 * areaCancerousNeighbors / (areaCells(i) + totalNeighborArea); % Calculate transition probability
                
                if rand() < p
                    newCellStatus(i) = 1; % Change to cancer state
                end
            end
        end

        % Rule 2: Cancer cell becomes complex cell
        if cellStatus(i) == 1
            neighborCellState = cellStatus(neighbors(i, :));
            if all(neighborCellState ~= 0) && rand() < p2  % If all neighbors are non-normal and probability condition is met
                newCellStatus(i) = 2;  % Change to complex state
            end
        end

        % Rule 3: Complex cell becomes necrotic cell
        if cellStatus(i) == 2 && rand() < p3
            newCellStatus(i) = 3;  % Change to necrotic state
        end
    end

    cellStatus = newCellStatus; % Update cell statuses at the end of each iteration

    % Visualization
    clf; % Clear the current figure
    clr = zeros(length(cells), 3); % Initialize array for element colors
    
    % Assign colors based on cell states
    for i = 1:length(cells)
        cellState = cellStatus(cells(i,1)); 
        if cellState == 0
            clr(i,:) = [0, 1, 0]; % Green for normal cells
        elseif cellState == 1
            clr(i,:) = [1, 0, 0]; % Red for cancer cells
        elseif cellState == 2
            clr(i,:) = [0, 0, 1]; % Blue for complex cells
        elseif cellState == 3
            clr(i,:) = [0.5, 0.5, 0.5]; % Gray for necrotic cells
        end
    end

    % Create a 3D surface plot of the cells
    trisurf(cells, coordinates(:,1), coordinates(:,2), coordinates(:,3), 'FaceVertexCData', clr, 'FaceColor', 'flat');
    daspect([1 1 1]); % Set aspect ratio to 1:1:1
    colormap('jet'); % Set colormap for visualization
    colorbar; % Add colorbar for reference
    drawnow; % Update the display

    % Print statistics for the current time step
    fprintf('Time Step %d: Normal=%d, Cancer=%d, Complex=%d, Necrotic=%d\n', ...
             t, sum(cellStatus == 0), sum(cellStatus == 1), ...
             sum(cellStatus == 2), sum(cellStatus == 3));
end
