
function draw_line(start_point, end_point, line_color)
    % DRAW_LINE_3D Draws a line from a start point to an end point using the first
    % three components of points in R^4, with the specified line color.
    %
    % Inputs:
    %   start_point - A 1x4 vector [x1, y1, z1, w1] representing the coordinates of the start point in R^4.
    %   end_point   - A 1x4 vector [x2, y2, z2, w2] representing the coordinates of the end point in R^4.
    %   line_color  - A character or string specifying the color of the line.
    %
    % Example:
    %   draw_line_3d([1, 2, 3, 4], [4, 5, 6, 7], 'b');


    if nargin < 3
        line_color = 'black';
    end

    % Extract the first three coordinates from input vectors
    x1 = start_point(1);
    y1 = start_point(2);
    z1 = start_point(3);
    
    x2 = end_point(1);
    y2 = end_point(2);
    z2 = end_point(3);

    % Create a figure and plot the line in 3D space
    plot3([x1, x2], [y1, y2], [z1, z2], 'Color', line_color, 'LineWidth', 1.5); % Draw the line with specified color

    % Adding labels and title for better visualization
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    title('Line from Start Point to End Point in R^4 (First Three Components)');
    grid on; % Turn on the grid for better visual reference

    % Highlight the start and end points
    plot3(x1, y1, z1, 'ro', 'MarkerSize', 5, 'DisplayName', 'Start Point', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black'); % Mark the start point in red
    plot3(x2, y2, z2, 'go', 'MarkerSize', 5, 'DisplayName', 'End Point', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black'); % Mark the end point in green
    % legend('Line', 'Start Point', 'End Point');

end
