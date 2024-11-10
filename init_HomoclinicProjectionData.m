function data_out = init_HomoclinicProjectionData(poFn, hom_idx)
    % hom_proj_init_data initializes data for projection boundary condition homoclinic continuation.
    %
    % ARGS:
    %   - poFn: Filename of the periodic orbit (PO) continuation run.
    %   - hom_idx: Index of the candidate homoclinic from the PO continuation run.
    %
    % OUTPUT:
    %   data_out: Structure containing initialized data for homoclinic continuation.

    % If hom_idx is not provided, use the last label in the PO continuation run as the approximate homoclinic
    if nargin < 2
        bd = coco_bd_read(poFn);           % Read the continuation database
        col = coco_bd_col(bd);             % Get column data
        hom_idx = col{end};                % Use the last column as homoclinic index
    end

    % Attempt to read the PO data; fallback to 'homo.po.orb' if 'po.orb' fails
    try
        [hom_po, ~] = coll_read_solution('po.orb', poFn, hom_idx);  % Read PO solution
    catch
        fprintf("Trying homo.po.orb \n")
        [hom_po, ~] = coll_read_solution('homo.po.orb', poFn, hom_idx);  % Read alternative PO solution
    end

    % Evaluate the function to find the equilibrium point closest to the PO
    fEval = AutoFuncODE(hom_po.xbp, repmat(hom_po.p, [1 size(hom_po.xbp, 1)])); 
    [~, idx] = min(norm(fEval));        % Find index with minimum norm
    x_ss = hom_po.xbp(idx, :);         % Equilibrium point

    % Initial scaling factors along each linear eigenspace
    eps_u_guess = 1e-4;
    eps_s_guess = 1e-4;
    [x_u, x_s] = locatePointsAtArcLengths(hom_po.xbp, x_ss, eps_s_guess, eps_u_guess);   
    eps_u = norm(x_ss - x_u);
    eps_s = norm(x_ss - x_s);
    epsilon0 = [eps_s, eps_u];         % Epsilon array to return

    % Remove parts of the homoclinic orbit outside the vicinity of the equilibrium
    whereXs = hom_po.xbp == x_s;
    row_indices_s = find(all(whereXs == 1, 2));      % Indices matching x_s
    whereXu = hom_po.xbp == x_u;
    row_indices_u = find(all(whereXu == 1, 2));      % Indices matching x_u

    hom_po.xbp(1:row_indices_u-1, :) = NaN;         % Remove before x_u
    hom_po.tbp(1:row_indices_u-1, :) = NaN;
    hom_po.xbp(row_indices_s+1:end, :) = NaN;       % Remove after x_s
    hom_po.tbp(row_indices_s+1:end, :) = NaN;

    hom_po.tbp(isnan(hom_po.tbp), :) = [];          % Remove NaNs from time
    hom_po.xbp(isnan(hom_po.xbp(:,1)), :) = [];     % Remove NaNs from state

    % Compute the new period of the truncated solution
    T = hom_po.tbp(end) - hom_po.tbp(1);
    hom_po.tbp = hom_po.tbp - hom_po.tbp(1);        % Normalize time

    % Plot the phase space trajectory
    figure
    hold on
        plot3(hom_po.xbp(:,1), hom_po.xbp(:,2), hom_po.xbp(:,3), 'linewidth', 1, 'color', 'blue')  % Trajectory
        scatter3(x_ss(1), x_ss(2), x_ss(3), 5, 'MarkerFaceColor', 'black')  % Equilibrium point
        scatter3(x_u(1), x_u(2), x_u(3), 5, 'MarkerFaceColor', 'red')      % Unstable point
        scatter3(x_s(1), x_s(2), x_s(3), 5, 'MarkerFaceColor', 'blue')     % Stable point
        draw_line(x_ss, x_u, 'red')  % Line from equilibrium to unstable
        draw_line(x_ss, x_s, 'blue') % Line from equilibrium to stable
        view(3)                      % 3D view
        xlim([x_ss(1)-0.1, x_ss(1)+0.1])  % Set x-axis limits
        ylim([x_ss(2)-0.1, x_ss(2)+0.1])  % Set y-axis limits
        zlim([x_ss(3)-0.1, x_ss(3)+0.1])  % Set z-axis limits
    hold off

    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%     Output     
    %%%%%%%%%%%%%%%%%%%%%%%%
    data_out.t_sol = hom_po.tbp;              % Time solution
    data_out.x_sol = hom_po.xbp;              % State solution
    data_out.p0    = hom_po.p;                 % Parameters
    data_out.pdim  = length(hom_po.p);         % Number of parameters
    data_out.xdim  = length(hom_po.xbp);       % Dimension of state space
    data_out.x_ss  = x_ss;                     % Equilibrium point
    data_out.T     = T;                         % Period
    data_out.pnames= {'mu', 'eta', 'k2', 'ep'};% Parameter names
    data_out.t0    = 0;                         % Initial time
    data_out.x_init_u = x_u';                  % Initial unstable point
    data_out.x_init_s = x_s';                  % Initial stable point
    data_out.eps = epsilon0;                    % Epsilon scaling factors

    %%% TESTING PLOTS (Commented Out)
    % figure
    % hold on
    %     plot(hom_po.tbp(:,1), hom_po.xbp(:,2))
    % hold off
end


function [x1, x2] = locatePointsAtArcLengths(C, x0, eps_s, eps_u)
    % find_points_at_arc_length finds points on a D-dimensional curve at specified arc lengths from x0.
    %
    % Inputs:
    %   - C: NxD array of curve points, ordered sequentially.
    %   - x0: 1xD array representing a point on the curve C.
    %   - eps_s: Arc length distance to the left of x0.
    %   - eps_u: Arc length distance to the right of x0.
    %
    % Outputs:
    %   - x1: Point on the curve at arc length eps_u to the right of x0.
    %   - x2: Point on the curve at arc length eps_s to the left of x0.

    % Validate inputs
    if size(C,2) < 2
        error('Curve C must have at least two dimensions (N x D, D >= 2).');
    end

    if ~isequal(size(x0), [1, size(C,2)])
        error('Point x0 must be a 1 x D array, matching the dimensions of C.');
    end

    % Find the index of x0 in C using a tolerance
    idx = find(ismembertol(C, x0, 1e-6, 'ByRows', true), 1);
    if isempty(idx)
        error('Point x0 is not found on the curve C.');
    end

    % Compute differences and cumulative arc length
    deltas = diff(C, 1, 1);                      % (N-1)xD differences between consecutive points
    segment_lengths = sqrt(sum(deltas.^2, 2));    % Euclidean distance for each segment
    cumulative_length = [0; cumsum(segment_lengths)]; % Cumulative arc length

    % Locate points at specified arc lengths using helper function
    x1 = locate_point(C, cumulative_length, idx, 1, eps_u);  % Right of x0
    C = flip(C);                                                     % Reverse curve for left search
    x2 = locate_point(C, cumulative_length, idx, 1, eps_s);  % Left of x0
end


function x = locate_point(C, cum_length, idx, direction, eps_dist)
    % locate_point_corrected locates a point at a given arc length in a specified direction.
    %
    % Inputs:
    %   - C: NxD array of curve points.
    %   - cum_length: Cumulative arc length array.
    %   - idx: Current index on the curve.
    %   - direction: 1 for forward, -1 for backward.
    %   - eps_dist: Desired arc length distance.
    %
    % Output:
    %   - x: Located point on the curve.

    % Calculate target length based on direction
    if direction == 1
        target_length = cum_length(idx) + eps_dist;
        target_length = min(target_length, cum_length(end) - 1e-10);  % Clamp within bounds
        seg_idx = find(cum_length > target_length, 1);              % Find segment index
    else
        target_length = cum_length(idx) - eps_dist;
        target_length = max(target_length, cum_length(1) + 1e-10);  % Clamp within bounds
        seg_idx = find(cum_length < target_length, 1, 'last');     % Find segment index
    end

    % Handle cases where no valid segment is found
    if isempty(seg_idx)
        if direction == 1
            x = C(end-1, :) + 1e-10 * (C(end, :) - C(end-1, :));  % Near end
        else
            x = C(2, :) - 1e-10 * (C(2, :) - C(1, :));            % Near start
        end
        return;
    end

    % Determine indices for interpolation
    if direction == 1
        i = seg_idx - 1;
    else
        i = seg_idx;
    end

    % Clamp index to valid range
    i = max(1, min(i, size(C, 1) - 1));

    % Perform linear interpolation to find the exact position
    l1 = cum_length(i);
    l2 = cum_length(i + 1);
    ratio = (target_length - l1) / (l2 - l1);
    ratio = max(0, min(1, ratio));  % Ensure ratio is within [0, 1]

    x_interpolated = C(i, :) + ratio * (C(i + 1, :) - C(i, :));  % Interpolated point

    % Find the closest actual point in C to the interpolated point
    distances = vecnorm(C - x_interpolated, 2, 2);  % Euclidean distances
    [~, closest_idx] = min(distances);             % Index of closest point
    x = C(closest_idx, :);                        % Assign closest point
end
