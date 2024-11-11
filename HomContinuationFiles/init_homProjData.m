function data_out = init_homProjData(poFn, hom_idx, eps)
    % init_homProjData Initializes data for projection boundary condition homoclinic continuation
    %
    % This function sets up the necessary data for performing homoclinic continuation
    % using projection boundary conditions. It processes the periodic orbit (PO) continuation
    % results to identify and initialize a homoclinic orbit.
    %
    % ARGS:
    %   - poFn (string): Filename of the PO continuation run.
    %   - hom_idx (integer, optional): Index of the candidate homoclinic orbit from the PO continuation run.
    %                                    If not provided, the last entry in the PO continuation is used.
    %   - eps (1x2 array, optional): Initial guesses for the distances along the stable and unstable manifolds.
    %                                  If not provided, default values are used.
    %
    % RETURNS:
    %   - data_out (struct): A structure containing initialized data for homoclinic continuation.
    %
    % Example:
    %   data = hom_expl_init_data('po_continuation.mat', 10, [1e-3, 1e-3]);

    % Default behavior: If homoclinic index is not provided, use the last entry from the PO continuation
    if nargin < 2
        bd = coco_bd_read(poFn);          % Read the boundary data from the PO continuation file
        col = coco_bd_col(bd);            % Extract column data from the boundary data
        hom_idx = col{end};               % Use the last index as the homoclinic candidate
    end

    % Attempt to read the periodic orbit (PO) data corresponding to hom_idx
    try
        % Try reading the solution from 'po.orb' file
        [hom_po, ~] = coll_read_solution('po.orb', poFn, hom_idx);
    catch
        % If reading 'po.orb' fails, attempt to read from 'x.coll' as a fallback
        fprintf("Trying homo.coll \n");
        [hom_po, ~] = po_read_solution('homo', poFn, hom_idx);
    end

    % Evaluate the function vector to find the equilibrium point closest to the PO
    fEval = AutoFuncODE(hom_po.xbp, repmat(hom_po.p, [1 size(hom_po.xbp, 1)])); 
    [~, idx] = min(norm(fEval));           % Find the index with the minimum norm of fEval
    x_ss = hom_po.xbp(idx, :);            % Extract the equilibrium point

    % Initialize scaling factors along the unstable and stable eigenspaces
    if nargin == 3
        % If eps is provided, use its elements as initial guesses
        eps_u_guess = eps(1);
        eps_s_guess = eps(2);
    else
        % If eps is not provided, use default small values
        eps_u_guess = 1e-7;
        eps_s_guess = 1e-7;
    end
    % Find points on the PO at specified arc lengths from the equilibrium
    [x_u, x_s] = find_points_at_arc_length(hom_po.xbp, x_ss, eps_s_guess, eps_u_guess);   

    % Calculate the actual distances from the equilibrium to the found points
    eps_u = norm(x_ss - x_u);
    eps_s = norm(x_ss - x_s);

    % Build the initial explicit boundary conditions using the equilibrium and parameters
    [v_un, v_st] = buildExplBC(x_ss', hom_po.p);

    % Compute coefficients for linear combination of start/end points in terms of eigenbasis
    coeff_s = compute_coefficients(x_ss', v_st', x_s');
    coeff_u = compute_coefficients(x_ss', v_un', x_u');

    % Scale the epsilon values by the computed coefficients
    eps_u  = eps_u * coeff_u(1:size(v_un,1), :)';
    eps_s  = eps_s * coeff_s(1:size(v_st,1), :)';

    % Combine the scaled epsilon values into an array to return
    epsilon0 = [eps_s, eps_u];

    % Construct initial boundary conditions for visual verification

    % Boundary condition for unstable manifold (Eu)
    X_u = x_ss;
    for i = 1:size(v_un, 1)
        X_u = X_u + v_un(i, :) * eps_u(i);
    end

    % Boundary condition for stable manifold (Es)
    X_s = x_ss;
    for i = 1:size(v_st, 1)
        X_s = X_s + v_st(i, :) * eps_s(i);
    end

    % Remove segments of the homoclinic orbit that lie between x_ss and x_u/x_s
    % Identify indices where x_s and x_u occur in the PO data
    whereXs = hom_po.xbp == x_s;
    row_indices_s = find(all(whereXs == 1, 2));
    whereXu = hom_po.xbp == x_u;
    row_indices_u = find(all(whereXu == 1, 2));

    % Set points before x_u to NaN to exclude them from the homoclinic orbit
    hom_po.xbp(1:row_indices_u-1, :) = NaN;
    hom_po.tbp(1:row_indices_u-1, :) = NaN;

    % Set points after x_s to NaN to exclude them from the homoclinic orbit
    hom_po.xbp(row_indices_s+1:end, :) = NaN;
    hom_po.tbp(row_indices_s+1:end, :) = NaN;

    % Remove all NaN entries from the PO data
    hom_po.tbp(isnan(hom_po.tbp), :) = [];
    hom_po.xbp(isnan(hom_po.xbp(:,2)), :) = [];

    % Compute the new period T of the truncated homoclinic solution
    T = hom_po.tbp(end) - hom_po.tbp(1);
    hom_po.tbp = hom_po.tbp - hom_po.tbp(1);  % Normalize time to start at 0

    % Plot the phase space trajectory for visual inspection
    figure;
    hold on;
        % Plot the 3D trajectory of the homoclinic orbit
        plot3(hom_po.xbp(:,1), hom_po.xbp(:,2), hom_po.xbp(:,3), 'linewidth', 0.5, 'color', 'black');

        % Scatter plot of the equilibrium point and the initial unstable and stable points
        scatter3(x_ss(1), x_ss(2), x_ss(3), 5, 'MarkerFaceColor', 'black');   % Equilibrium
        scatter3(x_u(1), x_u(2), x_u(3), 5, 'MarkerFaceColor', 'red');        % Unstable
        scatter3(x_s(1), x_s(2), x_s(3), 5, 'MarkerFaceColor', 'blue');       % Stable

        % Draw lines representing the unstable and stable manifold directions
        draw_line(x_ss, X_u, 'red');   % Unstable manifold direction
        draw_line(x_ss, X_s, 'blue');  % Stable manifold direction
        view(3);                        % Set 3D view

        % Set axis limits around the equilibrium point for better visualization
    hold off;

    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%     Output     
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Populate the output structure with relevant data
    data_out.t_sol = hom_po.tbp;                % Time solution
    data_out.x_sol = hom_po.xbp;                % State variables solution
    data_out.p0    = hom_po.p;                   % Parameter values
    data_out.pdim  = length(hom_po.p);           % Number of parameters
    data_out.xdim  = length(hom_po.xbp);         % Dimension of state space
    data_out.x_ss    = x_ss;                     % Equilibrium point
    data_out.T     = T;                           % Period of the homoclinic orbit
    data_out.pnames= {'mu', 'eta', 'k2', 'ep'};  % Names of the parameters
    data_out.t0    = 0;                           % Initial time
    data_out.x_init_u = x_u';                    % Initial unstable manifold point
    data_out.x_init_s = x_s';                    % Initial stable manifold point
    data_out.eps = epsilon0;                      % Initial epsilon values for manifolds

    %%% TESTING PLOTS
    % Uncomment the following lines to plot the trajectory in the TS plane for testing
    % figure;
    % hold on;
    %     plot(hom_po.tbp(:,1), hom_po.xbp(:,2));
    % hold off;
end

function [x1, x2] = find_points_at_arc_length(C, x0, eps_s, eps_u)
    % FIND_POINTS_AT_ARC_LENGTH Finds points on a curve at specified arc lengths from a reference point
    %
    % This function locates two points on the D-dimensional curve C that are at arc lengths
    % eps_s (stable direction) and eps_u (unstable direction) from the reference point x0.
    %
    % Inputs:
    %   - C (NxD array): Curve represented as N sequential points in D-dimensional space.
    %   - x0 (1xD array): Reference point on the curve C.
    %   - eps_s (scalar): Arc length distance to find the stable manifold point.
    %   - eps_u (scalar): Arc length distance to find the unstable manifold point.
    %
    % Outputs:
    %   - x1 (1xD array): Point on the curve at arc length eps_u to the right of x0.
    %   - x2 (1xD array): Point on the curve at arc length eps_s to the left of x0.
    %
    % Example:
    %   [xu, xs] = find_points_at_arc_length(curve, [0,0], 1.0, 1.0);

    % Validate that the curve has at least two dimensions
    if size(C,2) < 2
        error('Curve C must have at least two dimensions (N x D, D >= 2).');
    end

    % Validate that the reference point x0 matches the dimensions of C
    if ~isequal(size(x0), [1, size(C,2)])
        error('Point x0 must be a 1 x D array, matching the dimensions of C.');
    end

    % Locate the index of x0 in the curve C using a tolerance to account for numerical precision
    idx = find(ismembertol(C, x0, 1e-6, 'ByRows', true), 1);
    if isempty(idx)
        error('Point x0 is not found on the curve C.');
    end

    % Compute differences between consecutive points to calculate segment lengths
    deltas = diff(C, 1, 1);                    % (N-1) x D array of differences
    segment_lengths = sqrt(sum(deltas.^2, 2));  % Euclidean distance for each segment
    cumulative_length = [0; cumsum(segment_lengths)];  % Cumulative arc length from the start

    % Locate x1 (unstable manifold point) using a helper function
    x1 = locate_point_corrected(C, cumulative_length, idx, 1, eps_u);
    % Locate x2 (stable manifold point) by flipping the curve and using the helper function
    C = flip(C);
    x2 = locate_point_corrected(C, cumulative_length, idx, 1, eps_s);
end

function x = locate_point_corrected(C, cum_length, idx, direction, eps_dist)
    % LOCATE_POINT_CORRECTED Locates a point at a specified arc length from a reference index
    %
    % This helper function finds a point on the curve C that is at an arc length
    % eps_dist from the point at index idx in the specified direction.
    %
    % Inputs:
    %   - C (NxD array): Curve represented as N points in D-dimensional space.
    %   - cum_length (Nx1 array): Cumulative arc lengths corresponding to points in C.
    %   - idx (integer): Index of the reference point in C.
    %   - direction (integer): 1 for forward (right), -1 for backward (left) traversal.
    %   - eps_dist (scalar): Arc length distance to locate the target point.
    %
    % Output:
    %   - x (1xD array): Located point at the specified arc length.
    %
    % Example:
    %   point = locate_point_corrected(curve, cum_len, 10, 1, 1.0);

    % Determine the target cumulative length based on the direction
    if direction == 1
        % Forward traversal: target length is current length plus eps_dist
        target_length = cum_length(idx) + eps_dist;
        
        % Ensure the target length does not exceed the curve's total length
        target_length = min(target_length, cum_length(end) - 1e-10);  % Slightly less to stay inside
        
        % Find the segment index where the target length falls
        seg_idx = find(cum_length > target_length, 1);  % Use '>' to stay strictly inside
    else
        % Backward traversal: target length is current length minus eps_dist
        target_length = cum_length(idx) - eps_dist;
        
        % Ensure the target length does not go below the start of the curve
        target_length = max(target_length, cum_length(1) + 1e-10);  % Slightly greater to stay inside
        
        % Find the segment index where the target length falls
        seg_idx = find(cum_length < target_length, 1, 'last');  % Use '<' to stay strictly inside
    end

    if isempty(seg_idx)
        % If no valid segment is found, return a point very close to the nearest endpoint
        if direction == 1
            % Forward direction: near the end of the curve
            x = C(end-1, :) + 1e-10 * (C(end, :) - C(end-1, :)); 
        else
            % Backward direction: near the start of the curve
            x = C(2, :) - 1e-10 * (C(2, :) - C(1, :)); 
        end
        return;
    end

    % Determine the indices for interpolation within the identified segment
    if direction == 1
        i = seg_idx - 1;
    else
        i = seg_idx;
    end

    % Clamp the index to ensure it is within valid bounds for interpolation
    i = max(1, min(i, size(C, 1) - 1));

    % Retrieve the cumulative lengths at the segment boundaries
    l1 = cum_length(i);
    l2 = cum_length(i + 1);
    % Compute the ratio of the remaining distance within the segment
    ratio = (target_length - l1) / (l2 - l1);
    % Ensure the ratio is within [0, 1] to avoid extrapolation
    ratio = max(0, min(1, ratio)); 

    % Perform linear interpolation to find the exact location of the target point
    x_interpolated = C(i, :) + ratio * (C(i + 1, :) - C(i, :));

    % Find the closest actual point in C to the interpolated point
    distances = vecnorm(C - x_interpolated, 2, 2);  % Euclidean distances to all points
    [~, closest_idx] = min(distances);               % Index of the closest point

    % Assign the closest point as the located point
    x = C(closest_idx, :);
end
