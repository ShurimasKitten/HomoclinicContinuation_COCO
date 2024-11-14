function data_out = init_homProjData(poFn, hom_idx, continuationSettings)
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

    % Try reading the solution from 'po.orb' file
    [hom_po, ~] = coll_read_solution('x', poFn, hom_idx);
    [x_ss, ~] = ep_read_solution('x_ss', poFn, hom_idx);
    x_ss = x_ss.x;

    % Compute the new period T of the truncated homoclinic solution
    T = hom_po.tbp(end) - hom_po.tbp(1);
    % hom_po.tbp = hom_po.tbp - hom_po.tbp(1);  % Normalize time to start at 0


    %% Plot the phase space trajectory for visual inspection
    figure;
    hold on;
        % Plot the 3D trajectory of the homoclinic orbit
        plot3(hom_po.xbp(:,1), hom_po.xbp(:,2), hom_po.xbp(:,3), 'linewidth', 0.5, 'color', 'black');

        % Scatter plot of the equilibrium point and the initial unstable and stable points
        scatter3(x_ss(1), x_ss(2), x_ss(3), 5, 'MarkerFaceColor', 'black');   % Equilibrium

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
    data_out.f = continuationSettings.f;         % COCO vec.f
    data_out.counter = 0;                        % Counter used in homoclinic error handling
end
