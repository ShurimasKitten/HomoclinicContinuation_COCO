function hom_data = hom2hom_ProjData(poFn, hom_idx, continuationSettings)
    % init_homProjData Initializes data for projection boundary condition homoclinic continuation
    %    
    % ARGS:
    %   - poFn (string): Filename of the PO continuation run.
    %   - hom_idx (integer, optional): Index of the candidate homoclinic orbit from the PO continuation run.
    %                                  If not provided, the last entry in
    %                                  the PO continuation is used.
    %
    %   - hom_data: Hom data used in previous continuation run.
    %
    % RETURNS:
    %   - data_out (struct): A structure containing initialized data for homoclinic continuation.
    %
    % Example:
    %   data = hom_expl_init_data('po_continuation.mat', 10, [1e-3, 1e-3]);

    % Try reading the solution from 'po.orb' file
    [hom_po, ~] = coll_read_solution('hom.orb', poFn, hom_idx);
    [x_ss, ~] = ep_read_solution('hom', poFn, hom_idx);
    x_ss = real(x_ss.x);

    % Compute the new period T of the truncated homoclinic solution
    T = hom_po.tbp(end) - hom_po.tbp(1);
    hom_po.tbp = hom_po.tbp - hom_po.tbp(1);  % Normalize time to start at 0

    %%% Plot the phase space trajectory for visual inspection
    figure
    hold on
        plot3(hom_po.xbp(:,1), hom_po.xbp(:,2), hom_po.xbp(:,3), 'linewidth', 0.5, 'color', 'black');
        scatter3(x_ss(1), x_ss(2), x_ss(3), 5, 'MarkerFaceColor', 'black'); 
        view(3);                        
    hold off

    %%%     Output     
    hom_data.t_sol = hom_po.tbp;                 % Time solution
    hom_data.x_sol = hom_po.xbp;                 % State variables solution
    hom_data.p0    = hom_po.p;                   % Parameter values
    hom_data.pdim  = length(hom_po.p);           % Number of parameters
    hom_data.xdim  = length(hom_po.xbp);         % Dimension of state space
    hom_data.x_ss    = x_ss;                     % Equilibrium point
    hom_data.T     = T;                          % Period of the homoclinic orbit
    hom_data.pnames= {'mu', 'eta', 'k2', 'ep'};  % Names of the parameters
    hom_data.t0    = 0;                          % Initial time
    hom_data.f = continuationSettings.f;         % COCO vec.f
    hom_data.counter = 0;                        % Counter used in homoclinic error handling
    hom_data.hom_cid = coco_get_id('hom', 'orb');
    hom_data.ep_cid = coco_get_id('hom', '');

end
