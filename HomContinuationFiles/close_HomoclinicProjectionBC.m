function prob = close_HomoclinicProjectionBC(prob, data)
    % hom_close_proj encodes projection boundary conditions for homoclinic continuation.
    %
    % Inputs:
    %   - prob: COCO continuation problem structure.
    %   - data: Data structure containing initialization data.
    %
    % Output:
    %   - prob: Updated COCO continuation problem structure with added boundary conditions.

    % Extract 'coll' toolbox data for the homoclinic continuation
    [data_h, uidx_h] = coco_get_func_data(prob, 'x.coll', 'data', 'uidx');
    
    % Extract 'ep' toolbox data (equilibrium point)
    [data_0, uidx_0] = coco_get_func_data(prob, 'x_ss.ep', 'data', 'uidx');

    % Add state- and parameter-space dimensions to input data structure
    data.xdim = data_h.xdim;  % Dimension of the state space
    data.pdim = data_h.pdim;  % Dimension of the parameter space

    % Store mapping indices from 'coll' and 'ep' data
    maps_h  = data_h.coll_seg.maps;  % Collocation segment mappings
    maps_0 = data_0.ep_eqn;          % Equilibrium equation mappings

    % Initialize eigendata as empty
    data.eigdata = [];

    % Glue equilibrium ('ep') and homoclinic ('coll') problem parameters
    prob = coco_add_glue(prob, 'shared', ...
                         uidx_h(maps_h.p_idx), ...  % Parameter indices from 'coll'
                         uidx_0(maps_0.p_idx));    % Parameter indices from 'ep'

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%    Standard Homoclinic Boundary Conditions   %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Add a slot to update eigenvalues before each continuation step
    data.eq_tbid = coco_get_id('x_ss', 'ep');  % Get identifier for 'ep' equilibrium
    data = coco_func_data(data);               % Prepare function data
    prob = coco_add_slot(prob, data.eq_tbid, @eig_update, data, 'update');  % Add eigenvalue update slot

    % Add a slot to update the phase condition before each continuation step
    data.coll_tbid = coco_get_id('x', 'coll');  % Get identifier for 'coll' continuation
    data = coco_func_data(data);                % Prepare function data
    prob = coco_add_slot(prob, data.coll_tbid, @phase_update, data, 'update');  % Add phase condition update slot

    % Apply projection boundary conditions
    prob = coco_add_func(prob, 'bcs_proj', @BC_projConditions, data, 'zero', ...
                        'uidx', [uidx_h(maps_h.x1_idx); ...    % Indices for x1
                                 uidx_h(maps_h.x0_idx); ...    % Indices for x0
                                 uidx_0(maps_0.x_idx); ...     % Indices for equilibrium
                                 uidx_h(maps_h.p_idx)]);        % Indices for parameters
    prob = coco_add_slot(prob, 'bcs_proj', @coco_save_data, data, 'save_full');  % Add data saving slot

    % Apply phase condition
    prob = coco_add_func(prob, 'bcs_phase', @BC_phase, data, 'zero', ...
                        'uidx', uidx_h(maps_h.xbp_idx), ...          % Indices for phase condition
                        'remesh', @coco_remesh);                     % Remeshing function
    prob = coco_add_slot(prob, 'bcs_phase', @coco_save_data, data, 'save_full');  % Add data saving slot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%       Monitor Homoclinic Degeneracies        %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Define special points to detect during continuation
    SP_points = {'NSS', 'NSF', 'DRS', 'DRU', 'TLS', 'TLR', ...
                'NDS', 'NDU', 'OFS', 'OFU', 'IFS', 'IFU', 'RES', 'EqType'};

    % Add a homoclinic monitoring test function to detect special points
    prob = coco_add_func(prob, 'homTst', @homoclinicTestFunction, data, ...
                        'regular', SP_points, ...
                        'uidx', [uidx_h(maps_h.x1_idx); ...     % Indices for x1
                                 uidx_h(maps_h.x0_idx); ...     % Indices for x0
                                 uidx_0(maps_0.x_idx); ...      % Indices for equilibrium
                                 uidx_h(maps_h.p_idx); ...      % Indices for parameters
                                 uidx_h(maps_h.T_idx)]);         % Indices for period T

    prob = coco_add_slot(prob, 'homTst', @coco_save_data, [], 'save_full');  % Add data saving slot

    % Add events for each special point to be detected
    for i = 1:length(SP_points)
        event_name = strcat(SP_points{i}, '_');  % Create event name by appending "_"
        prob = coco_add_event(prob, event_name, SP_points{i}, 0);  % Add event to COCO
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%    Define Problem Parameters    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Add system parameters to the continuation problem
    prob = coco_add_pars(prob, 'pars_system', uidx_h(maps_h.p_idx), ...
                         {'mu', 'eta', 'k2', 'ep'});  % Parameter names

    % Add the half-period parameter (inactive by default)
    prob = coco_add_pars(prob, 'pars_T', uidx_h(maps_h.T_idx), 'T_hom', 'inactive');  % Half-period parameter
end

function data = eig_update(prob, data, cseg, varargin)
    % eig_update updates eigenvectors of the adjoint problem before each continuation step.
    %
    % Inputs:
    %   - prob: COCO continuation problem structure.
    %   - data: Data structure containing necessary information.
    %   - cseg: Current continuation segment data.
    %   - varargin: Additional arguments (unused).
    %
    % Output:
    %   - data: Updated data structure with new eigendata.

    % Retrieve equilibrium point data and indices
    [data_0, uidx] = coco_get_func_data(prob, data.eq_tbid, 'data', 'uidx');
    maps_0 = data_0.ep_eqn;                % Mappings for equilibrium equations
    u = cseg.src_chart.x(uidx);            % Extract relevant variables from continuation segment
    p0 = u(maps_0.p_idx);                  % Extract current parameters
    x0 = u(maps_0.x_idx);                  % Extract current equilibrium point
    
    data.eigdata = computeOrthogComplement(x0, p0, data);  % Recompute eigenvalue data based on current state and parameters
    data.p0 = p0;                           % Update parameters in data structure
end

function data_out = computeOrthogComplement(x0, p0, data)
    % computeOrthogComplement recomputes eigenspace information based on the current state and parameters.
    %
    % Inputs:
    %   - x0: Current equilibrium point.
    %   - p0: Current parameters.
    %
    % Output:
    %   - data_out: Cell array containing stable and unstable orthog. complements eigenvectors.

    f = data.f;  % Access global function structure (assumed to contain system definitions)
    J = f.dfdx(x0, p0);  % Compute Jacobian matrix at the equilibrium point

    % Calculate eigenvalues and eigenvectors of the Jacobian
    [eigvec, eigval] = eig(J);

    % Sort eigenvalues and corresponding eigenvectors
    [eigval,ind] = sort(real(diag(eigval)));
    eigvec = eigvec(:, ind);             % Rearrange eigenvectors accordingly

    % Identify unstable and stable eigenvectors based on eigenvalues
    unstable_index = eigval > 0;  % Indices of unstable eigenvalues (positive real part)
    stable_index = eigval < 0;    % Indices of stable eigenvalues (negative real part)

    % Extract and normalize unstable and stable eigenvectors
    v_un = eigvec(:, unstable_index);
    v_st = eigvec(:, stable_index);
    v_un = v_un ./ vecnorm(v_un, 2, 1);  % Normalize columns
    v_st = v_st ./ vecnorm(v_st, 2, 1);  % Normalize columns

    % Orthog. complement of eignspaces
    W_UN_star = null(v_un')';  % Unstable orthog. complement
    V_ST_star = null(v_st')';  % Stable orthog. complement

    % Store eigendata for output
    data_out = {V_ST_star, W_UN_star};
end

function data = phase_update(prob, data, cseg, varargin)
    % phase_update updates the phase condition before each continuation step.
    %
    % Inputs:
    %   - prob: COCO continuation problem structure.
    %   - data: Data structure containing necessary information.
    %   - cseg: Current continuation segment data.
    %   - varargin: Additional arguments (unused).
    %
    % Output:
    %   - data: Updated data structure with new phase condition information.

    % Retrieve collocation segment data and indices
    [fdata, uidx] = coco_get_func_data(prob, data.coll_tbid, 'data', 'uidx');
    maps = fdata.coll_seg.maps;  % Mappings for collocation segments
    u = cseg.src_chart.x;        % Extract relevant variables from continuation segment

    % Compute the initial phase condition value
    data.xp0 = u(uidx(maps.xbp_idx))' * data.intfac;
end
