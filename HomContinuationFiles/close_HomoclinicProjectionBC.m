function prob = close_HomoclinicProjectionBC(prob, data)
    % close_HomoclinicProjectionBC: Close the homoclinic continuation problem.
    %   - Applies distributed projection boundary conditions along with a remeshing function.
    %   - Implements projection boundary conditions.
    %
    % Inputs:
    %   - data: Contains complete information about the phase condition and
    %           the initial homoclinic solution to initiate continuation.
    %
    % Potential Issues:
    %   - Unlikely, but is `x1_idx` and `x0_idx` preserved under remeshing?
    %     If not, this could lead to significant issues.
    %

    %%% Set up 
    % Extract toolbox data 
    [data_h, uidx_h] = coco_get_func_data(prob, coco_get_id(data.hom_cid, 'coll'), 'data', 'uidx');    
    [data_0, uidx_0] = coco_get_func_data(prob, coco_get_id(data.ep_cid, 'ep'), 'data', 'uidx');

    % Dims
    data.xdim = data_h.xdim;
    data.pdim = data_h.pdim;  

    % Get collocation segment data f
    maps_h = data_h.coll_seg.maps; 
    maps_0 = data_0.ep_eqn;        

    % Initialize eigendata as empty
    data.eigdata = [];

    % Share problem parameters between 'coll' and 'ep' toolboxes
    prob = coco_add_glue(prob, 'shared', uidx_h(maps_h.p_idx), uidx_0(maps_0.p_idx));    

    %%% Standard Homoclinic Boundary Conditions  
    data = coco_func_data(data);  % Prep. COCO data.              
    prob = coco_add_slot(prob, coco_get_id('homo.orb', 'coll'), @eig_update, data, 'update');  % Add eigenvalue update slot

    % Add a slot to update the phase condition before each continuation step
    prob = coco_add_slot(prob, coco_get_id('hom', 'ep'), @phase_update, data, 'update');  

    % Apply projection boundary conditions
    prob = coco_add_func(prob, 'bcs_proj', @BC_projConditions, data, 'zero', ...
                        'uidx', [uidx_h(maps_h.x1_idx); ...    
                                 uidx_h(maps_h.x0_idx); ...    
                                 uidx_0(maps_0.x_idx); ...     
                                 uidx_h(maps_h.p_idx)]);       
    prob = coco_add_slot(prob, 'bcs_proj', @coco_save_data, data, 'save_full');  

    % Apply phase condition
    prob = coco_add_func(prob, 'bcs_phase', @BC_phase, data, 'zero', ...
                        'uidx', uidx_h(maps_h.xbp_idx), ...          % Indices for phase condition
                        'remesh', @hom_remesh);                      % Remeshing function
    prob = coco_add_slot(prob, 'bcs_phase', @coco_save_data, data, 'save_full'); 

    
    %%% Monitor Homoclinic Degeneracies           
    % Define special points to detect during continuation
    SP_points = {'NSS', 'NSF', 'DRS', 'DRU', 'TLS', 'TLR', ...
                'NDS', 'NDU', 'OFS', 'OFU', 'IFS', 'IFU', 'RES', 'EqType', 'END'};
            
    % Add a homoclinic monitoring test function to detect special points
    prob = coco_add_func(prob, 'homTst', @homoclinicTestFunction, data, ...
                        'regular', SP_points, ...
                        'uidx', [uidx_h(maps_h.x1_idx); ...    
                                 uidx_h(maps_h.x0_idx); ...    
                                 uidx_0(maps_0.x_idx); ...     
                                 uidx_h(maps_h.p_idx); ...    
                                 uidx_h(maps_h.T_idx)]);       
    prob = coco_add_slot(prob, 'homTst', @coco_save_data, [], 'save_full');  % Add data saving slot

    % Add events for each special point to be detected
    for i = 1:length(SP_points)
        event_name = strcat(SP_points{i}, '_');  % Create event name by appending "_"
        prob = coco_add_event(prob, event_name, SP_points{i}, 0);  % Add event to COCO
    end
    

    %%%    Define Problem Parameters       
    prob = coco_add_pars(prob, 'pars_system', uidx_h(maps_h.p_idx), ...
                         {'mu', 'eta', 'k2', 'ep'});  

    prob = coco_add_pars(prob, 'pars_T', uidx_h(maps_h.T_idx), 'T_hom', 'inactive');  % (Half) period parameter
end

function hom_data = eig_update(prob, hom_data, cseg, varargin)
    % eig_update updates hom_data with eignspace of the previous continuation step.
    %
    % Inputs:
    %   - cseg: Chart data.
    %
    % Output:
    %   - data: Updated data structure with new eigendata.

    % Set up
    [data_0, uidx] = coco_get_func_data(prob, coco_get_id(hom_data.ep_cid, 'ep'), 'data', 'uidx');
    maps_0 = data_0.ep_eqn;                
    u = cseg.src_chart.x(uidx);            
    p0 = u(maps_0.p_idx);                  
    x0 = u(maps_0.x_idx);                  
    
    % Store eignspaces 
    hom_data.eigdata = computeOrthogComplement(x0, p0, hom_data);
    hom_data.p0 = p0;                          
end

function data = computeOrthogComplement(x0, p0, hom_data)
    % computeOrthogComplement: recomputes orgonal complement eigenspace information based on the current state and parameters.
    %

    f = hom_data.f;  
    J = f.dfdx(x0, p0); 
    [eigvec, eigval] = eig(J);
    eigval = diag(eigval); 
  
    % Identify unstable and stable eigenvectors based on eigenvalues
    unstable_index = eigval > 0;  
    stable_index = eigval < 0;    

    % Extract and normalize unstable and stable eigenvectors
    v_un = eigvec(:, unstable_index);
    v_st = eigvec(:, stable_index);

    % Compute the orthogonal complements of the unstable and stable eigenvectors
    v_un_star = null(v_un')';  
    v_st_star = null(v_st')'; 

    % Normalise 
    for i = 1:size(v_st_star,1)
        v_st_star(i,:) = v_st_star(i,:) / norm(v_st_star(i,:));
    end
    for i = 1:size(v_un_star,1)
        v_un_star(i,:) = v_un_star(i,:) / norm(v_un_star(i,:));
    end

    % Eignspaces
    W_UN_star = v_un_star; 
    V_ST_star = v_st_star; 
    data = {V_ST_star, W_UN_star};
end

function data = phase_update(prob, data, cseg, varargin)
    % phase_update. Stores the phase condition before each continuation step.
    % 

    % Retrieve collocation segment data and indices
    [fdata, uidx] = coco_get_func_data(prob, coco_get_id(data.hom_cid, 'coll'), 'data', 'uidx');
    maps = fdata.coll_seg.maps;  
    u = cseg.src_chart.x;       

    % Compute the initial phase condition value
    data.xp0 = u(uidx(maps.xbp_idx))' * data.intfac;
end

function [prob, status, xtr] = hom_remesh(prob, data, chart, old_u, old_V) %#ok<INUSD>
    % Remeshing function related to the "BC_phase" condition. Updates
    % discrisiation data.
    %
    [fdata, uidx] = coco_get_func_data(prob, coco_get_id(data.hom_cid, 'coll'), 'data', 'uidx');
    xtr       = []; % No invariant indices
    data   = hom_init_data(prob, data.hom_cid, data);      % Rebuild toolbox data
    prob      = coco_change_func(prob, data, 'uidx', uidx(fdata.coll_seg.maps.xbp_idx));
    status    = 'success';
end





