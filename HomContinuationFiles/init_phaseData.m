function data_hom = init_phaseData(prob, tbid, data_hom)
    % phase_init_data initializes data for discretized phase conditions in homoclinic continuation.
    %
    % Inputs:
    %   - prob: COCO continuation problem structure.
    %   - tbid: Toolbox identifier string (e.g., 'x.coll') used to retrieve specific data.
    %   - data_hom: Data structure to be populated with phase condition initialization data.
    %
    % Output:
    %   - data_hom: Updated data structure containing initialization data for phase conditions.
    %

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%    Retrieve Function and Solution Data    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Generate the unique identifier for the collocation toolbox segment
    stbid = coco_get_id(tbid, 'coll');  % Example: 'x.coll'
    
    % Extract function data and solution vector 'u0' from the collocation segment
    [fdata, u0] = coco_get_func_data(prob, stbid, 'data', 'u0');
    
    % Extract mappings and mesh information from the collocation segment data
    maps  = fdata.coll_seg.maps;  % Structure containing various index mappings
    mesh = fdata.coll_seg.mesh;  % Structure containing mesh-related information
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%    Populate Initialization Data    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Store index mappings for the equilibrium points in the homoclinic orbit
    data_hom.x0_idx = maps.x0_idx;  % Index of x0 (starting point) in the state vector
    data_hom.x1_idx = maps.x1_idx;  % Index of x1 (ending point) in the state vector
    
    % Compute the integration factor used in the discretized phase condition
    % This factor is derived from the collocation mappings and mesh weights
    data_hom.intfac = maps.Wp' * mesh.wts2 * maps.W;  % Integration factor matrix
    
    % Calculate the initial phase condition value by projecting the state vector
    % at the phase condition point (xbp) using the integration factor
    data_hom.xp0 = u0(maps.xbp_idx)' * data_hom.intfac;  % Initial phase condition value
end
