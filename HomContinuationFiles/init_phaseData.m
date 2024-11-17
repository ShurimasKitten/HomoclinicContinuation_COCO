function data_hom = init_phaseData(prob, data_hom)
    % phase_init_data initializes data for discretized phase conditions in homoclinic continuation.
    %
    % Output:
    %   - data_hom: Updated data structure with phase discrisiation information
    %

    % Extract function data and solution vector 'u0' from the collocation segment
    [fdata, u0] = coco_get_func_data(prob, coco_get_id(data_hom.hom_cid, 'coll'), 'data', 'u0');
    
    % Extract mappings and mesh information from the collocation segment data
    maps  = fdata.coll_seg.maps;  % Structure containing various index mappings
    mesh = fdata.coll_seg.mesh;  % Structure containing mesh-related information
    
    % Populate initial phase condition    
    data_hom.intfac = maps.Wp' * mesh.wts2 * maps.W;  % Integration factor matrix
    data_hom.xp0 = u0(maps.xbp_idx)' * data_hom.intfac;  % Initial phase condition value
end
