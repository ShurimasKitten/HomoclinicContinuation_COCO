function prob = proj_isol2hom(fnPO, idx_hom, continuationSettings)
    % PROJ_ISOL2HOM Initializes a COCO problem for homoclinic continuation from a periodic solution 
    % using projection boundary conditions
    %
    % ARGS:
    %   - fnPO (string): Filename of the PO continuation run data file.
    %   - idx_hom (integer): Index of the candidate homoclinic orbit from the PO continuation run.
    %   - continuationSettings (struct): Nested structure containing COCO settings, including
    %                                     continuation, correction, and collocation parameters.
    %
    % RETURNS:
    %   - prob (struct): COCO problem structure initialized for homoclinic continuation.
    %
    % Example:
    %   settings.contSettings.norm = 'some_norm';
    %   settings.contSettings.PtMX = 100;
    %   settings.corrSettings.TOL = 1e-6;
    %   settings.collSettings.NTST = 10;
    %   settings.collSettings.NCOL = 4;
    %   prob = proj_isol2hom('po_data', 5, settings);
    %
    
    % Construct COCO problem
    prob = coco_prob();
    prob = configure_COCO_settings(prob, continuationSettings);
       
    % Initilize Hom problem data
    hom_data = init_homProjData(fnPO, idx_hom, continuationSettings);
   
    % Add ODE problems
    prob = ode_isol2coll(prob, hom_data.hom_cid, hom_data.f, hom_data.t_sol, hom_data.x_sol, hom_data.p0);    
    prob = ode_isol2ep(prob, hom_data.ep_cid, hom_data.f, hom_data.x_ss, hom_data.p0);
    
    % Initialize Phase Condition data   
    hom_data = init_phaseData(prob, hom_data);

    % Close BVP
    prob = close_HomoclinicProjectionBC(prob, hom_data);
end

function prob = configure_COCO_settings(prob, PS)
    % configure_COCO_settings Applies continuation settings to the COCO problem structure
    %
    % ARGS:
    %   - prob (struct): Existing COCO problem structure to be configured.
    %   - PS (struct): Nested structure containing COCO settings, organized as:
    %                 PS.contSettings: Continuation settings
    %                 PS.corrSettings: Correction settings
    %                 PS.collSettings: Collocation settings
    %
    % RETURNS:
    %   - prob (struct): Updated COCO problem structure with applied settings.
    %
    % Example:
    %   prob = configure_COCO_settings(prob, settings);
    %
    
    %% Apply Continuation Settings
    prob = coco_set(prob, 'cont', 'norm', PS.contSettings.norm);
    prob = coco_set(prob, 'cont', 'PtMX', PS.contSettings.PtMX);
    prob = coco_set(prob, 'cont', 'h0', PS.contSettings.h0); 
    prob = coco_set(prob, 'cont', 'h_min', PS.contSettings.h_min);
    prob = coco_set(prob, 'cont', 'h_max', PS.contSettings.h_max);
    prob = coco_set(prob, 'cont', 'NAdapt', PS.contSettings.NAdapt);
    % prob = coco_set(prob, 'cont', 'TOL', PS.corrSettings.mTOL);

    %% Apply Correction Settings
    prob = coco_set(prob, 'corr', 'TOL', PS.corrSettings.TOL);
    prob = coco_set(prob, 'corr', 'ItMX', PS.corrSettings.ItMX);

    %% Apply Collocation Settings
    prob = coco_set(prob, 'coll', 'NTST', PS.collSettings.NTST);
    prob = coco_set(prob, 'coll', 'NCOL', PS.collSettings.NCOL);
    prob = coco_set(prob, 'coll', 'MXCL', 'off');
    prob = coco_set(prob, 'coll', 'bifus', 'off');
    
    %% Additional Collocation Settings. Very intense settings, which may be relaxed for simpler problems
    prob = coco_set(prob, 'coll', 'NTSTMN', 80);    
    prob = coco_set(prob, 'coll', 'NTSTMX', 800);    
    prob = coco_set(prob, 'coll', 'TOLDEC', 1e-6);
    prob = coco_set(prob, 'coll', 'TOLINC', 1e-5);
end
