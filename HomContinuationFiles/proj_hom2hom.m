function prob = proj_hom2hom(fnPO, idx_hom, continuationSettings)
    % PROJ_ISOL2HOM Initializes a COCO problem for homoclinic continuation using projection boundary conditions
    %
    % This function sets up a COCO (Continuation COntinuation) problem structure tailored for
    % homoclinic orbit continuation. It leverages data from a periodic orbit (PO) continuation run
    % to initialize the necessary conditions and parameters for projecting boundary conditions.
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
    

    %% Initialize COCO Problem and Apply Settings
    
    % Create an initial (empty) COCO problem structure
    prob = coco_prob();
    
    % Apply the continuation settings to the COCO problem using a helper function
    prob = child_applySettings(prob, continuationSettings);
   
    %% Initialize Homoclinic Data from PO Continuation
    
    % Initialize data required for homoclinic continuation using data from the PO continuation run
    hom_data = hom2hom_ProjData(fnPO, idx_hom, continuationSettings);
   
    %% Construct COCO Collocation and Equilibrium Problems
    % Add an ODE collocation problem to the COCO problem structure
    prob = ode_isol2coll(prob, 'x', hom_data.f, hom_data.t_sol, hom_data.x_sol, hom_data.p0);
    
    % Add an equilibrium point problem to the COCO problem structure
    prob = ode_isol2ep(prob, 'x_ss', hom_data.f, hom_data.x_ss, hom_data.p0);
    
    % Initialize Phase Conditions    
    hom_data = init_phaseData(prob, 'x', hom_data);

    % Set Up Boundary Conditions
    prob = close_HomoclinicProjectionBC(prob, hom_data);
end

function prob = child_applySettings(prob, PS)
    % CHILD_APPLYSETTINGS Applies continuation settings to the COCO problem structure
    %
    % This helper function configures various continuation, correction, and collocation
    % settings for the COCO problem structure based on the provided settings structure.
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
    %   prob = child_applySettings(prob, settings);
    
    %% Apply Continuation Settings
    prob = coco_set(prob, 'cont', 'norm', PS.contSettings.norm);
    prob = coco_set(prob, 'cont', 'PtMX', PS.contSettings.PtMX);
    prob = coco_set(prob, 'cont', 'h0', PS.contSettings.h0); 
    prob = coco_set(prob, 'cont', 'h_min', PS.contSettings.h_min);
    prob = coco_set(prob, 'cont', 'h_max', PS.contSettings.h_max);
    prob = coco_set(prob, 'cont', 'NAdapt', 0); %PS.contSettings.NAdapt);
    % prob = coco_set(prob, 'cont', 'TOL', PS.corrSettings.TOL);
    prob = coco_set(prob, 'cont', 'NTST', PS.collSettings.NTST);
    prob = coco_set(prob, 'cont', 'NCOL', PS.collSettings.NCOL);

    %% Apply Correction Settings
    prob = coco_set(prob, 'corr', 'TOL', PS.corrSettings.TOL);
    prob = coco_set(prob, 'corr', 'ItMX', PS.corrSettings.ItMX);
    
    %% Apply Collocation Settings
    prob = coco_set(prob, 'coll', 'NTST', PS.collSettings.NTST);
    prob = coco_set(prob, 'coll', 'NCOL', PS.collSettings.NCOL);
    prob = coco_set(prob, 'coll', 'MXCL', 'off');
    prob = coco_set(prob, 'coll', 'bifus', 'off');
    
    %% Additional Collocation Settings
    % prob = coco_set(prob, 'coll', 'NTSTMN', 50);    
    % prob = coco_set(prob, 'coll', 'NTSTMX', 300);    
    % prob = coco_set(prob, 'coll', 'TOLDEC', 1e-7);
    % prob = coco_set(prob, 'coll', 'TOLINC', 1e-6);
end
