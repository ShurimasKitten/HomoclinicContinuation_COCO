function PO_hb2po(PS, branchName, labID, newName, checkBif)
    % Starts a periodic orbit (PO) continuation from a Hopf branch/point.
    %
    % Arguments:
    %   - PS: Struct array of COCO settings.
    %   - branchName: Name of the branch containing the Hopf bifurcation (HB) point.
    %   - labID: Label of the Hopf bifurcation point.
    %   - newName: Name for the new continuation data files.
    %   - checkBif: (Optional) Check for bifurcations ('on' or 'off'). Defaults to 'on'.
    %

    % Default the `checkBif` argument to 'on' if not provided
    if nargin < 5
        checkBif = 'on';
    end
    
    % Initialize COCO problem structure
    prob = coco_prob();
    
    % Set bifurcation monitoring for periodic orbits (if enabled)
    prob = coco_set(prob, 'po', 'bifus', checkBif);

    % Continuation settings
    prob = coco_set(prob, 'cont', 'norm', PS.contSettings.norm);          % Norm for continuation
    prob = coco_set(prob, 'cont', 'PtMX', PS.contSettings.PtMX);          % Maximum continuation points
    prob = coco_set(prob, 'cont', 'h0', PS.contSettings.h0);              % Initial step size
    prob = coco_set(prob, 'cont', 'h_min', PS.contSettings.h_min);        % Minimum step size
    prob = coco_set(prob, 'cont', 'h_max', PS.contSettings.h_max);        % Maximum step size
    prob = coco_set(prob, 'cont', 'NAdapt', PS.contSettings.NAdapt);      % Step-size adaptation frequency
    prob = coco_set(prob, 'cont', 'NCOL', PS.collSettings.NCOL);          % Collocation nodes per interval
    prob = coco_set(prob, 'cont', 'NTST', PS.collSettings.NTST);          % Number of mesh intervals

    % Corrector settings
    prob = coco_set(prob, 'corr', 'TOL', PS.corrSettings.TOL);            % Tolerance for convergence
    prob = coco_set(prob, 'corr', 'ItMX', PS.corrSettings.ItMX);          % Maximum number of iterations

    % Collocation settings
    prob = coco_set(prob, 'coll', 'NCOL', PS.collSettings.NCOL);          % Collocation nodes per interval
    prob = coco_set(prob, 'coll', 'NTST', PS.collSettings.NTST);          % Number of mesh intervals

    % Additional collocation settings for finer control
    prob = coco_set(prob, 'coll', 'NTSTMN', 80);                          % Minimum number of mesh intervals
    prob = coco_set(prob, 'coll', 'NTSTMX', 300);                         % Maximum number of mesh intervals
    prob = coco_set(prob, 'coll', 'TOLDEC', 1e-8);                        % Decrease tolerance threshold
    prob = coco_set(prob, 'coll', 'TOLINC', 1e-7);                        % Increase tolerance threshold

    % Construct the periodic orbit continuation problem from Hopf bifurcation
    prob = ode_HB2po(prob, '', branchName, labID);

    % Add an event for detecting bifurcations if enabled
    if checkBif ~= 'off'
        prob = coco_add_event(prob, 'TR', 'po.test.NS', 0);  % Event for Neimark-Sacker bifurcation
    end
    
    % Run continuation with specified parameters and save results
    coco(prob, newName, [], 1, ...
        [PS.contSettings.PO.contVar{1}, 'po.test.USTAB', 'po.test.NS'], ...
        PS.contSettings.PO.contVar{2});
end
