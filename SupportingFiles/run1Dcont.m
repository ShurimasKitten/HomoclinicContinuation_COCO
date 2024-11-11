function run1Dcont(PS, newFn, par, bifvar, varint, SP)
    % Run a one-parameter continuation from a given initial setup.
    %
    % ARGS:
    %   - PS: Struct array of COCO settings.
    %   - newFn: Name of the output continuation files.
    %   - par: Vector of parameter values [mu, eta, ep].
    %   - bifvar: The bifurcation parameter to vary ('mu', 'eta', or other).
    %   - varint: The interval for the continuation [min, max].
    %   - SP: (Optional) Special point to monitor for termination. Defaults to false.
    %
    % This function sets up and runs a one-parameter continuation using COCO,
    % starting from an equilibrium and tracking bifurcations.

    f = PS.f;  % Dynamical system function.

    % Default SP to false if not provided.
    if nargin < 6
        SP = false;
    end

    % Scaling constants and system parameters.
    SV = 10^6;              
    sig = 2.1 * 10^4 * SV; 
    mu = par(1);
    eta = par(2);
    ep = par(3);
    
    % Maximum number of continuation attempts.
    maxCheck = 50;

    % Continuation loop.
    for i = 1:maxCheck 
        % Initialize the COCO problem structure.
        prob = coco_prob();
        
        % Find initial equilibrium for the given parameter values.
        x0 = FindEq(mu, eta, 1)';  % Assumes `FindEq` computes the equilibrium.

        % CONTINUATION SETTINGS
        prob = coco_set(prob, 'cont', 'norm', PS.contSettings.norm);          % Norm for continuation
        prob = coco_set(prob, 'cont', 'PtMX', PS.contSettings.PtMX);          % Max number of continuation steps
        prob = coco_set(prob, 'cont', 'h0', PS.contSettings.h0);              % Initial step size
        prob = coco_set(prob, 'cont', 'h_min', PS.contSettings.h_min);        % Minimum step size
        prob = coco_set(prob, 'cont', 'h_max', PS.contSettings.h_max);        % Maximum step size

        % COLLOCATION SETTINGS
        prob = coco_set(prob, 'cont', 'NTST', PS.collSettings.NTST);          % Number of mesh intervals
        prob = coco_set(prob, 'cont', 'NCOL', PS.collSettings.NCOL);          % Number of collocation points
        prob = coco_set(prob, 'ode', 'hfac1', 1e-8);                          % Step-size reduction factor
        prob = coco_set(prob, 'ode', 'hfac2', 1e-5);                          % Step-size growth factor

        % CORRECTOR SETTINGS
        prob = coco_set(prob, 'corr', 'TOL', PS.corrSettings.TOL);            % Convergence tolerance
        prob = coco_set(prob, 'corr', 'ItMX', PS.corrSettings.ItMX);          % Max number of corrector iterations

        % Initialize the continuation problem with an equilibrium.
        prob = ode_isol2ep(prob, '', f, x0(:,1), ...
                           {'mu', 'eta', 'k2', 'ep'}, [mu; eta; 8.0 * SV / sig; ep]);
    
        % Add additional data (e.g., stream function for visualization).
        prob = ep_add_bddat(prob, '', 'psi', @plotStream);

        % Run the continuation based on the specified bifurcation parameter.
        if bifvar == "mu"
            coco(prob, newFn, [], 1, {'mu', 'eta'}, varint);  % Vary mu
        elseif bifvar == "eta"
            coco(prob, newFn, [], 1, {'eta', 'mu'}, varint);  % Vary eta
        else
            coco(prob, newFn, [], 1, bifvar, varint);         % Generic bifurcation variable
        end

        % Check for special points (SP) and terminate if conditions are met.
        if SP == false || sum(ismember(coco_bd_col(newFn, 'TYPE'), SP)) > 0
            return;
        end
    end
end
